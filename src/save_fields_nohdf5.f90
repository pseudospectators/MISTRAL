!-------------------------------------------------------------------------------
! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files. 
!-------------------------------------------------------------------------------
subroutine save_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use solid_model
  use insect_module
  implicit none

  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid),dimension(1:nBeams),intent(inout) :: beams
  type(diptera),intent(inout) :: Insect
  
  real(kind=pr) :: tmp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr) :: volume, t1
  character(len=6) :: name

  t1=MPI_wtime()
  
  !--Set up file name base    
  if ( save_only_one_period == "yes" ) then
    ! overwrite files from last period to save disk space
    ! i.e. t=1.05 is written to t=0.05, as well as 2.05 and 3.05
    write(name,'(i6.6)') floor( (time%time-real(floor(time%time/tsave_period)))*1000.d0 )
  else
    ! name is just the time
    write(name,'(i6.6)') floor(time%time*1000.d0) 
  endif
  
  if (mpirank == 0 ) then
    write(*,'("Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A," ...")',advance='no') & 
    time%time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
  endif

  ! Save the velocity
  if (isaveVelocity == 1) then
    call save_field_nohdf5(time,"./ux_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"ux")
    call save_field_nohdf5(time,"./uy_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"uy")
    call save_field_nohdf5(time,"./uz_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"uz")
  endif  
  
  ! Save the pressure
  if (isavePress == 1) then  
    call save_field_nohdf5(time,'./p_'//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4),"p")
  endif
     
  ! Save the Vorticity
  if (isaveVorticity==1) then
    !-- compute vorticity:
    call curl_x( u(:,:,:,1:3), tmp )
    call save_field_nohdf5(time,"./vorx_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"vorx")
    call save_field_nohdf5(time,"./vory_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"vory")
    call save_field_nohdf5(time,"./vorz_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"vorz")
  endif
      
  ! Save the Mask
  if (isaveMask == 1 .and. iPenalization == 1) then
    ! create mask at current time
    if(iMoving==1) call create_mask( time%time, mask, mask_color, us, Insect, beams, 0 )
    ! make sure mask is between 0 and 1 
    mask = mask*eps
    ! warn if we're about to store an empty mask
    call compute_mask_volume(mask,volume)
    if ((mpirank==0).and.(volume<1e-10)) write(*,*) "WARNING: saving empty mask"
    ! save the mask
    call save_field_nohdf5(time,'./mask_'//name,mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),"mask")
    mask = mask/eps
  endif
  
  ! save solid velocity
  if (isaveSolidVelocity == 1 .and. iPenalization == 1 .and. iMoving == 1) then
    call save_field_nohdf5(time,'./usx_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"usx")
    call save_field_nohdf5(time,'./usy_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"usy")
    call save_field_nohdf5(time,'./usz_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"usz")
  endif
  
  time_save = time_save + MPI_wtime() - t1
  if (mpirank==0) write(*,*) " ...DONE!"
end subroutine save_fields


! Write the field field_out to file filename.
subroutine save_field_nohdf5(time,filename,field_out,dsetname)
  use mpi
  use vars
  implicit none
  
  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  type(timetype), intent (in) :: time
  character(len=*), intent (in) :: filename
  character(len=*), intent (in) :: dsetname

  integer :: error  ! error flags

  ! HDF attribute variables
  integer, parameter :: arank = 1

  integer :: mpierror, i
  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  t1 = MPI_wtime()
  
  !!! Tell HDF5 how our  data is organized:
!  dimensions_file = (/nx,ny,nz/)
!  offset(1) = ra(1)
!  offset(2) = ra(2)
!  offset(3) = ra(3)
!  dimensions_local(1) = rb(1)-ra(1) +1
!  dimensions_local(2) = rb(2)-ra(2) +1
!  dimensions_local(3) = rb(3)-ra(3) +1
  ! each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
!  do i = 1, 3
!     call MPI_REDUCE(dimensions_local(i),chunking_dims(i),1, &
!          MPI_INTEGER8,MPI_MAX,0,&
!          MPI_COMM_WORLD,mpierror)
!     call MPI_BCAST(chunking_dims(i),1,MPI_INTEGER8,0, &
!          MPI_COMM_WORLD,mpierror )
!  enddo

  !!! Set up the HDF data structures:

  ! Initialize HDF5 library and Fortran interfaces.
!  call h5open_f(error)
  ! Setup file access property list with parallel I/O access.
  ! this sets up a property list (plist_id) with standard values for
  ! FILE_ACCESS
!  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! Modify the property list and store the MPI IO comminucator
  ! information in the file access property list
!  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file collectively. (existing files are overwritten)
!#ifdef TURING  
!  if ( index(filename,'.h5')==0 ) then
!    ! file does not contain *.h5 ending -> add suffix
!    call h5fcreate_f('bglockless:'//trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
!    file_id, error, access_prp = plist_id)
!  else
    ! field DOES contain .h5 ending -> just write
!    call h5fcreate_f('bglockless:'//trim(adjustl(filename)), H5F_ACC_TRUNC_F, &
!    file_id, error, access_prp = plist_id)
!  endif
!#else
!  if ( index(filename,'.h5')==0 ) then
!    ! file does not contain *.h5 ending -> add suffix
!    call h5fcreate_f(trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
!    file_id, error, access_prp = plist_id)
!  else
    ! field DOES contain .h5 ending -> just write
!    call h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, &
!    file_id, error, access_prp = plist_id)
!  endif
!#endif  


  ! this closes the property list plist_id (we'll re-use it)
!  call h5pclose_f(plist_id, error)

  ! Create the data space for the  dataset.
  ! Dataspace in the file: contains all data from all procs
!  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
!  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset.
  ! NB: chunking and hyperslab are unrelated
!  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
!  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)
  ! Output files are single-precition
!  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, &
!       dset_id, error, plist_id)
!  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
!  call h5dget_space_f(dset_id, filespace, error)
!  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
!       error, stride, dimensions_local)

  ! Create property list for collective dataset write
!  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively, double precision in memory
!  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, &
!       error, file_space_id = filespace, mem_space_id = memspace,&
!       xfer_prp = plist_id)

  !!! Write the attributes to the HDF files.
  ! The attributes written are time, penalisation parameter,
  ! computational resolution, and physical domain size.
!  adims = (/1/)
!  call write_attribute_dble(adims,"time",(/time%time/),1,dset_id)
!  call write_attribute_dble(adims,"epsi",(/eps/),1,dset_id)
!  adims = (/3/)
!  call write_attribute_dble(adims,"domain_size",(/xl,yl,zl/),3,dset_id)
!  call write_attribute_int(adims,"nxyz",(/nx,ny,nz/),3,dset_id)

  !!! Close dataspaces:
!  call h5sclose_f(filespace, error)
!  call h5sclose_f(memspace, error)
!  call h5dclose_f(dset_id, error) ! Close the dataset.
!  call h5pclose_f(plist_id, error) ! Close the property list.
!  call h5fclose_f(file_id, error) ! Close the file.
!  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.

  ! write the XMF data for all of the saved fields
!  if ((mpirank==0).and.(isaveXMF==1)) then
!     ! the filename contains a leading "./" which we must remove
!     call Write_XMF(time%time,&
!          trim(adjustl(filename(3:len(filename)))),&
!          trim(adjustl(dsetname))&
!          )
!  endif

  time_hdf5=time_hdf5 + MPI_wtime() - t1 ! performance analysis
end subroutine save_field_nohdf5



! write a runtime backup to disk. We're saving the fluid state vector 
! u=(/ux,uy,uz,p/) in one HDF5 file and some FSI/insect information in separate 
! files. The hdf5 file is written in double precision, to allow RK2/RK4 schemes
! to resume exactly (without loss of precision) where it left. AB2 schemes are
! different - they apply the startup step again (this associated with loos of 
! precision)
subroutine dump_runtime_backup(time,nbackup,u,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none
  
  type(timetype),intent(inout) :: time
  integer,intent(inout) :: nbackup
  real(kind=pr),intent(in)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  type(solid),dimension(1:nBeams),intent(in) :: beams
  type(diptera),intent(in) :: Insect

  character(len=18) :: filename
  real(kind=pr) :: tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: t1
  integer :: error  ! error flags

  t1=MPI_wtime() ! performance diagnostic

  ! Create name for the backup file. We keep at any time at most 2 sets of 
  ! backups, runtime_backupX.h5 with X=0,1  
  if(mpirank == 0) then
     write(*,'("Dumping runtime_backup",i1,".h5 (time=",es12.4,") to disk....")',&
     advance='no') nbackup, time%time
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1,".h5")') nbackup

  ! Initialize HDF5 library and Fortran interfaces:
!  call h5open_f(error)

  !!! Setup file access property list with parallel I/O access.
  ! Set up a property list ("plist_id") with standard values for
  ! FILE_ACCESS:
!  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! Modify the property list and store MPI IO comminucator information
  ! in the file access property list:
!  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file collectively. (existing files are overwritten)
!  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, &
!       access_prp = plist_id)
  ! Close the property list (we'll re-use it)
!  call h5pclose_f(plist_id, error)

  ! Write the fluid backup field:
  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1)
  call dump_field_backup(tmp,"ux",time%time,time%dt_old,time%dt_new,time%n1,time%it)
  
  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2)
  call dump_field_backup(tmp,"uy",time%time,time%dt_old,time%dt_new,time%n1,time%it)
  
  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3)
  call dump_field_backup(tmp,"uz",time%time,time%dt_old,time%dt_new,time%n1,time%it)
  
  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4)
  call dump_field_backup(tmp,"p",time%time,time%dt_old,time%dt_new,time%n1,time%it)

       
  ! Close the file:
!  call h5fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library:
!  call h5close_f(error)
  
  
  !-------------------------------------------------------------------------
  ! backup for the rigid body solver (free-flight insect)
  !-------------------------------------------------------------------------
  if ((method=="fsi").and.(mpirank==0)) then
  if (iMask=="Insect") then
  if ((Insect%BodyMotion=="takeoff").and.(Insect%KineFromFile=="simplified_dynamic")) then
    write (*,'(A)',advance="no") "insect bckp in "//filename//".rigidsolver"
    open(10, file=filename//".rigidsolver", form='formatted', status='replace') 
    write(10, *) SolidDyn%var_new, SolidDyn%var_this,&
                 SolidDyn%rhs_this, SolidDyn%rhs_old
    close(10)
  endif
  endif
  endif
  
  !-------------------------------------------------------------------------
  !-- backup solid solver, if doing active FSI
  !-------------------------------------------------------------------------
  if((use_solid_model=="yes").and.(method=="fsi")) then
    call dump_solid_backup( time%time, beams, nbackup )
  endif    
  
  nbackup = 1 - nbackup
  time_bckp=time_bckp + MPI_wtime() -t1 ! Performance diagnostic

  if(mpirank == 0) write(*,'(A)') "...DONE!"
end subroutine dump_runtime_backup


! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "file_id". Attributes are stores in one attribute
! "bckp" which contains 8 values
subroutine dump_field_backup(field,dsetname,time,dt0,dt1,n1,it)
  use mpi
  use vars
  implicit none

  real(kind=pr),intent(inout) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time,dt1,dt0
  character(len=*), intent (in) :: dsetname
  integer,intent(in) :: n1,it

  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  integer :: error, mpierror, i  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes

!  dimensions_file = (/nx,ny,nz/)
!  dimensions_local(1) = rb(1)-ra(1) +1
!  dimensions_local(2) = rb(2)-ra(2) +1
!  dimensions_local(3) = rb(3)-ra(3) +1

  ! Offsets
!  offset(1) = ra(1)
!  offset(2) = ra(2)
!  offset(3) = ra(3)

  ! Each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
!  do i = 1, 3
!     call MPI_REDUCE(dimensions_local(i), chunking_dims(i),1, &
!          MPI_INTEGER8, MPI_MAX,0,MPI_COMM_WORLD,mpierror)
!     call MPI_BCAST(chunking_dims(i), 1, MPI_INTEGER8, 0, &
!          MPI_COMM_WORLD, mpierror)
!  enddo

  ! -----------------------------------
  ! Create the data space for the  dataset.
  ! -----------------------------------
  ! Dataspace in the file: contains all data from all procs
!  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
!  call h5screate_simple_f(rank, dimensions_local, memspace, error)
  ! Create chunked dataset.
!  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
!  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)
!  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
!       dset_id, error, plist_id)
!  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
!  call h5dget_space_f(dset_id, filespace, error)
!  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
!       error , stride, dimensions_local)

  ! Create property list for collective dataset write
!  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively.
!  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dimensions_file, error, &
!       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

  ! ------
  ! Attributes (we save everything in one, all double. to be converted
  ! when reading (to integer)
  ! ------
!  adims = (/8/)
!  allocate (attributes(1:8))
!  aname = "bckp"
!  attributes = (/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/)
!  call write_attribute_dble(adims,aname,attributes,8,dset_id)

  ! Close dataspaces, dataset and property list
!  call h5sclose_f(filespace, error)
!  call h5sclose_f(memspace, error)
!  call h5dclose_f(dset_id, error)
!  call h5pclose_f(plist_id, error)

!  deallocate(attributes)
end subroutine dump_field_backup





! Read in a single file that follows the naming convention
! this is a serial routine (parallel version below)
! note you need to know what dimension the file has,
! call fetch_attributes first
subroutine read_single_file_serial(filename,field)
  use mpi
  use vars
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr), intent (out) :: field(0:nx-1,0:ny-1,0:nz-1)
  
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  character(len=80)             :: dsetname

  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  
  ! the dataset is named the same way as the file: (this is convention)
  dsetname = filename ( 1:index( filename, '_' )-1 )

  if (mpisize>1) then
    write (*,*) "this routine is currently serial only"
    call abort()
  endif
  
  ! check if file exist
  call check_file_exists( filename )

  ! Initialize HDF5 library and Fortran interfaces.
!  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
!  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
!  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
!  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
!  call h5pclose_f(plist_id, error)
  
  ! Definition of memory distribution
!  dimensions_file = (/nx,ny,nz/)
!  dimensions_local(1) = nx
!  dimensions_local(2) = ny
!  dimensions_local(3) = nz

!  offset(1) = 0
!  offset(2) = 0
!  offset(3) = 0
  
!  chunking_dims(1) = nx
!  chunking_dims(2) = ny
!  chunking_dims(3) = nz

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
!  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
!  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
!  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
!  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
!  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
!  call h5dget_space_f(dset_id, filespace, error)
!  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
!       error , stride, dimensions_local)

  ! Create property list for collective dataset read
!  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

!  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
!       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

!  call h5sclose_f(filespace, error)
!  call h5sclose_f(memspace, error)
!  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
!  call h5dclose_f(dset_id, error)
!  call h5fclose_f(file_id,error)
!  call H5close_f(error)
  
  
end subroutine read_single_file_serial




! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
subroutine Read_Single_File ( filename, field )
  use mpi
  use vars
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),&
  dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
  intent (out) :: field
  
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time, xl_file, yl_file, zl_file
  character(len=80)             :: dsetname
  integer                       :: nx_file, ny_file, nz_file, mpierror, i  

  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1

  ! the dataset is named the same way as the file: (this is convention)
  dsetname = filename ( 1:index( filename, '_' )-1 )
  if (mpirank==0) then
    write (*,'("Reading file ",A,"  .....")',advance='no') trim(adjustl(filename))
  endif
  
  !-----------------------------------------------------------------------------
  ! perform tests
  !-----------------------------------------------------------------------------
  call check_file_exists ( filename )
  
  ! fetch attributes from file to see if it is a good idea to load it
!  call Fetch_attributes( filename, dsetname,nx_file,ny_file,nz_file,& 
!                         xl_file,yl_file ,zl_file,time )
         
  ! if the resolutions do not match, yell and hang yourself       
  if ((nx.ne.nx_file).or.(ny.ne.ny_file).or.(nz.ne.nz_file)) then                         
    if (mpirank == 0) then
    write (*,'(A)') "read_single_file: ERROR " // trim(filename)
    write (*,'("nx=",i4,"ny=",i4,"nz=",i4)') nx,ny,nz
    write (*,'("but in file: nx=",i4,"ny=",i4,"nz=",i4)') nx_file,ny_file,nz_file
    call abort()
    endif
  endif
  
  ! if the domain size doesn't match, proceed, but yell.
  if ((xl.ne.xl_file).or.(yl.ne.yl_file).or.(zl.ne.zl_file)) then                         
    if (mpirank == 0) then
    write (*,'(A)') "read_single_file: WARNING " // trim(filename)
    write (*,'("xl=",es12.4,"yl=",es12.4,"zl=",es12.4)')&
      xl,yl,zl
    write (*,'("but in file: xl=",es12.4,"yl=",es12.4,"zl=",es12.4)') & 
      xl_file,yl_file,zl_file
    write (*,'(A)') "proceed, with fingers crossed."
    endif
  endif
  
  !-----------------------------------------------------------------------------
  ! load the file
  !-----------------------------------------------------------------------------  
  ! Initialize HDF5 library and Fortran interfaces.
!  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
!  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
!  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
!  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
!  call h5pclose_f(plist_id, error)
  
  ! Definition of memory distribution
!  dimensions_file = (/nx,ny,nz/)
!  dimensions_local(1) = rb(1)-ra(1) +1
!  dimensions_local(2) = rb(2)-ra(2) +1
!  dimensions_local(3) = rb(3)-ra(3) +1

!  offset(1) = ra(1)
!  offset(2) = ra(2)
!  offset(3) = ra(3)
  
  ! Each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
!  do i = 1, 3
!     call MPI_REDUCE ( dimensions_local(i), chunking_dims(i),1, &
!          MPI_INTEGER8, MPI_MAX,0,MPI_COMM_WORLD,mpierror)
!     call MPI_BCAST  ( chunking_dims(i), 1, MPI_INTEGER8, 0, &
!          MPI_COMM_WORLD, mpierror )
!  enddo

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
!  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
!  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
!  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
!  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
!  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
!  call h5dget_space_f(dset_id, filespace, error)
!  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
!       error , stride, dimensions_local)

  ! Create property list for collective dataset read
!  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

!  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
!       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

!  call h5sclose_f(filespace, error)
!  call h5sclose_f(memspace, error)
!  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
!  call h5dclose_f(dset_id, error)
!  call h5fclose_f(file_id,error)
!  call H5close_f(error)
  
  if (mpirank==0) then
    write (*,'("...DONE! ")',advance='yes')
  endif
  
end subroutine Read_Single_File




! Load backup data from disk to initialize run for restart
subroutine read_runtime_backup(filename,time,u,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  character(len=*),intent(in) :: filename
  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  type(solid),dimension(1:nBeams),intent(in) :: beams
  type(diptera),intent(in) :: Insect

  integer :: error  ! Error flag

  if(mpirank == 0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! I'm trying to resume a backup file: "//filename
  endif
  
  call check_file_exists ( filename )

  ! Initialize HDF5 library and Fortran interfaces.
!  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
!  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
!  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
!  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
!  call h5pclose_f(plist_id, error)

  ! Read fluid backup field:
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"ux",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it)
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"uy",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it)
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"uz",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it)
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4),"p",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it)       

  ! Close the file:
!  call h5fclose_f(file_id,error)
  ! Close FORTRAN interfaces and HDF5 library:
!  call H5close_f(error)

  if(mpirank == 0) then
     write(*,'("time=",es15.8," dt0=",es15.8)') time%time, time%dt_old
     write(*,'("!!! DONE READING BACKUP (thats good news!)")')
     write(*,'("---------")')
  endif

end subroutine read_runtime_backup


! This routine reads a single field "dsetname" from a backup file
! "file_id". the field has the attribute "attributes", which is an 8x1
! array containing scalar backup information
subroutine read_field_backup(field,dsetname,time,dt0,dt1,n1,it)
  use mpi
  use vars
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), &
       intent(out) :: field
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (out)  :: time,dt1,dt0
  character(len=*), intent (in) :: dsetname
  integer,intent(out)           :: n1,it
  integer                       :: nx_file,ny_file,nz_file, mpierror, i

  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes

  ! Definition of memory distribution
!  dimensions_file = (/nx,ny,nz/)
!  dimensions_local(1) = rb(1)-ra(1) +1
!  dimensions_local(2) = rb(2)-ra(2) +1
!  dimensions_local(3) = rb(3)-ra(3) +1

!  offset(1) = ra(1)
!  offset(2) = ra(2)
!  offset(3) = ra(3)

  ! Each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
!  do i = 1, 3
!     call MPI_REDUCE(dimensions_local(i),chunking_dims(i),1,MPI_INTEGER8,&
!          MPI_MAX,0,MPI_COMM_WORLD,mpierror)
!     call MPI_BCAST(chunking_dims(i),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierror)
!  enddo

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
!  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
!  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
!  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
!  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
!  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
!  call h5dget_space_f(dset_id, filespace, error)
!  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
!       error , stride, dimensions_local)

  ! Create property list for collective dataset read
!  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

!  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
!       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

!  call h5sclose_f(filespace, error)
!  call h5sclose_f(memspace, error)
!  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! attributes (we save everything in one, all double. to be converted
  ! when reading (to integer)
!  adims = (/8/)
!  allocate (attributes(1:8))
!  aname = "bckp"
!  call h5aopen_f(dset_id, aname, attr_id, error)

  ! Get dataspace and read
!  call h5aget_space_f(attr_id, aspace_id, error)
!  call h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attributes, adims, error)
!  call h5aclose_f(attr_id, error) ! Close the attribute.
!  call h5sclose_f(aspace_id, error) ! Terminate access to the data space.

!  time    = attributes(1)
!  dt1     = attributes(2)
!  dt0     = attributes(3)
!  n1      = int(attributes(4))
!  it      = int(attributes(5))
!  nx_file = int(attributes(6))
!  ny_file = int(attributes(7))
!  nz_file = int(attributes(8))

!  if ( (nx_file.ne.nx).or.(nx_file.ne.nx).or.(nx_file.ne.nx)) then
!     write (*,'(A)') "!!! Thats odd...the backup you're trying to resume doesn't have the same nx,ny,nz"
!     write (*,'(A)') "I'll leave you crying and commit suicide here."
!     call abort()
!  endif

!  deallocate (attributes)

  ! Close dataset
!  call h5dclose_f(dset_id, error)
end subroutine read_field_backup


! checks if a given file ("fname") exists. if not, code is stopped brutally
subroutine check_file_exists(fname)
  use vars
  use mpi
  implicit none
  
  character (len=*), intent(in) :: fname
  logical :: exist1
  integer :: mpicode

  if (mpirank == 0) then
    inquire ( file=fname, exist=exist1 )
    if ( exist1 .eqv. .false.) then
      write (*,'("ERROR! file: ",A," not found")') trim(fname) 
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    endif  
  endif
  
end subroutine check_file_exists


! overwrite and initialize file
subroutine init_empty_file( fname )
  use vars
  implicit none
  character (len=*), intent(in) :: fname
  
  open (15, file=fname,status='replace')
  close(15)
end subroutine


!-------------------------------------------------------------------------------
! This is a stub. Postprocessing tools are not available without HDF5 support.
!-------------------------------------------------------------------------------
subroutine postprocessing()
end subroutine
