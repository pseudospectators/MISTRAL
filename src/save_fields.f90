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

  real(kind=pr), allocatable, dimension(:,:,:,:) :: tmp
  real(kind=pr) :: volume, t1
  character(len=6) :: name

  allocate(tmp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3))

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
    write(*,'(80("~"))')
    write(*,'("Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A," ...")') &
    time%time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
  endif

  ! Save the velocity
  if (isaveVelocity == 1) then
    call save_field_hdf5(time%time,"ux_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1))
    call save_field_hdf5(time%time,"uy_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2))
    call save_field_hdf5(time%time,"uz_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3))
  endif

  ! Save the pressure
  if (isavePress == 1) then
    call save_field_hdf5(time%time,'p_'//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4))
  endif

  ! Save the Vorticity
  if (isaveVorticity==1) then
    !-- compute vorticity:
    call curl_x( u(:,:,:,1:3), tmp )
    call save_field_hdf5(time%time,"vorx_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1))
    call save_field_hdf5(time%time,"vory_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2))
    call save_field_hdf5(time%time,"vorz_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3))
  endif

  if (iSaveDivergence==1) then
    call divergence( u(:,:,:,1:3), tmp(:,:,:,1) )
    call save_field_hdf5(time%time,"divu_"//name,tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1))
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
    call save_field_hdf5(time%time,'mask_'//name,mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
    mask = mask/eps
  endif

  ! save solid velocity
  if (isaveSolidVelocity == 1 .and. iPenalization == 1 .and. iMoving == 1) then
    call save_field_hdf5(time%time,'usx_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1))
    call save_field_hdf5(time%time,'usy_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2))
    call save_field_hdf5(time%time,'usz_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3))
  endif

  deallocate(tmp)
  time_save = time_save + MPI_wtime() - t1
  if (mpirank==0) then
    write(*,'(80("~"))')
  endif
end subroutine save_fields


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
  real(kind=pr),dimension(:,:,:),allocatable :: tmp
  real(kind=pr) :: t1

  t1=MPI_wtime() ! performance diagnostic

  allocate (tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  ! Create name for the backup file. We keep at any time at most 2 sets of
  ! backups, runtime_backupX.h5 with X=0,1
  if(mpirank == 0) then
     write(*,'("Dumping runtime_backup",i1,".h5 (time=",es12.4,") to disk....")',&
     advance='no') nbackup, time%time
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1,".h5")') nbackup
  call init_empty_file(filename)


  ! Write the fluid backup field:
  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1)
  call dump_field_backup(filename,tmp,"ux",time%time,time%dt_old,time%dt_new,time%n1,time%it)

  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2)
  call dump_field_backup(filename,tmp,"uy",time%time,time%dt_old,time%dt_new,time%n1,time%it)

  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3)
  call dump_field_backup(filename,tmp,"uz",time%time,time%dt_old,time%dt_new,time%n1,time%it)

  tmp(:,:,:) = u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4)
  call dump_field_backup(filename,tmp,"p",time%time,time%dt_old,time%dt_new,time%n1,time%it)

  !-------------------------------------------------------------------------
  ! backup for the rigid body solver (free-flight insect)
  !-------------------------------------------------------------------------
  ! if ((method=="fsi").and.(mpirank==0)) then
  ! if (iMask=="Insect") then
  ! if ((Insect%BodyMotion=="takeoff").and.(Insect%KineFromFile=="simplified_dynamic")) then
  !   write (*,'(A)',advance="no") "insect bckp in "//filename//".rigidsolver"
  !   open(10, file=filename//".rigidsolver", form='formatted', status='replace')
  !   write(10, *) SolidDyn%var_new, SolidDyn%var_this,&
  !                SolidDyn%rhs_this, SolidDyn%rhs_old
  !   close(10)
  ! endif
  ! endif
  ! endif

  !-------------------------------------------------------------------------
  !-- backup solid solver, if doing active FSI
  !-------------------------------------------------------------------------
  if((use_solid_model=="yes").and.(method=="fsi")) then
    call dump_solid_backup( time%time, beams, nbackup )
  endif

  nbackup = 1 - nbackup
  time_bckp=time_bckp + MPI_wtime() -t1 ! Performance diagnostic

  deallocate(tmp)

  if(mpirank == 0) write(*,'(A)') "...DONE!"
end subroutine dump_runtime_backup


!-------------------------------------------------------------------------------
! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "filename". Attributes are stores in one attribute
! "bckp" which contains 8 values
!-------------------------------------------------------------------------------
subroutine dump_field_backup(filename,field,dsetname,time,dt0,dt1,n1,it)
  use vars
  use hdf5_wrapper
  implicit none

  integer,intent(in) :: n1,it
  real(kind=pr), intent (in) :: time,dt1,dt0
  real(kind=pr),intent(in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*), intent (in) :: dsetname, filename

  call write_field_hdf5( filename,dsetname, ra, rb, field, overwrite=.false.)
  call write_attribute( filename,dsetname,"bckp",(/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/) )
end subroutine dump_field_backup


!-------------------------------------------------------------------------------
! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
!-------------------------------------------------------------------------------
subroutine Read_Single_File ( filename, field )
  use vars
  use hdf5_wrapper
  use basic_operators, only : fieldmax, fieldmin, fieldmean, checknan
  use helpers, only : get_dsetname
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(out) :: field
  integer, dimension(1:3) :: nxyz
  real(kind=pr), dimension(1:3) :: domain
  real(kind=pr), dimension(1) :: ttime, viscosity_dummy
  real(kind=pr) :: fmax,fmin,favg,t1,mbyte,t2

  t1 = MPI_wtime()

  call check_file_exists ( filename )
  call read_attribute( filename,get_dsetname(filename),"nxyz",nxyz)
  call read_attribute( filename,get_dsetname(filename),"domain_size",domain)
  call read_attribute( filename,get_dsetname(filename),"time",ttime)
  call read_attribute( filename,get_dsetname(filename),"viscosity",viscosity_dummy)

  if (mpirank==0) then
    write(*,'(40("~"))')
    write(*,'("Reading from file ",A)') trim(adjustl(filename))
    write(*,'("dsetname=",A)') trim(adjustl(get_dsetname(filename)))
    write(*,'("nx=",i4," ny=",i4," nz=",i4," time=",g12.4," viscosity=",g16.4)') nxyz,ttime(1),viscosity_dummy(1)
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') domain

    ! if the domain size doesn't match, proceed, but yell.
    if ((xl.ne.domain(1)).or.(yl.ne.domain(2)).or.(zl.ne.domain(3))) then
        write (*,'(A)') " WARNING! Domain size mismatch."
        write (*,'("in memory:   xl=",es12.4,"yl=",es12.4,"zl=",es12.4)') xl,yl,zl
        write (*,'("but in file: xl=",es12.4,"yl=",es12.4,"zl=",es12.4)') domain
        write (*,'(A)') "proceed, with fingers crossed."
    endif

    ! if the resolutions do not match, yell and hang yourself
    if ((nx/=nxyz(1)).or.(ny/=nxyz(2)).or.(nz/=nxyz(3))) then
      write (*,'(A)') "ERROR! Resolution mismatch"
      write (*,'(A)') "This happens if ra(:) and rb(:) are not properly initialized."
      write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
      write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nxyz
      call abort(125)
    endif
  endif

  ! actual reading of file
  call read_field_hdf5 ( filename, get_dsetname(filename), ra,rb, field)
  ! check if field contains NaN
  ! call checknan(field,"recently loaded field")

  fmax = fieldmax(field)
  fmin = fieldmin(field)
  favg = fieldmean(field)
  mbyte = dble(nx)*dble(ny)*dble(nz)*4.d0/1.0d+6
  t2 = MPI_wtime() - t1

  if (mpirank==0) then
    write (*,'("read ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') mbyte,t2,mbyte/t2
    write (*,'("max=",g12.4," min=",g12.4," mean=",g12.4)') fmax,fmin,favg
    write (*,'("Done reading file, Elapsed time=",g12.4,"s")') MPI_wtime() - t1
    write(*,'(40("~"))')
  endif


end subroutine Read_Single_File




! Load backup data from disk to initialize run for restart
subroutine read_runtime_backup(filename,time,u,Insect,beams)
  use vars
  use hdf5_wrapper
  use solid_model
  use insect_module
  implicit none

  character(len=*),intent(in) :: filename
  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  type(solid),dimension(1:nBeams),intent(in) :: beams
  type(diptera),intent(in) :: Insect

  real(kind=pr), dimension(1:8) :: attributes
  integer :: nx_file,ny_file,nz_file

  if(mpirank == 0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! I'm trying to resume a backup file: "//filename
  endif

  call check_file_exists ( filename )

  ! read the attribute
  call read_attribute( filename, "ux", "bckp", attributes )
  ! and extract the values
  time%time    = attributes(1)
  time%dt_new  = attributes(2)
  time%dt_old  = attributes(3)
  time%n1      = int(attributes(4))
  time%it      = int(attributes(5))
  nx_file = int(attributes(6))
  ny_file = int(attributes(7))
  nz_file = int(attributes(8))

  if ((nx/=nx_file).or.(ny/=ny_file).or.(nz/=nz_file)) then
    write (*,'(A)') "ERROR! Resolution mismatch"
    write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
    write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nx_file,ny_file,nz_file
    call abort(77776)
  endif

  call read_field_backup( filename,"ux",u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1))
  call read_field_backup( filename,"uy",u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2))
  call read_field_backup( filename,"uz",u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3))
  call read_field_backup( filename,"p",u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4))

  if(mpirank == 0) then
     write(*,'("time=",es15.8," dt0=",es15.8)') time%time, time%dt_old
     write(*,'("!!! DONE READING BACKUP (thats good news!)")')
     write(*,'("---------")')
  endif

end subroutine read_runtime_backup


!-------------------------------------------------------------------------------
! This routine reads a single field "dsetname" from a backup file
! "file_id". the field has the attribute "attributes", which is an 8x1
! array containing scalar backup information
!-------------------------------------------------------------------------------
subroutine read_field_backup(filename,dsetname,field)
  use vars
  use hdf5_wrapper
  use basic_operators, only : checknan
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(out) :: field
  character(len=*), intent (in) :: dsetname, filename

  if (mpirank==0) then
    write(*,'("Reading ",A," from backup file ",A)') trim(adjustl(dsetname)),trim(adjustl(filename))
  endif

  call read_field_hdf5( filename, dsetname, ra, rb, field )
  call checknan(field,'recently read backup file!!')

end subroutine read_field_backup


!----------------------------------------------------
! This routine fetches the resolution, the domain size and the time
! form a *.h5 file
!----------------------------------------------------
! filename: a *.h5 file to read from.
! dsetname: a dataset inside the file. in our system, we only have one per file
!           and this matches the prefix: mask_00010.h5  --> dsetname = "mask"
! note:
!           the file must contain the dataset
!           but especially the attributes "nxyz", "time", "domain_size"
!----------------------------------------------------
subroutine Fetch_attributes( filename, nx, ny, nz, xl, yl ,zl, time, viscosity )
  use hdf5_wrapper
  use helpers, only : get_dsetname
  use mpi
  implicit none

  character(len=*), intent(in) :: filename  ! file name
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time, viscosity

  real(kind=pr),dimension(1) :: attr_data1, attr_data0
  real(kind=pr),dimension(1:3) :: attr_data2
  integer,dimension(1:3) :: attr_data3

  call check_file_exists ( filename )
  call read_attribute( filename, get_dsetname(filename), "time", attr_data0)
  call read_attribute( filename, get_dsetname(filename), "viscosity", attr_data1)
  call read_attribute( filename, get_dsetname(filename), "domain_size", attr_data2)
  call read_attribute( filename, get_dsetname(filename), "nxyz", attr_data3)

  time = attr_data0(1)
  viscosity = attr_data1(1)
  xl = attr_data2(1)
  yl = attr_data2(2)
  zl = attr_data2(3)
  nx = attr_data3(1)
  ny = attr_data3(2)
  nz = attr_data3(3)
end subroutine Fetch_attributes


!-------------------------------------------------------------------------------
! Standart wrapper for the HDF5 library, saves a single array (e.g. one vector component)
! to a single HDF5 file. A bunch of useful attributes (resolution, domain size,
! penalization, viscosity, time) are stored as well.
!-------------------------------------------------------------------------------
subroutine save_field_hdf5(time,filename,field_out)
  use vars
  use hdf5_wrapper
  use helpers
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename
  character(len=strlen) :: fname,fname2
  integer, dimension(1:3) :: tmp
  real(kind=pr) :: t1, mbyte

  t1 = MPI_wtime()
  ! check if the file name contains the suffix *.h5
  ! if not, add it
  if (index(filename,'.h5')==0 ) then
    ! file does not contain *.h5 ending -> add suffix
    fname = trim(adjustl(filename))//'.h5'
  else
    fname = trim(adjustl(filename))
  endif

  ! header
  if (root) then
    write(*,'("Writing to ",A," dset=",A," stride=",i1," ...")',advance='no') &
    trim(adjustl(fname)), trim(adjustl(get_dsetname(fname))), striding
  endif

  if (striding<1) call abort("Striding value is bad, exit!")

  if (striding==1) then
    ! save the entire field to disk (no striding)
    ! write actual field to the file
    call write_field_hdf5( fname, get_dsetname(fname), ra, rb, field_out)
    ! append some useful attributes to the field in the file
    call write_attribute( fname, get_dsetname(fname), "time",(/time/))
    call write_attribute( fname, get_dsetname(fname), "viscosity",(/nu/))
    call write_attribute( fname, get_dsetname(fname), "epsi",(/eps/))
    call write_attribute( fname, get_dsetname(fname), "domain_size",(/xl,yl,zl/))
    call write_attribute( fname, get_dsetname(fname), "nxyz",(/nx,ny,nz/))
  else
    ! save strided field to disk
    call save_field_hdf5_strided(time,fname,field_out)
  endif

  ! footer
  t1 = MPI_wtime() - t1
  mbyte = dble(nx)*dble(ny)*dble(nz)*4.d0/1.0d+6
  if (root) write(*,'(".. wrote ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') &
  mbyte, t1, mbyte/t1
end subroutine save_field_hdf5





!-------------------------------------------------------------------------------
! save a strided field to HDF5. this is certainly not the most elegant way to do
! it, nor the most general, but it works.
!-------------------------------------------------------------------------------
! we figure out what the array bounds of the downsampled array on the local CPU
! are, then copy the data to the smaller field, and then pass both to the HDF5
! wrapper and write it to disk
!-------------------------------------------------------------------------------
subroutine save_field_hdf5_strided(time,fname,field_out)
  use helpers
  use vars
  use hdf5_wrapper
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: fname
  real(kind=pr),dimension(:,:,:), allocatable :: field_red

  integer :: ix,iy,iz,ixred,iyred,izred
  integer,dimension(1:3) :: rared, rbred
  integer :: ixmin,ixmax,ixstride,iystride,izstride,iymin,izmin,iymax,izmax

  ! do not touch lower/upper bounds, code does NOT work on arbitrary subsets.
  ! ONLY striding of the FULL field is possible.
  ixmin=0; ixmax=nx-1; ixstride=striding;
  iymin=0; iymax=ny-1; iystride=striding;
  izmin=0; izmax=nz-1; izstride=striding;

  ixred=0; rared(1) = nx*99; rbred(1) = -99
  iyred=0; rared(2) = ny*99; rbred(2) = -99
  izred=0; rared(3) = nz*99; rbred(3) = -99

  do ix = ixmin,ixmax,ixstride
    ! check if this x-coordinate is on my local memory
    if ( on_proc((/ix,ra(2),ra(3)/)) ) then
      rared(1) = min(rared(1),ixred)
      rbred(1) = max(rbred(1),ixred)
    endif
    ixred=ixred+1
  enddo

  do iy = iymin,iymax,iystride
    ! check if this y-coordinate is on my local memory
    if ( on_proc((/ra(1),iy,ra(3)/)) ) then
      rared(2) = min(rared(2),iyred)
      rbred(2) = max(rbred(2),iyred)
    endif
    iyred=iyred+1
  enddo

  do iz = izmin,izmax,izstride
    ! check if this z-coordinate is on my local memory
    if ( on_proc((/ra(1),ra(2),iz/)) ) then
      rared(3) = min(rared(3),izred)
      rbred(3) = max(rbred(3),izred)
    endif
    izred=izred+1
  enddo

  allocate( field_red(rared(1):rbred(1),rared(2):rbred(2),rared(3):rbred(3)) )

  ! copy
  do ixred = rared(1),rbred(1)
    do iyred = rared(2),rbred(2)
      do izred = rared(3),rbred(3)
        field_red(ixred,iyred,izred) = field_out(ixstride*ixred,iystride*iyred,izstride*izred)
      enddo
    enddo
  enddo

  call write_field_hdf5( fname, get_dsetname(fname), rared, rbred, field_red)
  ! append some useful attributes to the field in the file
  call write_attribute( fname, get_dsetname(fname), "time",(/time/))
  call write_attribute( fname, get_dsetname(fname), "viscosity",(/nu/))
  call write_attribute( fname, get_dsetname(fname), "epsi",(/eps/))
  call write_attribute( fname, get_dsetname(fname), "domain_size",(/xl,yl,zl/))
  call write_attribute( fname, get_dsetname(fname), "nxyz",(/nx/2,ny/2,nz/2/))

  deallocate (field_red)
end subroutine save_field_hdf5_strided
