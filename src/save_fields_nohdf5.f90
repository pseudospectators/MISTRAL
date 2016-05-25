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
  integer :: sz_out(1:3) ! local array size
  !integer :: error  ! error flags
  !integer :: mpierror, i

  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  t1 = MPI_wtime()
  
  ! Array bounds and sizes
  sz_out(1) = rb(1)-ra(1) +1
  sz_out(2) = rb(2)-ra(2) +1
  sz_out(3) = rb(3)-ra(3) +1

  ! Write Fortran binary file
  open(10, file = trim(adjustl(filename)), form='unformatted', access='direct', recl=sz_out(1)*sz_out(2)*sz_out(3))
  write (10,rec=1) field_out
  close (10)

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

  character(len=19) :: filename
  real(kind=pr) :: t1
  integer :: error  ! error flags
  character(len=6) :: suffix

  t1=MPI_wtime() ! performance diagnostic

  ! Create current filename:
  write(filename,'("runtime_backup",i1,".bin")') nbackup

  ! Create name for the backup file. We keep at any time at most 2 sets of 
  ! backups, runtime_backupX.h5 with X=0,1  
  if(mpirank == 0) then
     write(*,'("Dumping runtime_backup",i1,".bin (time=",es12.4,") to disk....")',&
     advance='no') nbackup, time%time
     ! Write meta data
     open (15, file = trim(adjustl(filename)), form='formatted', status='replace')     
     write (15,*) "time  dt_old  dt_new  n1  it"
     write (15,*) time%time,time%dt_old,time%dt_new,time%n1,time%it
  endif

  ! Create current filename on this mpi rank
  write(suffix,'(i0.6)') mpirank

  !-------------------------------------------------------------------------
  ! Write the fluid backup field
  !-------------------------------------------------------------------------
  open (15, file = trim(adjustl(filename))//"."//trim(adjustl(suffix)), form='unformatted', status='replace')     
  write (15) time%time,time%dt_old,time%dt_new,time%n1,time%it,u
  close (15)  

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

  if(mpirank == 0) write(*,'(A)') "...DONE!"
end subroutine dump_runtime_backup



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
  integer :: sz_out(1:3) ! local array size

  ! Array bounds and sizes
  sz_out(1) = rb(1)-ra(1) +1
  sz_out(2) = rb(2)-ra(2) +1
  sz_out(3) = rb(3)-ra(3) +1

  ! Verbose
  if (mpirank==0) then
    write (*,'("Reading file ",A,"  .....")',advance='no') trim(adjustl(filename))
  endif
  
  !-----------------------------------------------------------------------------
  ! perform tests
  !-----------------------------------------------------------------------------
  call check_file_exists ( filename )
  
  !-----------------------------------------------------------------------------
  ! load the file
  !-----------------------------------------------------------------------------  
  open(10, file = trim(adjustl(filename)), form='unformatted', access='direct', recl=sz_out(1)*sz_out(2)*sz_out(3))
  read (10,rec=1) field
  close (10)
  
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
  character(len=7) :: suffix

  if(mpirank == 0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! I'm trying to resume a backup file: "//filename
  endif
  
  ! Create current filename for this mpi rank:
  write(suffix,'(i0.6)') mpirank

  call check_file_exists ( trim(adjustl(filename))//"."//trim(adjustl(suffix)) )

  !-------------------------------------------------------------------------
  ! Read the fluid backup field
  !-------------------------------------------------------------------------
  open (15, file = trim(adjustl(filename))//"."//trim(adjustl(suffix)), form='unformatted', status='old')     
  read (15) time%time,time%dt_old,time%dt_new,time%n1,time%it,u
  close (15)  

  if(mpirank == 0) then
     write(*,'("time=",es15.8," dt0=",es15.8)') time%time, time%dt_old
     write(*,'("!!! DONE READING BACKUP (thats good news!)")')
     write(*,'("---------")')
  endif

end subroutine read_runtime_backup



!-------------------------------------------------------------------------------
! This is a stub. Postprocessing tools are not available without HDF5 support.
!-------------------------------------------------------------------------------
subroutine postprocessing()
end subroutine
