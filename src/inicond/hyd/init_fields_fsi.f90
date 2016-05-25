! Set initial conditions for fsi code.
subroutine init_fields_fsi(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use p3dfft_wrapper
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
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect

  integer :: ix,iy,iz
  real (kind=pr) :: x,y,z,r,gamma0,x00,r00,omega,a,b
  real (kind=pr) :: R1,R2,uu,r12,r22,Om,y1,y2,z1,z2,r02

  ! Assign zero values
  time%time = 0.0d0
  time%dt_new  = tsave
  time%it = 0

  u = 0.d0
  us = 0.d0
  nlk = 0.d0
  work = 0.d0
  mask = 0.d0
  mask_color = 0

  select case(inicond)
  case ("couette")
    !--------------------------------------------------
    ! couette flow
    !--------------------------------------------------
    R1=0.4d0
    R2=1.0d0
    omega=1.25d0

    a = omega*(-R1**2 / (R2**2 - R1**2))
    b = omega*(R1**2 * R2**2) / (R2**2 - R1**2)

    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy - 0.5d0*yl
        z = dble(iz)*dz - 0.5d0*zl
        R = dsqrt(y**2 + z**2)

        if ((R>R1).and.(R<R2)) then
          ! fluid domain
          uu = a*R + b/R
          u(:,iy,iz,1) = 0.d0
          u(:,iy,iz,2) =+uu*z/R
          u(:,iy,iz,3) =-uu*y/R
        elseif (R>=R2) then
          ! outer cylinder
          u(:,iy,iz,1) = 0.d0
          u(:,iy,iz,2) = 0.d0
          u(:,iy,iz,3) = 0.d0
        elseif (R<=R1) then
          ! inner cylinder
          u(:,iy,iz,1) = 0.d0
          u(:,iy,iz,2) = +omega*z
          u(:,iy,iz,3) = -omega*y
        endif
      enddo
    enddo

  case ("taylor_green_2d")
    !--------------------------------------------------
    ! taylor green vortices
    !--------------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy
        z = dble(iz)*dz
        u(:,iy,iz,1) = 0.0d0
        u(:,iy,iz,2) = dsin( y ) * dcos( z )
        u(:,iy,iz,3) =-dcos( y ) * dsin( z )
        ! remember that p is the total pressure, i.e. p+0.5u*u
        ! note: we changed the formulation, so now we have p static pressure
        u(:,iy,iz,4) =0.25d0*(dcos(2.d0*y)+dcos(2.d0*z)) !+ &
                      !0.5d0*(u(:,iy,iz,2)*u(:,iy,iz,2)+u(:,iy,iz,3)*u(:,iy,iz,3))
      enddo
    enddo

  case ("guermond_minev")
    !--------------------------------------------------
    ! Guermond-Minev Stokes flow
    !--------------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy
        z = dble(iz)*dz
        u(:,iy,iz,1) = 0.0d0
        u(:,iy,iz,2) = dsin( y ) * dsin( z )
        u(:,iy,iz,3) = dcos( y ) * dcos( z )
        ! remember that p is the total pressure
        u(:,iy,iz,4) = dcos( y ) * dsin( z )
      enddo
    enddo

  case ("dipole_keetels")
    !--------------------------------------------------
    ! Keetels et al JCP dipole
    !--------------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        ! Parameters
        Om = 20d0
        r02 = 0.2d0**2
        y1 = 1.0d0
        z1 = 1.0d0-0.3d0
        y2 = 1.0d0
        z2 = 1.0d0+0.3d0
        y = dble(iy)*dy
        z = dble(iz)*dz
        ! Distances
        r12 = (y-y1)*(y-y1)+(z-z1)*(z-z1)
        r22 = (y-y2)*(y-y2)+(z-z2)*(z-z2)
        ! Velocity components
        u(:,iy,iz,1) = 0.0d0
        u(:,iy,iz,2) = -0.5d0*Om*(z-z1)*dexp(-r12/r02)+0.5d0*Om*(z-z2)*dexp(-r22/r02)
        u(:,iy,iz,3) = 0.5d0*Om*(y-y1)*dexp(-r12/r02)-0.5d0*Om*(y-y2)*dexp(-r22/r02)
        ! remember that p is the total pressure, i.e. p+0.5u*u
        u(:,iy,iz,4) =0.0d0 + 0.5d0*(u(:,iy,iz,2)*u(:,iy,iz,2)+u(:,iy,iz,3)*u(:,iy,iz,3))
      enddo
    enddo

  case("infile")
     !--------------------------------------------------
     ! read HDF5 files
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: reading infiles"
     call Read_Single_File ( file_ux, u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )
     call Read_Single_File ( file_uy, u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2) )
     call Read_Single_File ( file_uz, u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3) )
     call Read_Single_File ( file_p , u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4) )
     if (mpirank==0) write (*,*) "*** done reading infiles"

  case("MeanFlow")
     !--------------------------------------------------
     ! mean flow only
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: mean flow"

     u(:,:,:,1) = Uxmean
     u(:,:,:,2) = Uymean
     u(:,:,:,3) = Uzmean

  case("quiescent")
     !--------------------------------------------------
     ! fluid at rest
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: fluid at rest"
     u = 0.d0

  case default
     if(inicond(1:8) == "backup::") then
        !--------------------------------------------------
        ! read from backup
        !--------------------------------------------------
        if (mpirank==0) write (*,*) "*** inicond: retaking backup " // &
             inicond(9:len(inicond))
        call read_runtime_backup(inicond(9:len(inicond)),time,u,Insect,beams)

     else
        !--------------------------------------------------
        ! unknown inicond : error
        !--------------------------------------------------
        if (mpirank==0) write (*,*) inicond
        if (mpirank==0) write (*,*) '??? ERROR: Invalid initial condition'
        call abort()
     endif
  end select

  !-----------------------------------------------------------------------------
  ! Synchronize ghost points
  !-----------------------------------------------------------------------------
  call synchronize_ghosts(u,neq)

  !-----------------------------------------------------------------------------
  ! If module is in use, initialize also the solid solver
  !-----------------------------------------------------------------------------
  if (use_solid_model=="yes") then
    if(mpirank==0) write(*,*) "Initializing solid solver and testing..."
    call init_beams( beams )
    call show_beam(beams(1))
    call surface_interpolation_testing( time%time, beams(1),work(:,:,:,1),mask,mask_color,us )
    call init_beams( beams )
  endif

tstart=time%time

end subroutine init_fields_fsi
