!-------------------------------------------------------------------------------
! Warpper calling different individual time steppers
!-------------------------------------------------------------------------------
subroutine FluidTimestep(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect

  real(kind=pr)::t1
  t1=MPI_wtime()

  !-----------------------------------------------------------------------------
  ! Call fluid advancement subroutines.
  !-----------------------------------------------------------------------------
  select case(iTimeMethodFluid)
  case("RK2")
      call RK2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("RK4")
      call RK4(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("semiimplicit")
      call semiimplicit_time_stepping(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("AB2")
      if(time%it == 0) then
        call EE1(time,u,nlk,work,mask,mask_color,us,Insect,beams)
      else
        call AB2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
      endif
  case("FSI_RK2_semiimplicit")
      call FSI_RK2_semiimplicit(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("FSI_RK4_semiimplicit")
      call FSI_RK4_semiimplicit(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case default
      if (root) write(*,*) "Error! iTimeMethodFluid unknown. Abort."
      call abort()
  end select

  ! Force zero mode for mean flow
  call set_mean_flow(u,time%time)

  !-----------------------------------------------------------------------------
  ! compute unsteady corrections in every time step.
  !-----------------------------------------------------------------------------
  if (unst_corrections==1) then
    call cal_unst_corrections ( time, mask, mask_color, us, Insect )  
  endif

  time_fluid=time_fluid + MPI_wtime() - t1
end subroutine FluidTimestep




!-------------------------------------------------------------------------------
! semi implicit explicit staggered scheme for FSI simulations
! We first interpolate the pressure at the old timelevel t^n
! then advance the fluid to the next time level t^n+1 by using a suitable
! integrator. The interpolated pressure at the new time level is then the input
! to the solid solver which is advanced then
! during the fluid time step, the mask function is held constant
!-------------------------------------------------------------------------------
subroutine FSI_RK2_semiimplicit(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  ! useful error messages
  if (use_solid_model/="yes") then 
   write(*,*) "using FSI_RK2_semiimplicit without solid model?"
   call abort()
  endif
  
  !---------------------------------------------------------------------------
  ! get forces at old time level
  !---------------------------------------------------------------------------  
  ! since the BDF scheme evaluates the beam RHS only at the new time level (and
  ! not at the old one at all) the pressure at t^n is not required when using BDF2
  if ( TimeMethodSolid .ne. "BDF2") then
    call get_surface_pressure_jump (time%time, beams(1), u(:,:,:,4), timelevel="old")
  endif
  
  !---------------------------------------------------------------------------
  ! advance fluid to from (n) to (n+1)
  !---------------------------------------------------------------------------
  call RK2 (time,u,nlk,work,mask,mask_color,us,Insect,beams)
  
  !---------------------------------------------------------------------------
  ! get forces at new time level
  !---------------------------------------------------------------------------  
  call get_surface_pressure_jump (time%time, beams(1), u(:,:,:,4), timelevel="new")
  
  !---------------------------------------------------------------------------  
  ! advance solid model from (n) to (n+1)
  !---------------------------------------------------------------------------
  call SolidSolverWrapper ( time%time, time%dt_new, beams )
    
end subroutine FSI_RK2_semiimplicit


!-------------------------------------------------------------------------------
! semi implicit explicit staggered scheme for FSI simulations
! We first interpolate the pressure at the old timelevel t^n
! then advance the fluid to the next time level t^n+1 by using a suitable
! integrator. The interpolated pressure at the new time level is then the input
! to the solid solver which is advanced then
! during the fluid time step, the mask function is held constant
!-------------------------------------------------------------------------------
subroutine FSI_RK4_semiimplicit(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  ! useful error messages
  if (use_solid_model/="yes") then 
   write(*,*) "using FSI_RK2_semiimplicit without solid model?"
   call abort()
  endif
  
  !---------------------------------------------------------------------------
  ! get forces at old time level
  !---------------------------------------------------------------------------  
  ! since the BDF scheme evaluates the beam RHS only at the new time level (and
  ! not at the old one at all) the pressure at t^n is not required when using BDF2
  if ( TimeMethodSolid .ne. "BDF2") then
    call get_surface_pressure_jump (time%time, beams(1), u(:,:,:,4), timelevel="old")
  endif
  
  !---------------------------------------------------------------------------
  ! advance fluid to from (n) to (n+1)
  !---------------------------------------------------------------------------
  call RK4 (time,u,nlk,work,mask,mask_color,us,Insect,beams)
  
  !---------------------------------------------------------------------------
  ! get forces at new time level
  !---------------------------------------------------------------------------  
  call get_surface_pressure_jump (time%time, beams(1), u(:,:,:,4), timelevel="new")
  
  !---------------------------------------------------------------------------  
  ! advance solid model from (n) to (n+1)
  !---------------------------------------------------------------------------
  call SolidSolverWrapper ( time%time, time%dt_new, beams )
    
end subroutine FSI_RK4_semiimplicit

!-------------------------------------------------------------------------------
! Standard RK2 scheme for integrating the RHS defined in cal_nlk
! The scheme is also referred to as "midpoint rule"
! INPUT:
!       time: a derived datatype that contains current time, iteration number
!             time step an so on (definition see vars.f90)
!       u: The solution vector with "neq" components at time t
!       work: work array, currently unused
!       mask: the mask function (geometry of the problem)
!-------------------------------------------------------------------------------
subroutine RK2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  type(timetype) :: t
  t=time
    
  ! define the time step we want to use
  call adjust_dt(time%time,u,time%dt_new)
  
  ! now NLK2 is old velocity
  nlk(:,:,:,:,2) = u
  
  ! compute rhs at old time and make half a time step into the future
  call cal_nlk(time, u, nlk(:,:,:,:,1), work, mask, mask_color, us, Insect, beams)
  u = u + 0.5d0*time%dt_new * nlk(:,:,:,:,1)
  t%time = time%time + 0.5d0*time%dt_new
  
  
  ! now NLK1 is old velocity
  nlk(:,:,:,:,1) = nlk(:,:,:,:,2)
  
  ! using the solution at t^(n+1/2), we compute the RHs again
  call cal_nlk(t, u, nlk(:,:,:,:,2), work, mask, mask_color, us, Insect, beams)

  ! final step, advancing the midpoint derivative by a whole time step
  u = nlk(:,:,:,:,1) + time%dt_new * nlk(:,:,:,:,2)
end subroutine RK2

!-------------------------------------------------------------------------------

subroutine RK4(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  use basic_operators
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  type(timetype) :: t
  real(kind=pr)::tmp1,tmp2,tmp3,tmp4
  t = time
  
  call adjust_dt(time%time,u,time%dt_new)
  
  ! NLK 5th register holds old velocity
  nlk(:,:,:,:,5) = u
  
  !-- Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  call cal_nlk(time,u,nlk(:,:,:,:,1),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + 0.5d0*time%dt_new*nlk(:,:,:,:,1)
  t%time = time%time + 0.5d0*time%dt_new
  call cal_nlk(t,u,nlk(:,:,:,:,2),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + 0.5d0*time%dt_new*nlk(:,:,:,:,2)
  t%time = time%time + 0.5d0*time%dt_new
  call cal_nlk(t,u,nlk(:,:,:,:,3),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + time%dt_new * nlk(:,:,:,:,3)
  t%time = time%time + time%dt_new
  call cal_nlk(t,u,nlk(:,:,:,:,4),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + time%dt_new/6.d0*(nlk(:,:,:,:,1)+2.d0*nlk(:,:,:,:,2)&
      +2.d0*nlk(:,:,:,:,3)+nlk(:,:,:,:,4))
  
end subroutine RK4

!-------------------------------------------------------------------------------

subroutine EE1(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 

  call adjust_dt(time%time,u,time%dt_new)
  
  !-- Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  call cal_nlk(time,u,nlk(:,:,:,:,1),work,mask,mask_color,us,Insect,beams)
  
  !-- advance in time
  u = u + time%dt_new * nlk(:,:,:,:,1)
  
end subroutine EE1

!-------------------------------------------------------------------------------

subroutine AB2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  real(kind=pr)::b10,b11,dt1,dt0
  integer::i
  
  if (time%n0==time%n1 .or. time%n0<1 .or. time%n0>2 .or. time%n1<1 .or. time%n1>2) then
    if(mpirank==0) write(*,*) time%n0, time%n1
    if(mpirank==0) write(*,*) "In AB2, the registers n0,n1 are not well defned"
    call abort()
  endif
  
  call adjust_dt(time%time,u,time%dt_new)

  !-- Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,u,nlk(:,:,:,:,time%n0),work,mask,mask_color,us,Insect,beams)

  dt1 = time%dt_new
  dt0 = time%dt_old
  
  b10 = dt1/dt0*(0.5d0*dt1 + dt0)
  b11 = -0.5d0*dt1*dt1/dt0

  !-- Advance in time
  u = u + b10*nlk(:,:,:,:,time%n0) + b11*nlk(:,:,:,:,time%n1)
  
end subroutine AB2


!-------------------------------------------------------------------------------
! Set the time step based on the CFL condition and penalization
! stability contidion. The following limitations exist:
! 1 - CFL condition (for umax and the speed of sound, c_0)
! 2 - fixed time step dt_fixed, ignoring all other constraints, if set in params
! 3 - penalization restriction dt<eps
! 4 - maximum time step dt_max, if set in params
! 5 - dt is smaller than tsave and tintegral
!-------------------------------------------------------------------------------
subroutine adjust_dt(time,u,dt1)
  use vars
  use mpi
  use basic_operators
  implicit none

  real(kind=pr),intent(in)::time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr), intent(out)::dt1
  integer::mpicode
  real(kind=pr)::umax
  !real(kind=pr)::tmp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)

  if (dt_fixed>0.0) then
     !-- fix the time step no matter what. the result may be unstable.
     dt1=dt_fixed
     !-- stop exactly at tmax
     if (dt1 > tmax-time .and. tmax-time>0.d0) dt1=tmax-time
  else
     !-- FSI runs just need to respect CFL for velocity
     !tmp(:,:,:,:) = u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
     !umax = fieldmaxabs3(tmp)
     umax = fieldmaxabs3(u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3))
     
     !-- Adjust time step at 0th process
     if(mpirank == 0) then
        if(is_nan(umax)) then
           write(*,*) "Evolved field contains a NAN: aborting run."
           call abort()
        endif
     
        !-- Impose the CFL condition.
        if (umax >= 1.0d-8) then
           if (nx==1) then
             ! 2D case
             dt1=min(dy,dz)*cfl/umax
           else
             ! 3D case
             dt1=min(dx,dy,dz)*cfl/umax
           endif
        else
           !-- umax is very very small
           dt1=1.0d-3
        endif

        !-- Impose penalty stability condition: dt cannot be larger than eps
        if (iPenalization > 0) dt1=min(0.99*eps,dt1) 
        
        ! Don't jump past save-points: if the time-step is larger than
        ! the time interval between outputs, decrease the time-step.
        if(tsave > 0.d0 .and. dt1 > tsave) then
           dt1=tsave
        endif
        if(tintegral > 0.d0 .and. dt1 > tintegral) then
           dt1=tintegral
        endif

        ! CFL condition for speed of sound
        if (nx==1) then
          ! 2D case
          dt1 = min( dt1, min(dy,dz)*cfl/c_0 )
        else
          ! 3D case
          dt1 = min( dt1, min(dx,dy,dz)*cfl/c_0 )
        endif
        
        !-- impose max dt, if specified
        if (dt_max>0.d0) dt1=min(dt1,dt_max)
        if (dt1 > tmax-time .and. tmax-time>0.d0) dt1=tmax-time
     endif

     ! Broadcast time step to all processes
     call MPI_BCAST(dt1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
  endif
  
end subroutine adjust_dt


!-------------------------------------------------------------------------------
! Force the mean flow
!-------------------------------------------------------------------------------
subroutine set_mean_flow(u,time)
  use vars
  use ghosts
  use basic_operators
  implicit none
  
  real(kind=pr),intent(inout) :: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout) :: time
  real(kind=pr) :: uxmean_tmp,uymean_tmp,uzmean_tmp

  ! TODO: volume_integral should be modified for stretched grid
  if (iMeanFlow_x=="fixed") then 
    uxmean_tmp = volume_integral(u(:,:,:,1)) / (xl*yl*zl)
    u(:,:,:,1) = u(:,:,:,1) - uxmean_tmp + Uxmean
  endif
  if (iMeanFlow_y=="fixed") then
    uymean_tmp = volume_integral(u(:,:,:,2)) / (xl*yl*zl)
    u(:,:,:,2) = u(:,:,:,2) - uymean_tmp + Uymean
  endif
  if (iMeanFlow_z=="fixed") then
    uzmean_tmp = volume_integral(u(:,:,:,3)) / (xl*yl*zl)
    u(:,:,:,3) = u(:,:,:,3) - uzmean_tmp + Uzmean
  endif
end subroutine set_mean_flow


