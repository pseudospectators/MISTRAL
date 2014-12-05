!-------------------------------------------------------------------------------
! Semi-implicit time stepping routine
! INPUT:
!       time: struct containing the current time, time step and so on
!       mask: mask function containg geometry
!       mask_color: solid parts 'color' markers
!       us: velocity inside solid
!       work: work array
!       nlk: work array
!       Insect: insect parameters
!       beam: beam parameters
! INPUT/OUTPUT:
!       u:  solution vector, 1:3 - velocity, 4 - pressure
!           input at time n, output at time n+1
!-------------------------------------------------------------------------------
subroutine semiimplicit_time_stepping(time,u,nlk,work,mask,mask_color,us,Insect,beams)
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
  integer :: chflag

  ! Set new time step size
  time%dt_old = time%dt_new
  call adjust_dt_impl(time%time,u(:,:,:,1:3),time%dt_new,chflag)

  ! Zero time step init
  if (time%it==0) then
    ! Synchronize ghost points
    call synchronize_ghosts(u,neq)
    ! p^{n-1/2} = p^0 : static pressure
    ! p is stored in u(:,:,:,4)
    ! Solve A phi^{-1/2} = -1/dt * div u^0
    ! phi is stored in nlk(:,:,:,4,2)
    call GM_pressure(time%dt_new,nlk(:,:,:,4,2),u(:,:,:,1:3),work)
  endif

  ! Pressure predictor
  ! p^{*,n+1/2} = p^{n-1/2} + phi^{n-1/2}
  nlk(:,:,:,4,3) = u(:,:,:,4) + nlk(:,:,:,4,2)

  ! Velocity update
  ! Predictor step
  ! Viscous and pressure gradient terms
  ! Stored in nlk(:,:,:,1:3,3)
  call cal_lin_expl(u(:,:,:,1:neq),nlk(:,:,:,1:neq,3))

  ! Nonlinear, penalty and forcing terms at time t (predictor)
  ! Stored in nlk(:,:,:,1:3,1)
  call cal_nl_expl(time%time,u(:,:,:,1:neq),nlk(:,:,:,1:neq,1),&
                    work,mask,mask_color,us(:,:,:,1:neq),Insect,beams)

  ! Update nlk(:,:,:,1:3,2) to store u*^{n+1/2}
  nlk(:,:,:,1:3,2) = u(:,:,:,1:3) + 0.5d0*time%dt_new*(nlk(:,:,:,1:3,1)+nlk(:,:,:,1:3,3))

  ! First stage of Douglas splitting: compute xi^{n+1}
  ! Stored in nlk(:,:,:,1:3,1)
  call cal_nl_expl(time%time+0.5d0*time%dt_new,nlk(:,:,:,1:neq,2),nlk(:,:,:,1:neq,1),&
                    work,mask,mask_color,us(:,:,:,1:neq),Insect,beams)

  ! Store u^n in nlk(:,:,:,1:3,2)
  nlk(:,:,:,1:3,2) = u(:,:,:,1:3)
  ! Update u(:,:,:,1:3) to store xi^{n+1}
  u(:,:,:,1:3) = u(:,:,:,1:3) + time%dt_new*(nlk(:,:,:,1:3,1)+nlk(:,:,:,1:3,3))

  ! Remaining 3 stages of Douglas splitting
  ! u^{n+1} is stored in u(:,:,:,1:3)
  ! u^n remains in nlk(:,:,:,1:3,2)
  call impl_lin(time%time,time%dt_new,u(:,:,:,1:neq),nlk(:,:,:,1:neq,2))

  ! Pressure corrector
  ! Solve A phi^{n+1/2} = -1/dt * div u^{n+1}
  ! phi is stored in nlk(:,:,:,4,2) at output
  ! u^{n+1} is remains in u(:,:,:,1:3)
  call GM_pressure(time%dt_new,nlk(:,:,:,4,2),u(:,:,:,1:3),work)

  ! Pressure update
  ! p^{n+1/2} = p^{n-1/2} + phi^{n+1/2} - 1/2 chi nu div(u^{n+1}+u^n)
  ! u^n is taken from nlk(:,:,:,1:3,2)
  ! u^{n+1} is taken from u(:,:,:,1:3)
  ! phi^{n+1/2} is taken from in nlk(:,:,:,4,2)
  ! p^{n-1/2} is taken from u(:,:,:,4)
  ! p^{n+1/2} is stored in u(:,:,:,4) at output
  call GM_pressure_update(nlk(:,:,:,1:neq,2),u(:,:,:,1:neq))

end subroutine semiimplicit_time_stepping


!-------------------------------------------------------------------------------
! Set the time step based on the CFL condition and penalization
! stability contidion. The following limitations exist:
! 1 - CFL condition (for umax and the speed of sound, c_0)
! 2 - fixed time step dt_fixed, ignoring all other constraints, if set in params
! 3 - penalization restriction dt<eps
! 4 - maximum time step dt_max, if set in params
! 5 - dt is smaller than tsave and tintegral
! 7 - time+dt<=tmax
! INPUT:
!      time: time variable
!      u: velocity & pressure field
! OUTPUT:
!      dt1: new time step size
!      chflag: flag, 0: dt unchanged; 1: dt changed
!-------------------------------------------------------------------------------
subroutine adjust_dt_impl(time,u,dt1,chflag)
  use vars
  use mpi
  use basic_operators
  implicit none

  real(kind=pr),intent(in)::time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr), intent(out)::dt1
  integer,intent(inout)::chflag
  integer::mpicode
  real(kind=pr)::umax

  ! Initialize time step change flag
  chflag = 0

  if (dt_fixed>0.0d0) then
     !-- fix the time step no matter what. the result may be unstable.
     dt1=dt_fixed
     !-- stop exactly at tmax
     if (dt1 > tmax-time .and. tmax-time>0.d0) then 
       dt1=tmax-time
       chflag = 1
     endif
     ! Broadcast time step to all processes
     call MPI_BCAST(dt1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
  else
     !-- here we assume that time step changes
     chflag = 1

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
        if (iPenalization > 0) dt1=min(0.99d0*eps,dt1) 
        
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
  
end subroutine adjust_dt_impl


!-------------------------------------------------------------------------------
! Linear terms of the momentum equation required for Douglas splitting scheme
! rhs = - grad(p) + nu * nabla^2 u
! computed using SECOND order finite differences.
!
! INPUT:
!       u: solution vector, u=(/ux,uy,uz,p/)
! OUTPUT:
!       rhs: linear rhs terms
!-------------------------------------------------------------------------------
subroutine cal_lin_expl(u,rhs)
  use p3dfft_wrapper
  use basic_operators
  use insect_module
  use solid_model
  use vars
  implicit none

  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)

  integer::ix,iy,iz
  real(kind=pr)::ux,uy,uz,pdx,pdy,pdz,&
  uxdxdx,uxdydy,uxdzdz,uydxdx,uydydy,uydzdz,uzdxdx,uzdydy,uzdzdz,&
  dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv

  ! Grid step inverse
  dxinv = 0.5d0/dx
  dyinv = 0.5d0/dy
  dzinv = 0.5d0/dz
  
  dx2inv = 1.d0/(dx*dx)
  dy2inv = 1.d0/(dy*dy)
  dz2inv = 1.d0/(dz*dz)
  
  if (nx==1) then
    ! 2d case
    ix=0
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        ! Velocity components
        uy = u(ix,iy,iz,2)
        uz = u(ix,iy,iz,3)
        
        ! Pressure gradient 
        pdy  = (rhs(ix,iy+1,iz,4) - rhs(ix,iy-1,iz,4))*dyinv
        pdz  = (rhs(ix,iy,iz+1,4) - rhs(ix,iy,iz-1,4))*dzinv
      
        ! Laplacian 
        uydydy = (u(ix,iy-1,iz,2)-2.d0*u(ix,iy,iz,2)+u(ix,iy+1,iz,2))*dy2inv 
        uydzdz = (u(ix,iy,iz-1,2)-2.d0*u(ix,iy,iz,2)+u(ix,iy,iz+1,2))*dz2inv 
      
        uzdydy = (u(ix,iy-1,iz,3)-2.d0*u(ix,iy,iz,3)+u(ix,iy+1,iz,3))*dy2inv 
        uzdzdz = (u(ix,iy,iz-1,3)-2.d0*u(ix,iy,iz,3)+u(ix,iy,iz+1,3))*dz2inv 

        ! RHS
        rhs(ix,iy,iz,1) = 0.d0
        rhs(ix,iy,iz,2) = - pdy + nu*(uydydy+uydzdz) 
        rhs(ix,iy,iz,3) = - pdz + nu*(uzdydy+uzdzdz)
      enddo
    enddo     
  else
    ! 3d case 
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          ! Velocity components
          ux = u(ix,iy,iz,1)
          uy = u(ix,iy,iz,2)
          uz = u(ix,iy,iz,3)
        
          ! Pressure gradient
          pdx = (rhs(ix+1,iy,iz,4) - rhs(ix-1,iy,iz,4))*dxinv
          pdy = (rhs(ix,iy+1,iz,4) - rhs(ix,iy-1,iz,4))*dyinv
          pdz = (rhs(ix,iy,iz+1,4) - rhs(ix,iy,iz-1,4))*dzinv
        
          ! Laplacian
          uxdxdx = (u(ix-1,iy,iz,1)-2.d0*u(ix,iy,iz,1)+u(ix+1,iy,iz,1))*dx2inv 
          uxdydy = (u(ix,iy-1,iz,1)-2.d0*u(ix,iy,iz,1)+u(ix,iy+1,iz,1))*dy2inv 
          uxdzdz = (u(ix,iy,iz-1,1)-2.d0*u(ix,iy,iz,1)+u(ix,iy,iz+1,1))*dz2inv 
        
          uydxdx = (u(ix-1,iy,iz,2)-2.d0*u(ix,iy,iz,2)+u(ix+1,iy,iz,2))*dx2inv 
          uydydy = (u(ix,iy-1,iz,2)-2.d0*u(ix,iy,iz,2)+u(ix,iy+1,iz,2))*dy2inv 
          uydzdz = (u(ix,iy,iz-1,2)-2.d0*u(ix,iy,iz,2)+u(ix,iy,iz+1,2))*dz2inv 
        
          uzdxdx = (u(ix-1,iy,iz,3)-2.d0*u(ix,iy,iz,3)+u(ix+1,iy,iz,3))*dx2inv 
          uzdydy = (u(ix,iy-1,iz,3)-2.d0*u(ix,iy,iz,3)+u(ix,iy+1,iz,3))*dy2inv 
          uzdzdz = (u(ix,iy,iz-1,3)-2.d0*u(ix,iy,iz,3)+u(ix,iy,iz+1,3))*dz2inv 

          ! RHS
          rhs(ix,iy,iz,1) = - pdx + nu*(uxdxdx+uxdydy+uxdzdz) 
          rhs(ix,iy,iz,2) = - pdy + nu*(uydxdx+uydydy+uydzdz) 
          rhs(ix,iy,iz,3) = - pdz + nu*(uzdxdx+uzdydy+uzdzdz) 
        enddo
      enddo
    enddo     
  endif
  
  ! Synchronize ghost points
  call synchronize_ghosts(rhs(:,:,:,:),3)
 
end subroutine cal_lin_expl


!-------------------------------------------------------------------------------
! Nonlinear, penality and forcing terms of the momentum equation, explicit
! rhs = vor \cross u - chi/eta*(u-us) + f
! Tis is a wrapper for separate 2d and 3d subroutines
! and 2nd or 4th order schemes for the nonlinear term.
!
! INPUT:
!       time: time variable
!       u: solution vector, u=(/ux,uy,uz,p/)
!       mask: mask function containg geometry
!       mask_color: solid parts 'color' markers
!       us: velocity inside solid
!       work: work array
!       Insect: insect parameters
!       beam: beam parameters
! OUTPUT:
!       rhs: rhs terms
!-------------------------------------------------------------------------------
subroutine cal_nl_expl(time,u,rhs,work,mask,mask_color,us,Insect,beams)
  use p3dfft_wrapper
  use basic_operators
  use insect_module
  use solid_model
  use vars
  implicit none

  real(kind=pr),intent(in)::time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout)::beams
  type(diptera), intent(inout)::Insect 

  ! Update mask function to ensure it is at the right time
  if ((iMoving==1).and.(iPenalization==1)) then
    call create_mask(time,mask,mask_color,us,Insect,beams)
  endif

  ! Call specialized rhs subroutines
  if (nx==1) then
    ! 2d case, 2nd order finite differences
    call rhs_impl_2nd_2D(time,u,rhs,work,mask,mask_color,us,Insect,beams)
  else 
    ! 3d case, 2nd order finite differences
    call rhs_impl_2nd(time,u,rhs,work,mask,mask_color,us,Insect,beams)
  endif
  
  ! Synchronize ghost points
  call synchronize_ghosts(rhs(:,:,:,:),3)
 
end subroutine cal_nl_expl

!-------------------------------------------------------------------------------
! Nonlinear, penality and forcing terms of the momentum equation, explicit
! 3d case, 2nd order finite differences
! INPUT/OUTPUT as in cal_nl_expl 
!-------------------------------------------------------------------------------
subroutine rhs_impl_2nd(time,u,rhs,work,mask,mask_color,us,Insect,beams)
  use vars
  use insect_module
  use solid_model
  use ghosts
  
  implicit none
  real(kind=pr),intent(in)::time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid),dimension(1:nBeams),intent(inout)::beams
  type(diptera),intent(inout)::Insect 
  
  integer::ix,iy,iz
  real(kind=pr)::ux,uy,uz,vorx,vory,vorz,uxdy,uxdz,uydx,uydz,uzdx,uzdy,&
  dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv,penalx,penaly,penalz
  real(kind=pr)::forcing(1:3)
  type(timetype)::t
  
  ! fetch forcing term used to accelerate the mean flow
  t%time = time
  t%dt_new = 0.d0
  call forcing_term(t,u,forcing)
 
  ! Inverse grid step 
  dxinv = 1.d0/(2.d0*dx)
  dyinv = 1.d0/(2.d0*dy)
  dzinv = 1.d0/(2.d0*dz)
  
  dx2inv = 1.d0/(dx*dx)
  dy2inv = 1.d0/(dy*dy)
  dz2inv = 1.d0/(dz*dz)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ! Velocity components
        ux = u(ix,iy,iz,1)
        uy = u(ix,iy,iz,2)
        uz = u(ix,iy,iz,3)
        
        ! Velocity gradient components
        uxdy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dyinv
        uxdz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dzinv
        
        uydx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dxinv
        uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
        
        uzdx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dxinv
        uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv
        
        ! Vorticity
        vorx = uzdy - uydz
        vory = uxdz - uzdx
        vorz = uydx - uxdy
        
        ! Penalty term
        penalx = -mask(ix,iy,iz)*(ux-us(ix,iy,iz,1))
        penaly = -mask(ix,iy,iz)*(uy-us(ix,iy,iz,2))
        penalz = -mask(ix,iy,iz)*(uz-us(ix,iy,iz,3))
      
        ! Nonlinear, penalty and forcing terms
        rhs(ix,iy,iz,1) = uy*vorz -uz*vory + penalx + forcing(1)
        rhs(ix,iy,iz,2) = uz*vorx -ux*vorz + penaly + forcing(2)
        rhs(ix,iy,iz,3) = ux*vory -uy*vorx + penalz + forcing(3)
      enddo
    enddo
  enddo     
        
end subroutine rhs_impl_2nd


!-------------------------------------------------------------------------------
! Nonlinear, penality and forcing terms of the momentum equation, explicit
! 2d case, 2nd order finite differences
! INPUT/OUTPUT as in cal_nl_expl 
!-------------------------------------------------------------------------------
subroutine rhs_impl_2nd_2d(time,u,rhs,work,mask,mask_color,us,Insect,beams)
  use vars
  use insect_module
  use solid_model
  use ghosts
  
  implicit none
  real(kind=pr),intent(in)::time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid),dimension(1:nBeams),intent(inout)::beams
  type(diptera),intent(inout)::Insect 
  
  integer::ix,iy,iz
  real(kind=pr)::uy,uz,vorx,uydz,uzdy,&
  dyinv,dzinv,dy2inv,dz2inv,penaly,penalz
  real(kind=pr)::forcing(1:3)
  type(timetype)::t

  ! fetch forcing term used to accelerate the mean flow
  t%time = time
  t%dt_new = 0.0d0
  call forcing_term(t,u,forcing)

  ! Inverse grid step
  dyinv = 1.d0/(2.d0*dy)
  dzinv = 1.d0/(2.d0*dz)
  
  dy2inv = 1.d0/(dy*dy)
  dz2inv = 1.d0/(dz*dz)
  
  ix=0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      ! Velocity components
      uy = u(ix,iy,iz,2)
      uz = u(ix,iy,iz,3)
      
      ! Velocity gradient
      uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
      uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv
      
      ! Vorticity
      vorx = uzdy - uydz
      
      ! Penalty term
      penaly = -mask(ix,iy,iz)*(uy-us(ix,iy,iz,2))
      penalz = -mask(ix,iy,iz)*(uz-us(ix,iy,iz,3))
 
      ! Sum of nonlinear, penality and forcing terms
      rhs(ix,iy,iz,1) = 0.d0
      rhs(ix,iy,iz,2) = +uz*vorx + penaly + forcing(2)
      rhs(ix,iy,iz,3) = -uy*vorx + penalz + forcing(3)
    enddo
  enddo     

end subroutine rhs_impl_2nd_2d


!-------------------------------------------------------------------------------
! Implicit directional splitting momentum equation solver
! INPUT:
!      time: current time
!      dt: time step size
!      u0: velocity field at the beginning of Douglas splitting
! INPUT/OUTPUT:
!      u: velocity field, 3 components
!-------------------------------------------------------------------------------
subroutine impl_lin(time,dt,u,u0)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(in) :: time,dt
  real(kind=pr),intent(inout) :: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout) :: u0(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  ! Local variables
  integer :: ix,iy,iz,ic,mpicommdir,mpiszdir,radir,rbdir,gadir,gbdir,mpirankdir
  integer :: bcipivy(2*mpidims(2)),bcipivz(2*mpidims(1))
  real(kind=pr) :: h2inv,difcoef,det0
  real(kind=pr) :: vly(ra(2):rb(2)),vry(ra(2):rb(2)),&
                   utmpy(ga(2):gb(2)),u0tmpy(ga(2):gb(2)),&
                   bcmaty(2*mpidims(2),2*mpidims(2)),&
                   cndiagy(ra(2):rb(2),1:2)
  real(kind=pr) :: vlz(ra(3):rb(3)),vrz(ra(3):rb(3)),&
                   utmpz(ga(3):gb(3)),u0tmpz(ga(3):gb(3)),&
                   bcmatz(2*mpidims(1),2*mpidims(1)),&
                   cndiagz(ra(3):rb(3),1:2)
  real(kind=pr) :: vlx(ra(1):rb(1)),vrx(ra(1):rb(1)),&
                   utmpx(ga(1):gb(1)),u0tmpx(ga(1):gb(1)),&
                   cndiagx(ra(1):rb(1),1:2)

  ! This subroutine assumes that the domain decomposition is 2D
  ! The domain is NOT split in x direction
  ! In y direction, it is split in mpidims(2) parts
  ! In z direction, it is split in mpidims(1) parts

  ! Set viscosity
  difcoef = nu

  ! X DIRECTION
  ! Only if 3d
  if (nx.ne.1) then
    ! Set up parameters in x direction
    h2inv =  1.d0/(dx*dx)
    radir = ra(1)
    rbdir = rb(1)
    gadir = ga(1)
    gbdir = gb(1)
    ! Serial 1d solver init
    call impl_lin_1d_serial_init (h2inv,radir,rbdir,dt,&
                                  difcoef,det0,cndiagx,vlx,vrx)
    ! Loop for all fields
    do ic = 1,3
      ! Loop for all lines x=const. This is local.
      do iz=ga(3),gb(3)
        !zz = dble(iz)*dz
        do iy=ga(2),gb(2)
          !xx = dble(ix)*dx 
          ! Vector to be processed
          utmpx(:) = u(gadir:gbdir,iy,iz,ic)
          ! Constant vector
          u0tmpx(:) = u0(gadir:gbdir,iy,iz,ic)
          ! Solve linear system
          call impl_lin_1d_serial_solver (h2inv,radir,rbdir,gadir,gbdir,dt,&
                                          difcoef,det0,cndiagx,vlx,vrx,utmpx,u0tmpx)
          ! Vector returned
          u(radir:rbdir,iy,iz,ic) = utmpx(radir:rbdir)
        enddo
      enddo    
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_x_serial (u(:,:,:,:),3)
  endif

  ! Y DIRECTION
  ! Set up parameters in y direction
  mpicommdir = mpicommy
  mpiszdir = mpidims(2)
  h2inv =  1.d0/(dy*dy)
  radir = ra(2)
  rbdir = rb(2)
  gadir = ga(2)
  gbdir = gb(2)
  ! Cases if # subdomains = 1 or >=2
  if (mpiszdir>1) then 
    ! Parallel 1d solver init
    call impl_lin_1d_mpi_init (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,dt,&
                               difcoef,bcmaty,cndiagy,bcipivy,vly,vry)
    ! Loop for 3 velocity components
    do ic = 1,3
      ! Loop for all lines y=const
      do iz = ga(3),gb(3)
        !zz = dble(iz)*dz
        do ix = ga(1),gb(1)
          !xx = dble(ix)*dx 
          ! Vector to be processed
          utmpy(:) = u(ix,gadir:gbdir,iz,ic)
          ! Constant vector
          u0tmpy(:) = u0(ix,gadir:gbdir,iz,ic)
          ! Solve linear system
          call impl_lin_1d_mpi_solver (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,gadir,gbdir,dt,&
                                       difcoef,bcmaty,cndiagy,bcipivy,vly,vry,utmpy,u0tmpy)
          ! Vector returned
          u(ix,radir:rbdir,iz,ic) = utmpy(radir:rbdir)
        enddo 
      enddo
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_y_mpi (u(:,:,:,:),3)
  else
    ! Serial 1d solver init
    call impl_lin_1d_serial_init (h2inv,radir,rbdir,dt,&
                                  difcoef,det0,cndiagy,vly,vry)
    ! Loop for 3 velocity components
    do ic = 1,3
      ! Loop for all lines y=const
      do iz=ga(3),gb(3)
        !zz = dble(iz)*dz
        do ix=ga(1),gb(1)
          !xx = dble(ix)*dx 
          ! Vector to be processed
          utmpy(:) = u(ix,gadir:gbdir,iz,ic)
          ! Constant vector
          u0tmpy(:) = u0(ix,gadir:gbdir,iz,ic)
          ! Solve linear system
          call impl_lin_1d_serial_solver (h2inv,radir,rbdir,gadir,gbdir,dt,&
                                          difcoef,det0,cndiagy,vly,vry,utmpy,u0tmpy)
          ! Vector returned
          u(ix,radir:rbdir,iz,ic) = utmpy(radir:rbdir)
        enddo
      enddo    
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_y_serial (u(:,:,:,:),3)
  endif

  ! Z DIRECTION
  ! Set up parameters in z direction
  mpicommdir = mpicommz
  mpiszdir = mpidims(1)
  h2inv =  1.d0/(dz*dz)
  radir = ra(3)
  rbdir = rb(3)
  gadir = ga(3)
  gbdir = gb(3)
  ! Cases if # subdomains = 1 or >=2
  if (mpiszdir>1) then 
    ! Parallel 1d solver init
    call impl_lin_1d_mpi_init (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,dt,&
                               difcoef,bcmatz,cndiagz,bcipivz,vlz,vrz)
    ! Loop for 3 velocity components
    do ic = 1,3
      ! Loop for all lines z=const
      do iy=ga(2),gb(2)
        !yy = dble(iy)*dy
        do ix=ga(1),gb(1)
          !xx = dble(ix)*dx 
          ! Vector to be processed
          utmpz(:) = u(ix,iy,gadir:gbdir,ic)
          ! Constant vector
          u0tmpz(:) = u0(ix,iy,gadir:gbdir,ic)
          ! Solve linear system
          call impl_lin_1d_mpi_solver (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,gadir,gbdir,dt,&
                                       difcoef,bcmatz,cndiagz,bcipivz,vlz,vrz,utmpz,u0tmpz)
          ! Vector returned
          u(ix,iy,radir:rbdir,ic) = utmpz(radir:rbdir)
        enddo 
      enddo
    enddo
  else
    ! Serial 1d solver init
    call impl_lin_1d_serial_init (h2inv,radir,rbdir,dt,&
                                  difcoef,det0,cndiagz,vlz,vrz)
    ! Loop for all fields
    do ic = 1,3
      ! Loop for all lines z=const
      do iy=ga(2),gb(2)
        !yy = dble(iy)*dy
        do ix=ga(1),gb(1)
          !xx = dble(ix)*dx 
          ! Vector to be processed
          utmpz(:) = u(ix,iy,gadir:gbdir,ic)
          ! Constant vector
          u0tmpz(:) = u0(ix,iy,gadir:gbdir,ic)
          ! Solve linear system
          call impl_lin_1d_serial_solver (h2inv,radir,rbdir,gadir,gbdir,dt,&
                                          difcoef,det0,cndiagz,vlz,vrz,utmpz,u0tmpz)
          ! Vector returned
          u(ix,iy,radir:rbdir,ic) = utmpz(radir:rbdir)
        enddo
      enddo    
    enddo
  endif

  ! Synchronize ghost points
  call synchronize_ghosts (u(:,:,:,:),3)

end subroutine impl_lin


!-------------------------------------------------------------------------------
! Factorize full matrix
! INPUT:
!      nn: matrix size
! INPUT/OUTPUT:
!      mat: nn x nn, square full matrix (in) or factorization (out)
! OUTPUT:
!      ipiv: pivot element indices
!-------------------------------------------------------------------------------
subroutine factorize_loc1d(mat,ipiv,nn)
  use vars
  implicit none
  integer,intent(in) :: nn
  integer,intent(inout) :: ipiv(1:nn)
  real(kind=pr),intent(inout) :: mat(1:nn,1:nn)
  real(kind=pr) :: t0
  integer :: error

  t0 = MPI_wtime()
  call dgetrf (nn,nn,mat,nn,ipiv,error)
  if (error .ne. 0) then
    write(*,*) "!!! Fatal: dgetrf error.", error
    call abort()
  endif
  time_LAPACK = time_LAPACK + MPI_wtime() - t0
end subroutine factorize_loc1d


!-------------------------------------------------------------------------------
! Solve linear system with a full matrix, factorized
! INPUT:
!      nn: matrix size
!      mat: nn x nn, square full matrix LDU factorization
!      ipiv: pivot element indices
!      rhs: right-hand size
! OUTPUT:
!      x: solution of mat*x=rhs 
!-------------------------------------------------------------------------------
subroutine solve_loc1d (mat,ipiv,rhs,x,nn)
  use vars
  implicit none
  integer,intent(in) :: nn
  integer,intent(inout) :: ipiv(1:nn)
  real(kind=pr),intent(inout) :: mat(1:nn,1:nn)
  real(kind=pr),intent(inout) :: x(1:nn)
  real(kind=pr),intent(inout) :: rhs(1:nn)
  real(kind=pr) :: t0
  integer :: error 

  t0 = MPI_wtime()
  x = rhs
  call dgetrs ('N',nn,1,mat,nn,ipiv,x,nn,error)
  if (error .ne. 0) then 
    write(*,*) "!!! Fatal: dgetrs error.", error
    call abort()
  endif
  time_LAPACK = time_LAPACK + MPI_wtime() - t0
end subroutine solve_loc1d


!-------------------------------------------------------------------------------
! Factorize symmetric tridiagonal matrix
! INPUT:
!      nn: matrix size
! INPUT/OUTPUT:
!      diag: matrix (in) or factorization (out)
!-------------------------------------------------------------------------------
subroutine factorize_tri_loc1d(diag,nn)
  use vars
  implicit none
  integer,intent(in) :: nn
  real(kind=pr),intent(inout) :: diag(1:nn,1:2)
  real(kind=pr) :: t0
  integer :: error

  t0 = MPI_wtime()
  call dpttrf (nn,diag(1:nn,1),diag(1:nn-1,2),error)
  if (error .ne. 0) then
    write(*,*) "!!! Fatal: dpttrf error.", error
    call abort()
  endif
  time_LAPACK = time_LAPACK + MPI_wtime() - t0
end subroutine factorize_tri_loc1d


!-------------------------------------------------------------------------------
! Solve linear system with a symmetric tridiagonal matrix, factorized
! INPUT:
!      nn: matrix size
!      diag: matrix factorization
!      rhs: right-hand size
! OUTPUT:
!      x: solution of diag*x=rhs 
!-------------------------------------------------------------------------------
subroutine solve_tri_loc1d (diag,rhs,x,nn)
  use vars
  implicit none
  integer,intent(in) :: nn
  real(kind=pr),intent(inout) :: diag(1:nn,1:2)
  real(kind=pr),intent(inout) :: x(1:nn)
  real(kind=pr),intent(inout) :: rhs(1:nn)
  real(kind=pr) :: t0
  integer :: error 

  t0 = MPI_wtime()
  x = rhs
  call dpttrs (nn,1,diag(1:nn,1),diag(1:nn-1,2),x,nn,error)
  if (error .ne. 0) then 
    write(*,*) "!!! Fatal: dpttrs error.", error
    call abort()
  endif
  time_LAPACK = time_LAPACK + MPI_wtime() - t0
end subroutine solve_tri_loc1d


!-------------------------------------------------------------------------------
! Momentum equation LOD splitting. Initialization of the 1d implicit MPI solver
!-------------------------------------------------------------------------------
subroutine impl_lin_1d_mpi_init(mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,dt,&
                               difcoef,bcmat,cndiag,bcipiv,vl,vr)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: mpicommdir,mpiszdir,mpirankdir,radir,rbdir
  integer,intent(inout) :: bcipiv(2*mpiszdir)
  real(kind=pr),intent(in) :: h2inv,dt,difcoef
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),bcmat(2*mpiszdir,2*mpiszdir),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! Local variables
  integer :: nn,j,mpicode
  real(kind=pr) :: sendfoo(4),recvfoo(4*mpiszdir)
  real(kind=pr) :: vl1(mpiszdir),vlN(mpiszdir),vr1(mpiszdir),vrN(mpiszdir)
  real(kind=pr) :: rhs(radir:rbdir)
  ! Get local ranks in the line
  call MPI_COMM_RANK(mpicommdir,mpirankdir,mpicode)
  ! Crank-Nicolson matrix in x direction
  cndiag(:,:) = 0.d0
  do j = radir,rbdir
    cndiag(j,1) = 1.d0 + dt*difcoef*h2inv
  enddo
  do j = radir,rbdir-1
    cndiag(j,2) = - 0.5d0*dt*difcoef*h2inv
  enddo
  ! Factorize the CN matrix
  nn = rbdir-radir+1
  call factorize_tri_loc1d (cndiag,nn)
  ! Boundary conditions for domain decomposition
  ! BC influence basis
  rhs(:) = 0.d0
  rhs(radir) = 1.d0
  call solve_tri_loc1d (cndiag,rhs,vl,nn)
  vl(:) = (-0.5d0*dt*difcoef*h2inv)*vl(:)
  rhs(rbdir) = 1.d0
  rhs(radir) = 0.d0
  call solve_tri_loc1d (cndiag,rhs,vr,nn)
  vr(:) = (-0.5d0*dt*difcoef*h2inv)*vr(:)
  ! BC influence matrix
  ! It is only stored by one process
  ! Communicate values at the interface to rank 0
  sendfoo(1) = vl(radir)
  sendfoo(2) = vl(rbdir)
  sendfoo(3) = vr(radir)
  sendfoo(4) = vr(rbdir)
  call MPI_GATHER (sendfoo,4,MPI_DOUBLE_PRECISION,recvfoo,4,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  do j = 1,mpiszdir
    vl1(j) = recvfoo(4*j-3)
    vlN(j) = recvfoo(4*j-2)
    vr1(j) = recvfoo(4*j-1)
    vrN(j) = recvfoo(4*j)
  enddo
  ! BC influence matrix is only stored by one process
  if (mpirankdir == 0) then
    bcmat(:,:) = 0.d0
    do j = 1,mpiszdir
        bcmat(2*j-1,2*j-1) = 1.d0
        bcmat(2*j-1,modulo(2*j-3,2*mpiszdir)+1) = vl1(j)
        bcmat(2*j-1,modulo(2*j,2*mpiszdir)+1) = vr1(j)
        bcmat(2*j,2*j) = 1.d0
        bcmat(2*j,modulo(2*j-3,2*mpiszdir)+1) = vlN(j)
        bcmat(2*j,modulo(2*j,2*mpiszdir)+1) = vrN(j)
    enddo   
    ! Factorize the BC influence matrix
    call factorize_loc1d (bcmat,bcipiv,2*mpiszdir)
  endif
end subroutine impl_lin_1d_mpi_init


!-------------------------------------------------------------------------------
! Momentum equation LOD splitting. 1d MPI solver
!-------------------------------------------------------------------------------
subroutine impl_lin_1d_mpi_solver(mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,gadir,gbdir,dt,&
                                 difcoef,bcmat,cndiag,bcipiv,vl,vr,utmp,u0tmp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: mpicommdir,mpiszdir,mpirankdir,radir,rbdir,gadir,gbdir
  integer,intent(inout) :: bcipiv(2*mpiszdir)
  real(kind=pr),intent(in) :: h2inv,dt,difcoef
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),bcmat(2*mpiszdir,2*mpiszdir),&
                                 vl(radir:rbdir),vr(radir:rbdir),utmp(gadir:gbdir),u0tmp(gadir:gbdir)
  ! local variables
  integer :: j,mpicode
  real(kind=pr) :: bcxl,bcxr,shortfoo(2*mpiszdir),longfoo(4*mpiszdir)
  real(kind=pr) :: bcrhs(2*mpiszdir),bcx(2*mpiszdir),&
                   rhs(radir:rbdir),vf(radir:rbdir),bcxls(mpiszdir),&
                   bcxrs(mpiszdir),vf1(mpiszdir),vfN(mpiszdir)
  
  ! Crank-Nicolson explicit part
  rhs(:) = utmp(radir:rbdir)-0.5d0*dt*difcoef*(u0tmp((radir-1):(rbdir-1))-2.d0*u0tmp(radir:rbdir)+u0tmp((radir+1):(rbdir+1)))*h2inv
  ! Solve local system
  call solve_tri_loc1d (cndiag,rhs,vf,rbdir-radir+1)
  ! Communicate rhs to rank 0 in the line
  shortfoo(1) = vf(radir)
  shortfoo(2) = vf(rbdir)
  call MPI_GATHER (shortfoo,2,MPI_DOUBLE_PRECISION,longfoo,2,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  do j = 1,mpiszdir
    vf1(j) = longfoo(2*j-1)
    vfN(j) = longfoo(2*j)
  enddo
  ! BC influence RHS
  if (mpirankdir == 0) then
    do j = 1,mpiszdir
      bcrhs(2*j-1) = vf1(j)
      bcrhs(2*j) = vfN(j)
    enddo
    ! Solve BC influence system
    call solve_loc1d (bcmat,bcipiv,bcrhs,bcx,2*mpiszdir)
    ! Rearrange for mpi scatter
    do j = 1,mpiszdir
      bcxls(j) = bcx(modulo(2*j-3,2*mpiszdir)+1)
      bcxrs(j) = bcx(modulo(2*j,2*mpiszdir)+1)
    enddo
  endif
  ! Scatter from rank 0 in the line to all ranks
  do j = 1,mpiszdir
    longfoo(2*j-1) = bcxls(j)
    longfoo(2*j) = bcxrs(j)
  enddo
  call MPI_SCATTER (longfoo,2,MPI_DOUBLE_PRECISION,shortfoo,2,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  bcxl = shortfoo(1)
  bcxr = shortfoo(2)
  ! Superpose local solution and BC influence
  utmp(radir:rbdir) = vf(:)-bcxl*vl(:)-bcxr*vr(:)
end subroutine impl_lin_1d_mpi_solver


!-------------------------------------------------------------------------------
! Momentum equation LOD splitting. 1d serial solver initialization
!-------------------------------------------------------------------------------
subroutine impl_lin_1d_serial_init(h2inv,radir,rbdir,dt,&
                                  difcoef,det0,cndiag,vl,vr)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: radir,rbdir
  real(kind=pr),intent(in) :: h2inv,dt,difcoef
  real(kind=pr),intent(out) :: det0
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! local variables
  integer :: j,nn
  real(kind=pr) :: rhs(radir:rbdir)

  ! Crank-Nicolson matrix in x direction
  cndiag(:,:) = 0.d0
  do j = radir,rbdir
    cndiag(j,1) = 1.d0 + dt*difcoef*h2inv
  enddo
  do j = radir,rbdir-1
    cndiag(j,2) = - 0.5d0*dt*difcoef*h2inv
  enddo
  ! Factorize the CN matrix
  nn = rbdir-radir+1
  call factorize_tri_loc1d (cndiag,nn)
  ! Boundary conditions for domain decomposition
  ! BC influence basis
  rhs(:) = 0.d0
  rhs(radir) = 1.d0
  call solve_tri_loc1d (cndiag,rhs,vl,nn)
  vl(:) = (-0.5d0*dt*difcoef*h2inv)*vl(:)
  rhs(rbdir) = 1.d0
  rhs(radir) = 0.d0
  call solve_tri_loc1d (cndiag,rhs,vr,nn)
  vr(:) = (-0.5d0*dt*difcoef*h2inv)*vr(:)
  ! Compute determinant of the BC matrix
  det0 = (1.0d0+vr(radir))*(1.0d0+vl(rbdir))-vr(rbdir)*vl(radir)
end subroutine impl_lin_1d_serial_init


!-------------------------------------------------------------------------------
! Momentum equation LOD splitting. 1d serial solver
!-------------------------------------------------------------------------------
subroutine impl_lin_1d_serial_solver(h2inv,radir,rbdir,gadir,gbdir,dt,&
                                    difcoef,det0,cndiag,vl,vr,utmp,u0tmp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: radir,rbdir,gadir,gbdir
  real(kind=pr),intent(in) :: h2inv,dt,difcoef,det0
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),utmp(gadir:gbdir),u0tmp(gadir:gbdir),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! local variables
  real(kind=pr) :: detr,detl
  real(kind=pr) :: rhs(radir:rbdir),vf(radir:rbdir)

  ! Crank-Nicolson explicit part
  rhs(:) = utmp(radir:rbdir)-0.5d0*dt*difcoef*(u0tmp((radir-1):(rbdir-1))-2.d0*u0tmp(radir:rbdir)+u0tmp((radir+1):(rbdir+1)))*h2inv
  ! Solve local system
  call solve_tri_loc1d (cndiag,rhs,vf,rbdir-radir+1)
  ! Solve boundary correction
  detr = vf(radir)*(1+vl(rbdir))-vf(rbdir)*vl(radir)
  detl = vf(rbdir)*(1+vr(radir))-vf(radir)*vr(rbdir)
  ! Corrected vector
  utmp(radir:rbdir) = vf - (detl*vl+detr*vr)/det0
end subroutine impl_lin_1d_serial_solver


!-------------------------------------------------------------------------------
! Pressure update
! p^{n+1/2} = p^{n-1/2} + phi^{n+1/2} - 1/2 chi nu div(u^{n+1}+u^n)
! u^n is taken from uold uold(:,:,:,1:3)
! u^{n+1} is taken from u(:,:,:,1:3)
! phi^{n+1/2} is taken from uold(:,:,:,4)
! p^{n-1/2} is taken from u(:,:,:,4)
! p^{n+1/2} is returned in u(:,:,:,4)
!-------------------------------------------------------------------------------
subroutine GM_pressure_update(uold,u)
  use p3dfft_wrapper
  use basic_operators
  use insect_module
  use solid_model
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout) :: uold(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout) :: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  ! Local variables
  integer :: ix,iy,iz
  real(kind=pr) :: dxinv,dyinv,dzinv,coef

  ! 1/2 chi nu
  ! assuming chi=0.5
  coef = 0.25d0*nu

  ! u^{n+1}+u^n
  uold(:,:,:,1:3) = uold(:,:,:,1:3) + u(:,:,:,1:3)

  ! Inverse of grid step
  dxinv = 0.5d0/dx
  dyinv = 0.5d0/dy
  dzinv = 0.5d0/dz

  ! P plus phi plus divergence of the velocity terms
  if (nx==1) then
    ! 2D case
    ix=ra(1)
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        u(ix,iy,iz,4) = u(ix,iy,iz,4) + uold(ix,iy,iz,4) - coef *&
                    ((uold(ix,iy+1,iz,2) - uold(ix,iy-1,iz,2))*dyinv +&
                     (uold(ix,iy,iz+1,3) - uold(ix,iy,iz-1,3))*dzinv)
      enddo
    enddo
  else
    ! 3D case
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          u(ix,iy,iz,4) = u(ix,iy,iz,4) + uold(ix,iy,iz,4) - coef *&
                      ((uold(ix+1,iy,iz,1) - uold(ix-1,iy,iz,1))*dxinv +&
                       (uold(ix,iy+1,iz,2) - uold(ix,iy-1,iz,2))*dyinv +&
                       (uold(ix,iy,iz+1,3) - uold(ix,iy,iz-1,3))*dzinv)
        enddo
      enddo
    enddo
  endif

  ! Synchronize ghost points
  call synchronize_ghosts(u(:,:,:,4))
 
end subroutine GM_pressure_update


!-------------------------------------------------------------------------------
! Guermond-Minev equation for pressure correction phi
! (I-Dxx)(I-Dyy)(I-Dzz)phi = -1/dt* div u
! time : current time
! fld : unknown pressure correction field phi
! u : 3 velocity components
! work : work array
!-------------------------------------------------------------------------------
subroutine GM_pressure(dt,fld,u,work)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(in) :: dt
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout) :: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout) :: work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  ! Local variables
  integer :: ix,iy,iz,mpicommdir,mpiszdir,radir,rbdir,gadir,gbdir,mpirankdir
  integer :: bcipivy(2*mpidims(2)),bcipivz(2*mpidims(1))
  real(kind=pr) :: h2inv,dtinv,det0,dxinv,dyinv,dzinv
  real(kind=pr) :: vly(ra(2):rb(2)),vry(ra(2):rb(2)),&
                   utmpy(ga(2):gb(2)),bcmaty(2*mpidims(2),2*mpidims(2)),&
                   cndiagy(ra(2):rb(2),1:2)
  real(kind=pr) :: vlz(ra(3):rb(3)),vrz(ra(3):rb(3)),&
                   utmpz(ga(3):gb(3)),bcmatz(2*mpidims(1),2*mpidims(1)),&
                   cndiagz(ra(3):rb(3),1:2)
  real(kind=pr) :: vlx(ra(1):rb(1)),vrx(ra(1):rb(1)),&
                   utmpx(ga(1):gb(1)),cndiagx(ra(1):rb(1),1:2)

  ! This subroutine assumes that the domain decomposition is 2D
  ! The domain is NOT split in x direction
  ! In y direction, it is split in mpidims(2) parts
  ! In z direction, it is split in mpidims(1) parts

  ! Inverse of time step
  dtinv = 1.0d0/dt

  ! Half the inverse grid step
  dxinv = 0.5d0/dx
  dyinv = 0.5d0/dy
  dzinv = 0.5d0/dz

  ! Compute RHS
  ! Divergence of the velocity
  if (nx==1) then
    ! 2D case
    ix=ra(1)
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        work(ix,iy,iz,1) = (u(ix,iy+1,iz,2) - u(ix,iy-1,iz,2))*dyinv +&
                           (u(ix,iy,iz+1,3) - u(ix,iy,iz-1,3))*dzinv
      enddo
    enddo
  else
    ! 3D case
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          work(ix,iy,iz,1) = (u(ix+1,iy,iz,1) - u(ix-1,iy,iz,1))*dxinv +&
                             (u(ix,iy+1,iz,2) - u(ix,iy-1,iz,2))*dyinv +&
                             (u(ix,iy,iz+1,3) - u(ix,iy,iz-1,3))*dzinv
        enddo
      enddo
    enddo
  endif
  ! Multiply by -1/dt
  work(:,:,:,1) = -work(:,:,:,1)*dtinv
  ! Synchronize ghost points for the velocity
  call synchronize_ghosts(work,1)

  ! X DIRECTION
  ! Only if 3d
  if (nx.ne.1) then
    ! Set up parameters in x direction
    h2inv =  1.d0/(dx*dx)
    radir = ra(1)
    rbdir = rb(1)
    gadir = ga(1)
    gbdir = gb(1)
    ! Serial 1d solver init
    call gm_1d_serial_init (h2inv,radir,rbdir,&
                            det0,cndiagx,vlx,vrx)
    ! Loop for all lines x=const. This is local.
    do iz=ga(3),gb(3)
      !zz = dble(iz)*dz
      do iy=ga(2),gb(2)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpx(:) = work(gadir:gbdir,iy,iz,1)
        ! Solve linear system
        call gm_1d_serial_solver (radir,rbdir,gadir,gbdir,&
                                  det0,cndiagx,vlx,vrx,utmpx)
        ! Vector returned
        work(radir:rbdir,iy,iz,1) = utmpx(radir:rbdir)
      enddo    
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_x_serial (work(:,:,:,:),1)
  endif

  ! Y DIRECTION
  ! Set up parameters in y direction
  mpicommdir = mpicommy
  mpiszdir = mpidims(2)
  h2inv =  1.d0/(dy*dy)
  radir = ra(2)
  rbdir = rb(2)
  gadir = ga(2)
  gbdir = gb(2)
  ! Cases if # subdomains = 1 or >=2
  if (mpiszdir>1) then 
    ! Parallel 1d solver init
    call gm_1d_mpi_init (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,&
                         bcmaty,cndiagy,bcipivy,vly,vry)
    ! Loop for all lines y=const
    do iz = ga(3),gb(3)
      !zz = dble(iz)*dz
      do ix = ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpy(:) = work(ix,gadir:gbdir,iz,1)
        ! Solve linear system
        call gm_1d_mpi_solver (mpicommdir,mpiszdir,mpirankdir,radir,rbdir,gadir,gbdir,&
                               bcmaty,cndiagy,bcipivy,vly,vry,utmpy)
        ! Vector returned
        work(ix,radir:rbdir,iz,1) = utmpy(radir:rbdir)
      enddo
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_y_mpi (work(:,:,:,:),1)
  else
    ! Serial 1d solver init
    call gm_1d_serial_init (h2inv,radir,rbdir,&
                            det0,cndiagy,vly,vry)
    ! Loop for all lines y=const
    do iz=ga(3),gb(3)
      !zz = dble(iz)*dz
      do ix=ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpy(:) = work(ix,gadir:gbdir,iz,1)
        ! Solve linear system
        call gm_1d_serial_solver (radir,rbdir,gadir,gbdir,&
                                  det0,cndiagy,vly,vry,utmpy)
        ! Vector returned
        work(ix,radir:rbdir,iz,1) = utmpy(radir:rbdir)
      enddo    
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_y_serial (work(:,:,:,:),1)
  endif

  ! Z DIRECTION
  ! Set up parameters in z direction
  mpicommdir = mpicommz
  mpiszdir = mpidims(1)
  h2inv =  1.d0/(dz*dz)
  radir = ra(3)
  rbdir = rb(3)
  gadir = ga(3)
  gbdir = gb(3)
  ! Cases if # subdomains = 1 or >=2
  if (mpiszdir>1) then 
    ! Parallel 1d solver init
    call gm_1d_mpi_init (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,&
                         bcmatz,cndiagz,bcipivz,vlz,vrz)
    ! Loop for all lines z=const
    do iy=ga(2),gb(2)
      !yy = dble(iy)*dy
      do ix=ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpz(:) = work(ix,iy,gadir:gbdir,1)
        call gm_1d_mpi_solver (mpicommdir,mpiszdir,mpirankdir,radir,rbdir,gadir,gbdir,&
                               bcmatz,cndiagz,bcipivz,vlz,vrz,utmpz)
        ! Vector returned
        work(ix,iy,radir:rbdir,1) = utmpz(radir:rbdir)
      enddo
    enddo
  else
    ! Serial 1d solver init
    call gm_1d_serial_init (h2inv,radir,rbdir,&
                            det0,cndiagz,vlz,vrz)
    ! Loop for all lines z=const
    do iy=ga(2),gb(2)
      !yy = dble(iy)*dy
      do ix=ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpz(:) = work(ix,iy,gadir:gbdir,1)
        ! Solve linear system
        call gm_1d_serial_solver (radir,rbdir,gadir,gbdir,&
                                  det0,cndiagz,vlz,vrz,utmpz)
        ! Vector returned
        work(ix,iy,radir:rbdir,1) = utmpz(radir:rbdir)
      enddo    
    enddo
  endif

  ! Return pressure field
  fld(:,:,:) = work(:,:,:,1)

  ! Synchronize ghost points
  call synchronize_ghosts (fld(:,:,:))

end subroutine GM_pressure


!-------------------------------------------------------------------------------
! Guermond-Minev fractional step. Initialization of the 1d implicit MPI solver
!-------------------------------------------------------------------------------
subroutine gm_1d_mpi_init(mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,&
                          bcmat,cndiag,bcipiv,vl,vr)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: mpicommdir,mpiszdir,mpirankdir,radir,rbdir
  integer,intent(inout) :: bcipiv(2*mpiszdir)
  real(kind=pr),intent(inout) :: h2inv
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),bcmat(2*mpiszdir,2*mpiszdir),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! Local variables
  integer :: nn,j,mpicode
  real(kind=pr) :: sendfoo(4),recvfoo(4*mpiszdir)
  real(kind=pr) :: vl1(mpiszdir),vlN(mpiszdir),vr1(mpiszdir),vrN(mpiszdir)
  real(kind=pr) :: rhs(radir:rbdir)
  ! Get local ranks in the line
  call MPI_COMM_RANK(mpicommdir,mpirankdir,mpicode)
  ! Crank-Nicolson matrix in x direction
  cndiag(:,:) = 0.d0
  do j = radir,rbdir
    cndiag(j,1) = 1.d0 + 2.d0*h2inv
  enddo
  do j = radir,rbdir-1
    cndiag(j,2) = - h2inv
  enddo
  ! Factorize the CN matrix
  nn = rbdir-radir+1
  call factorize_tri_loc1d (cndiag,nn)
  ! Boundary conditions for domain decomposition
  ! BC influence basis
  rhs(:) = 0.d0
  rhs(radir) = 1.d0
  call solve_tri_loc1d (cndiag,rhs,vl,nn)
  vl(:) = -h2inv*vl(:)
  rhs(rbdir) = 1.d0
  rhs(radir) = 0.d0
  call solve_tri_loc1d (cndiag,rhs,vr,nn)
  vr(:) = -h2inv*vr(:)
  ! BC influence matrix
  ! It is only stored by one process
  ! Communicate values at the interface to rank 0
  sendfoo(1) = vl(radir)
  sendfoo(2) = vl(rbdir)
  sendfoo(3) = vr(radir)
  sendfoo(4) = vr(rbdir)
  call MPI_GATHER (sendfoo,4,MPI_DOUBLE_PRECISION,recvfoo,4,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  do j = 1,mpiszdir
    vl1(j) = recvfoo(4*j-3)
    vlN(j) = recvfoo(4*j-2)
    vr1(j) = recvfoo(4*j-1)
    vrN(j) = recvfoo(4*j)
  enddo
  ! BC influence matrix is only stored by one process
  if (mpirankdir == 0) then
    bcmat(:,:) = 0.d0
    do j = 1,mpiszdir
        bcmat(2*j-1,2*j-1) = 1.d0
        bcmat(2*j-1,modulo(2*j-3,2*mpiszdir)+1) = vl1(j)
        bcmat(2*j-1,modulo(2*j,2*mpiszdir)+1) = vr1(j)
        bcmat(2*j,2*j) = 1.d0
        bcmat(2*j,modulo(2*j-3,2*mpiszdir)+1) = vlN(j)
        bcmat(2*j,modulo(2*j,2*mpiszdir)+1) = vrN(j)
    enddo   
    ! Factorize the BC influence matrix
    call factorize_loc1d (bcmat,bcipiv,2*mpiszdir)
  endif
end subroutine gm_1d_mpi_init


!-------------------------------------------------------------------------------
! Guermond-Minev fractional step. 1d MPI solver
!-------------------------------------------------------------------------------
subroutine gm_1d_mpi_solver(mpicommdir,mpiszdir,mpirankdir,radir,rbdir,gadir,gbdir,&
                            bcmat,cndiag,bcipiv,vl,vr,utmp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: mpicommdir,mpiszdir,mpirankdir,radir,rbdir,gadir,gbdir
  integer,intent(inout) :: bcipiv(2*mpiszdir)
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),bcmat(2*mpiszdir,2*mpiszdir),&
                                 vl(radir:rbdir),vr(radir:rbdir),utmp(gadir:gbdir)
  ! local variables
  integer :: j,mpicode
  real(kind=pr) :: bcxl,bcxr,shortfoo(2*mpiszdir),longfoo(4*mpiszdir)
  real(kind=pr) :: bcrhs(2*mpiszdir),bcx(2*mpiszdir),&
                   rhs(radir:rbdir),vf(radir:rbdir),bcxls(mpiszdir),&
                   bcxrs(mpiszdir),vf1(mpiszdir),vfN(mpiszdir)
  
  ! Set RHS
  rhs(:) = utmp(radir:rbdir)
  ! Solve local system
  call solve_tri_loc1d (cndiag,rhs,vf,rbdir-radir+1)
  ! Communicate rhs to rank 0 in the line
  shortfoo(1) = vf(radir)
  shortfoo(2) = vf(rbdir)
  call MPI_GATHER (shortfoo,2,MPI_DOUBLE_PRECISION,longfoo,2,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  do j = 1,mpiszdir
    vf1(j) = longfoo(2*j-1)
    vfN(j) = longfoo(2*j)
  enddo
  ! BC influence RHS
  if (mpirankdir == 0) then
    do j = 1,mpiszdir
      bcrhs(2*j-1) = vf1(j)
      bcrhs(2*j) = vfN(j)
    enddo
    ! Solve BC influence system
    call solve_loc1d (bcmat,bcipiv,bcrhs,bcx,2*mpiszdir)
    ! Rearrange for mpi scatter
    do j = 1,mpiszdir
      bcxls(j) = bcx(modulo(2*j-3,2*mpiszdir)+1)
      bcxrs(j) = bcx(modulo(2*j,2*mpiszdir)+1)
    enddo
  endif
  ! Scatter from rank 0 in the line to all ranks
  do j = 1,mpiszdir
    longfoo(2*j-1) = bcxls(j)
    longfoo(2*j) = bcxrs(j)
  enddo
  call MPI_SCATTER (longfoo,2,MPI_DOUBLE_PRECISION,shortfoo,2,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  bcxl = shortfoo(1)
  bcxr = shortfoo(2)
  ! Superpose local solution and BC influence
  utmp(radir:rbdir) = vf(:)-bcxl*vl(:)-bcxr*vr(:)
end subroutine gm_1d_mpi_solver


!-------------------------------------------------------------------------------
! Guermond-Minev fractional step. 1d serial solver initialization
!-------------------------------------------------------------------------------
subroutine gm_1d_serial_init(h2inv,radir,rbdir,&
                             det0,cndiag,vl,vr)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: radir,rbdir
  real(kind=pr),intent(in) :: h2inv
  real(kind=pr),intent(out) :: det0
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! local variables
  integer :: j,nn
  real(kind=pr) :: rhs(radir:rbdir)

  ! Crank-Nicolson matrix in x direction
  cndiag(:,:) = 0.d0
  do j = radir,rbdir
    cndiag(j,1) = 1.d0 + 2.d0*h2inv
  enddo
  do j = radir,rbdir-1
    cndiag(j,2) = - h2inv
  enddo
  ! Factorize the CN matrix
  nn = rbdir-radir+1
  call factorize_tri_loc1d (cndiag,nn)
  ! Boundary conditions for domain decomposition
  ! BC influence basis
  rhs(:) = 0.d0
  rhs(radir) = 1.d0
  call solve_tri_loc1d (cndiag,rhs,vl,nn)
  vl(:) = -h2inv*vl(:)
  rhs(rbdir) = 1.d0
  rhs(radir) = 0.d0
  call solve_tri_loc1d (cndiag,rhs,vr,nn)
  vr(:) = -h2inv*vr(:)
  ! Compute determinant of the BC matrix
  det0 = (1.0d0+vr(radir))*(1.0d0+vl(rbdir))-vr(rbdir)*vl(radir)
end subroutine gm_1d_serial_init


!-------------------------------------------------------------------------------
! Guermond-Minev fractional step. 1d serial solver
!-------------------------------------------------------------------------------
subroutine gm_1d_serial_solver(radir,rbdir,gadir,gbdir,&
                               det0,cndiag,vl,vr,utmp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: radir,rbdir,gadir,gbdir
  real(kind=pr),intent(in) :: det0
  real(kind=pr),intent(inout) :: cndiag(radir:rbdir,1:2),utmp(gadir:gbdir),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! local variables
  real(kind=pr) :: detr,detl
  real(kind=pr) :: rhs(radir:rbdir),vf(radir:rbdir)

  ! Set RHS
  rhs(:) = utmp(radir:rbdir)
  ! Solve local system
  call solve_tri_loc1d (cndiag,rhs,vf,rbdir-radir+1)
  ! Solve boundary correction
  detr = vf(radir)*(1.0d0+vl(rbdir))-vf(rbdir)*vl(radir)
  detl = vf(rbdir)*(1.0d0+vr(radir))-vf(radir)*vr(rbdir)
  ! Corrected vector
  utmp(radir:rbdir) = vf - (detl*vl+detr*vr)/det0
end subroutine gm_1d_serial_solver

