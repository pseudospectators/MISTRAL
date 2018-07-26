!-------------------------------------------------------------------------------
! Warpper calling different individual time steppers
!-------------------------------------------------------------------------------
subroutine FluidTimestep(time,dt,u,nlk,work,mask,mask_color,us,Insect)
    use vars
    use insect_module
    implicit none

    real(kind=pr), intent(inout) :: time, dt
    real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
    real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
    real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
    real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
    real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
    integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

    type(diptera), intent(inout) :: Insect

    real(kind=pr)::t1
    t1=MPI_wtime()

    !-----------------------------------------------------------------------------
    ! Call fluid advancement subroutines.
    !-----------------------------------------------------------------------------
    select case(iTimeMethodFluid)
    case("RK4")
        call RK4(time,dt,u,nlk,work,mask,mask_color,us,Insect)

    case default
        call abort(1234,"Error! iTimeMethodFluid unknown. Abort.")

    end select

    ! Force zero mode for mean flow
    call set_mean_flow(u,time)

    !-----------------------------------------------------------------------------
    ! compute unsteady corrections in every time step.
    !-----------------------------------------------------------------------------
    if (unst_corrections==1) then
        call cal_unst_corrections ( time, mask, mask_color, us, Insect )
    endif

    time_fluid=time_fluid + MPI_wtime() - t1
end subroutine FluidTimestep

!-------------------------------------------------------------------------------

subroutine RK4(time,dt,u,nlk,work,mask,mask_color,us,Insect)
    use vars
    use insect_module
    use basic_operators
    implicit none

    real(kind=pr), intent(inout) :: time, dt
    real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
    real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
    real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
    real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
    real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
    integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

    type(diptera), intent(inout) :: Insect
    real(kind=pr) :: t
    t = time

    call adjust_dt(time,u,dt)

    ! NLK 5th register holds old velocity
    nlk(:,:,:,:,5) = u


    call cal_nlk(time,u,nlk(:,:,:,:,1),work,mask,mask_color,us,Insect)

    u = nlk(:,:,:,:,5) + 0.5d0*dt*nlk(:,:,:,:,1)
    t = time + 0.5d0*dt
    call cal_nlk(t,u,nlk(:,:,:,:,2),work,mask,mask_color,us,Insect)

    u = nlk(:,:,:,:,5) + 0.5d0*dt*nlk(:,:,:,:,2)
    t = time + 0.5d0*dt
    call cal_nlk(t,u,nlk(:,:,:,:,3),work,mask,mask_color,us,Insect)

    u = nlk(:,:,:,:,5) + dt * nlk(:,:,:,:,3)
    t = time + dt
    call cal_nlk(t,u,nlk(:,:,:,:,4),work,mask,mask_color,us,Insect)

    u = nlk(:,:,:,:,5) + dt/6.d0*(nlk(:,:,:,:,1)+2.d0*nlk(:,:,:,:,2)&
    +2.d0*nlk(:,:,:,:,3)+nlk(:,:,:,:,4))

end subroutine RK4


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

    real(kind=pr), intent(in) :: time
    real(kind=pr), intent(inout) :: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
    real(kind=pr), intent(out) :: dt1
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
