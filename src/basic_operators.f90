!-------------------------------------------------------------------------------
! this module contains basic elementary operators, and as it is a module
! function overloading can be used.
! List of functions:
!       * curl
!       * divergence
!       * max/min operations on real fields
!       * check fields for NaN
!-------------------------------------------------------------------------------
module basic_operators

  use ghosts
  use vars

  !-- interface for curl operators
  interface curl
    module procedure curl_x
  end interface

  !-- interface for maxima of fields
  interface fieldmaxabs
    module procedure fieldmaxabs, fieldmaxabs3
  end interface

  !-- check fields for NaN
  interface checknan
    module procedure checknan_cmplx, checknan_real
  end interface

  ! divergence (in x- or k- space)
  interface divergence
    module procedure divergence_x
  end interface


!!!!!!!!!!!
 contains
!!!!!!!!!!!


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl_x( u, rotu )
  use p3dfft_wrapper
  implicit none

  ! input/output field in x-space
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout)::rotu(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,dxinv,dyinv,dzinv,a1,a2,a4,a5
  real(kind=pr) :: uxdy,uxdz,uydx,uydz,uzdx,uzdy
  real(kind=pr)::a(-3:+3)


  select case(method)
  case('centered_2nd')
      call synchronize_ghosts( u, 3 )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
      !-------------------------------------------------------------------------
      dxinv = 1.d0/(2.d0*dx)
      dyinv = 1.d0/(2.d0*dy)
      dzinv = 1.d0/(2.d0*dz)

      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dyinv
              uxdz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dzinv
              uydx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dxinv
              uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
              uzdx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dxinv
              uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv

              rotu(ix,iy,iz,1) = uzdy - uydz
              rotu(ix,iy,iz,2) = uxdz - uzdx
              rotu(ix,iy,iz,3) = uydx - uxdy
            enddo
          enddo
        enddo
      elseif (nx==1) then
        ! two-dimensional simulation
        ix=0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
            uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv

            rotu(ix,iy,iz,1) = uzdy - uydz
            rotu(ix,iy,iz,2) = 0.d0
            rotu(ix,iy,iz,3) = 0.d0
          enddo
        enddo
      endif

  case('centered_4th')
      call synchronize_ghosts( u, 3 )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
      !-------------------------------------------------------------------------
      ! these coefficients are the th order optimized scheme by Tamm&Webb
      a=(/-0.02651995d0, +0.18941314d0, -0.79926643d0, 0.0d0, &
           0.79926643d0, -0.18941314d0, 0.02651995d0/)

      dxinv = 1.d0/dx
      dyinv = 1.d0/dy
      dzinv = 1.d0/dz

      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdy = (a(-3)*u(ix,iy-3,iz,1)+a(-2)*u(ix,iy-2,iz,1)+a(-1)*u(ix,iy-1,iz,1)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix,iy+3,iz,1)+a(+2)*u(ix,iy+2,iz,1)+a(+1)*u(ix,iy+1,iz,1))*dyinv
              uxdz = (a(-3)*u(ix,iy,iz-3,1)+a(-2)*u(ix,iy,iz-2,1)+a(-1)*u(ix,iy,iz-1,1)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix,iy,iz+3,1)+a(+2)*u(ix,iy,iz+2,1)+a(+1)*u(ix,iy,iz+1,1))*dzinv

              uydx = (a(-3)*u(ix-3,iy,iz,2)+a(-2)*u(ix-2,iy,iz,2)+a(-1)*u(ix-1,iy,iz,2)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix+3,iy,iz,2)+a(+2)*u(ix+2,iy,iz,2)+a(+1)*u(ix+1,iy,iz,2))*dxinv
              uydz = (a(-3)*u(ix,iy,iz-3,2)+a(-2)*u(ix,iy,iz-2,2)+a(-1)*u(ix,iy,iz-1,2)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix,iy,iz+3,2)+a(+2)*u(ix,iy,iz+2,2)+a(+1)*u(ix,iy,iz+1,2))*dzinv

              uzdx = (a(-3)*u(ix-3,iy,iz,3)+a(-2)*u(ix-2,iy,iz,3)+a(-1)*u(ix-1,iy,iz,3)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix+3,iy,iz,3)+a(+2)*u(ix+2,iy,iz,3)+a(+1)*u(ix+1,iy,iz,3))*dxinv
              uzdy = (a(-3)*u(ix,iy-3,iz,3)+a(-2)*u(ix,iy-2,iz,3)+a(-1)*u(ix,iy-1,iz,3)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix,iy+3,iz,3)+a(+2)*u(ix,iy+2,iz,3)+a(+1)*u(ix,iy+1,iz,3))*dyinv


              rotu(ix,iy,iz,1) = uzdy - uydz
              rotu(ix,iy,iz,2) = uxdz - uzdx
              rotu(ix,iy,iz,3) = uydx - uxdy
            enddo
          enddo
        enddo
      elseif (nx==1) then
        ! two-dimensional simulation
        ix=0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)

            uydz = (a(-3)*u(ix,iy,iz-3,2)+a(-2)*u(ix,iy,iz-2,2)+a(-1)*u(ix,iy,iz-1,2)+a(0)*u(ix,iy,iz,1)&
                   +a(+3)*u(ix,iy,iz+3,2)+a(+2)*u(ix,iy,iz+2,2)+a(+1)*u(ix,iy,iz+1,2))*dzinv
            uzdy = (a(-3)*u(ix,iy-3,iz,3)+a(-2)*u(ix,iy-2,iz,3)+a(-1)*u(ix,iy-1,iz,3)+a(0)*u(ix,iy,iz,1)&
                   +a(+3)*u(ix,iy+3,iz,3)+a(+2)*u(ix,iy+2,iz,3)+a(+1)*u(ix,iy+1,iz,3))*dyinv

            rotu(ix,iy,iz,1) = uzdy - uydz
            rotu(ix,iy,iz,2) = 0.d0
            rotu(ix,iy,iz,3) = 0.d0
          enddo
        enddo
      endif

  case default
      call abort('invalid METHOD in curl_x:'//method)
  end select
end subroutine curl_x


!-------------------------------------------------------------------------------
! compute the divergence of the input field u and return it in divu
! depending on the value of "method", different discretization
! is used.
!-------------------------------------------------------------------------------
subroutine divergence_x( u, divu )
  use p3dfft_wrapper
  implicit none
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout)::divu(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  integer :: ix,iy,iz
  real(kind=pr) :: dxinv,dyinv,dzinv,uxdx,uydy,uzdz,a1,a2,a4,a5
  real(kind=pr)::a(-3:+3)

  select case(method)
  case('centered_2nd')
      call synchronize_ghosts( u, 3 )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
      !-------------------------------------------------------------------------
      dxinv = 1.d0/(2.d0*dx)
      dyinv = 1.d0/(2.d0*dy)
      dzinv = 1.d0/(2.d0*dz)

      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdx = dxinv*(u(ix+1,iy,iz,1)-u(ix-1,iy,iz,1))
              uydy = dyinv*(u(ix,iy+1,iz,2)-u(ix,iy-1,iz,2))
              uzdz = dzinv*(u(ix,iy,iz+1,3)-u(ix,iy,iz-1,3))
              divu(ix,iy,iz) = uxdx + uydy + uzdz
            enddo
          enddo
        enddo
      elseif (nx==1) then
        ! two-dimensional simulation
        ix=0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydy = dyinv*(u(ix,iy+1,iz,2)-u(ix,iy-1,iz,2))
            uzdz = dzinv*(u(ix,iy,iz+1,3)-u(ix,iy,iz-1,3))
            divu(ix,iy,iz) = uydy + uzdz
          enddo
        enddo
      endif

  case('centered_4th')
      call synchronize_ghosts( u, 3 )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
      !-------------------------------------------------------------------------
      ! these coefficients are the th order optimized scheme by Tamm&Webb
      a=(/-0.02651995d0, +0.18941314d0, -0.79926643d0, 0.0d0, &
           0.79926643d0, -0.18941314d0, 0.02651995d0/)

      dxinv = 1.d0/dx
      dyinv = 1.d0/dy
      dzinv = 1.d0/dz

      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdx = (a(-3)*u(ix-3,iy,iz,1)+a(-2)*u(ix-2,iy,iz,1)+a(-1)*u(ix-1,iy,iz,1)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix+3,iy,iz,1)+a(+2)*u(ix+2,iy,iz,1)+a(+1)*u(ix+1,iy,iz,1))*dxinv
              uydy = (a(-3)*u(ix,iy-3,iz,2)+a(-2)*u(ix,iy-2,iz,2)+a(-1)*u(ix,iy-1,iz,2)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix,iy+3,iz,2)+a(+2)*u(ix,iy+2,iz,2)+a(+1)*u(ix,iy+1,iz,2))*dyinv
              uzdz = (a(-3)*u(ix,iy,iz-3,3)+a(-2)*u(ix,iy,iz-2,3)+a(-1)*u(ix,iy,iz-1,3)+a(0)*u(ix,iy,iz,1)&
                     +a(+3)*u(ix,iy,iz+3,3)+a(+2)*u(ix,iy,iz+2,3)+a(+1)*u(ix,iy,iz+1,3))*dzinv

              divu(ix,iy,iz) = uxdx + uydy + uzdz
            enddo
          enddo
        enddo
      elseif (nx==1) then
        ! two-dimensional simulation
        ix = 0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydy = (a(-3)*u(ix,iy-3,iz,2)+a(-2)*u(ix,iy-2,iz,2)+a(-1)*u(ix,iy-1,iz,2)+a(0)*u(ix,iy,iz,1)&
                   +a(+3)*u(ix,iy+3,iz,2)+a(+2)*u(ix,iy+2,iz,2)+a(+1)*u(ix,iy+1,iz,2))*dyinv
            uzdz = (a(-3)*u(ix,iy,iz-3,3)+a(-2)*u(ix,iy,iz-2,3)+a(-1)*u(ix,iy,iz-1,3)+a(0)*u(ix,iy,iz,1)&
                   +a(+3)*u(ix,iy,iz+3,3)+a(+2)*u(ix,iy,iz+2,3)+a(+1)*u(ix,iy,iz+1,3))*dzinv

            divu(ix,iy,iz) = uydy + uzdz
          enddo
        enddo
      endif

  case default
    call abort('invalid METHOD in divergence:'//method)
  end select
end subroutine divergence_x



!-------------------------------------------------------------------------------
! returns the globally largest entry of a given (real) field
! only inner points considered
!-------------------------------------------------------------------------------
real(kind=pr) function fieldmax( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(inout):: inx(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr) :: max_local, max_global
  integer :: mpicode

  max_local = maxval(inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmax = max_global
end function fieldmax

!-------------------------------------------------------------------------------
! returns the globally smallest entry of a given (real) field
!-------------------------------------------------------------------------------
real(kind=pr) function fieldmin( inx )
  implicit none
  real(kind=pr),intent(inout):: inx(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr) :: min_local, min_global
  integer :: mpicode

  min_local = minval(inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call MPI_ALLREDUCE (min_local,min_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmin = min_global
end function fieldmin


!-------------------------------------------------------------------------------
! returns the globally largest entry of a given vector field
! (L2-norm)
!-------------------------------------------------------------------------------
real(kind=pr) function fieldmaxabs3( inx )
  implicit none
  real(kind=pr),intent(inout) :: inx(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr) :: max_local, max_global, value
  integer :: mpicode
  integer :: ix, iy, iz
  max_local = 0.d0
  do ix = ra(1), rb(1)
    do iy = ra(2), rb(2)
       do iz = ra(3), rb(3)
        value = inx(ix,iy,iz,1)*inx(ix,iy,iz,1) + inx(ix,iy,iz,2)*inx(ix,iy,iz,2) &
              + inx(ix,iy,iz,3)*inx(ix,iy,iz,3)
        if (max_local<value) max_local=value
       enddo
    enddo
  enddo

  max_local = dsqrt( max_local )
  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmaxabs3 = max_global
end function fieldmaxabs3


!-------------------------------------------------------------------------------
! returns the globally largest entry of a given vector field
! (L2-norm)
!-------------------------------------------------------------------------------
real(kind=pr) function fieldmaxabs( inx1, inx2, inx3 )
  implicit none
  real(kind=pr),intent(inout):: inx1(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout):: inx2(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout):: inx3(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr) :: max_local, max_global
  integer :: mpicode

  max_local = maxval( inx1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) * &
                      inx1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) + &
                      inx2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) * &
                      inx2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) + &
                      inx3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) * &
                      inx3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  max_local = dsqrt( max_local )

  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmaxabs = max_global
end function fieldmaxabs


!-------------------------------------------------------------------------------
! check a real valued field for NaNs and display warning if found
!-------------------------------------------------------------------------------
subroutine checknan_real( field, msg )
  implicit none
  real(kind=pr),intent(inout)::field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  character(len=*),intent(inout)::msg
  integer :: foundnan,foundnans,mpicode,ix,iy,iz
  foundnan = 0
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)
        if (is_nan(field(ix,iy,iz))) foundnan = 1
      enddo
    enddo
  enddo

  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)

  if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans
end subroutine checknan_real


!-------------------------------------------------------------------------------
! check a complex field for NaN's, display warning if found
!-------------------------------------------------------------------------------
subroutine checknan_cmplx( field, msg )
  implicit none
  complex(kind=pr),intent(inout)::field(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  character(len=*),intent(inout)::msg
  integer :: foundnan,foundnans,mpicode,ix,iy,iz
  foundnan = 0
  do iz=ca(1),cb(1)
    do iy=ca(2),cb(2)
      do ix=ca(3),cb(3)
        if (is_nan(aimag(field(iz,iy,ix)))) foundnan = 1
        if (is_nan(real (field(iz,iy,ix)))) foundnan = 1
      enddo
    enddo
  enddo

  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)

  if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans
end subroutine checknan_cmplx

!-------------------------------------------------------------------------------
! computes the volume integral of the scalar quantity u
!-------------------------------------------------------------------------------
real(kind=pr) function  volume_integral( u )
  use p3dfft_wrapper
  implicit none

  ! input/output field in x-space
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  integer::ix,iy,iz,mpicode
  real(kind=pr)::int_local,dxyz

  int_local = 0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        int_local = int_local + u(ix,iy,iz)
      enddo
    enddo
  enddo

  dxyz = dx*dy*dz
  int_local = int_local*dxyz

  call MPI_ALLREDUCE (int_local,volume_integral,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
end function volume_integral


end module basic_operators
