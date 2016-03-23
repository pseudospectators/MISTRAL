!====================================================================
!====================================================================
!
!     This module is derived from the P3DFFT interface in FLUSI
!     In MISTRAL it only performs the domain decomposition
!     No FFT installation is required
!
!====================================================================
!====================================================================

module p3dfft_wrapper
  use vars
  implicit none
  
  contains

!----------------------------------------------------------------
! set up communicators for y and z directions separately
!----------------------------------------------------------------
subroutine setup_cart_groups
  ! Setup 1d communicators
  use mpi
  use vars

  implicit none
  integer :: mpicolor,mpikey,mpicode
  integer :: mpicommtmp1,mpicommtmp2
  logical :: mpiperiods(2),periods(1),reorder
  integer :: one=1,two=2,dims(1) ! Required for MPI_CART_GET, MPI_CART_CREATE in openmpi
  ! Set parameters
  periods(1)=.true. ! This should be an array - if not, openmpi fails
  reorder=.false.
  ! Get Cartesian topology information
  call MPI_CART_GET(mpicommcart,two,mpidims,mpiperiods,mpicoords,mpicode)
  ! Communicator for line in y direction
  mpicolor = mpicoords(2) 
  mpikey = mpicoords(1)
  call MPI_COMM_SPLIT (mpicommcart,mpicolor,mpikey,mpicommtmp1,mpicode)
  dims(1) = mpidims(1)
  call MPI_CART_CREATE(mpicommtmp1,one,dims,periods,reorder,mpicommz,mpicode)
  ! Communicator for line in z direction
  mpicolor = mpicoords(1) 
  mpikey = mpicoords(2)
  call MPI_COMM_SPLIT (mpicommcart,mpicolor,mpikey,mpicommtmp2,mpicode)
  dims(1) = mpidims(2)
  call MPI_CART_CREATE(mpicommtmp2,one,dims,periods,reorder,mpicommy,mpicode)
end subroutine setup_cart_groups


!-----------------------------------------------------------------------------
!     Initialize domain decomposition
!-----------------------------------------------------------------------------
subroutine decomposition_initialize
  use mpi ! Module incapsulates mpif.
  use vars
  implicit none

  integer,parameter :: nmpidims = 2
  integer :: mpicode,idir,L,n
  integer,dimension(1:3) :: ka,kb,ks,kat,kbt,kst
  logical,dimension(2) :: subcart
  real(kind=pr),dimension(:,:),allocatable :: f,ft

  ! default case, no decomposition
  decomposition="none"
  
  !-- Set up dimensions. It is very important that mpidims(2) > mpidims(1)
  ! because P3Dfft crashes otherwise. This means a 1D decomposition is always
  ! along the z direction in real space. 
  if (nz>mpisize.and.nx==1) then
     mpidims(1) = 1             ! due to p3dfft, 1D decomposition is always the
     mpidims(2) = mpisize       ! 3rd index in real space.
     decomposition="1D"
  else
     ! unfortunately, 2D data decomposition does not work with 2D code (nx==1)
     mpidims = 0
     call MPI_Dims_create(mpisize,nmpidims,mpidims,mpicode)
     if(mpidims(1) > mpidims(2)) then
        mpidims(1) = mpidims(2)
        mpidims(2) = mpisize / mpidims(1)
     endif
     decomposition="2D"
  endif
  
  if (root) write(*,'("mpidims= ",i3,1x,i3)') mpidims
  if (root) write(*,'("Using ",A," decomposition!")') trim(adjustl(decomposition))

  !-- Check dimensions
  if(mpidims(1)*mpidims(2)/=mpisize) then
     print *, 'wrong mpidims: change mpisize'
     call abort()
  endif

  !-- Set subdomain bounds
  !-- Get Cartesian topology info
  !-- Get local sizes
  call p3dfft_stub(mpidims,nx,ny,nz,MPI_COMM_WORLD,mpitaskid,mpitasks,mpicommcart,ra,rb,rs) 
  ra(:) = ra(:) - 1
  rb(:) = rb(:) - 1
  
  !-- extents of real arrays that have ghost points. We add ghosts in all 
  !-- directions, including the periodic ones.
  ga=ra-ng
  gb=rb+ng
  
  if (nx==1) then
    ga(1)=0
    gb(1)=0
  endif
  
  if ( rb(2)-ra(2)+1<2*ng .or. rb(3)-ra(3)+1<2*ng ) then
    if (mpirank==0) write(*,*) "Too many CPUs: the ghosts span more than one CPU"
    if (mpirank==0) write(*,*) "y", rb(2)-ra(2)+1, "z", rb(3)-ra(3)+1
    call abort()
  endif
  
  !-- Allocate domain partitioning tables and gather sizes from all processes 
  !-- (only for real arrays)
  ! TODO: These tables are currently not used for communication between subdomains,
  ! but may be still useful for development/debugging purposes.
  allocate ( ra_table(1:3,0:mpisize-1), rb_table(1:3,0:mpisize-1) )
  call MPI_ALLGATHER (ra, 3, MPI_INTEGER, ra_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)
  call MPI_ALLGATHER (rb, 3, MPI_INTEGER, rb_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)

end subroutine decomposition_initialize


!----------------------------------------------------------------
! wavenumber functions: return the kx,ky,kz wavenumbers
! as a function of the array index
!----------------------------------------------------------------
real(kind=pr) function wave_x( ix )
  use vars ! for scale and precision statement
  implicit none
  integer, intent (in) :: ix
  wave_x = scalex*dble(ix)
end function

real(kind=pr) function wave_y( iy )
  use vars ! for scale and precision statement
  implicit none
  integer, intent (in) :: iy
  wave_y = scaley*dble(modulo(iy+ny/2,ny)-ny/2)
end function

real(kind=pr) function wave_z( iz )
  use vars ! for scale and precision statement
  implicit none
  integer, intent (in) :: iz
  wave_z = scalez*dble(modulo(iz+nz/2,nz)-nz/2)
end function


!----------------------------------------------------------------
! domain decomposition routine derived from P3DFFT
!----------------------------------------------------------------
subroutine p3dfft_stub(dims_in,nx,ny,nz,mpi_comm_in,mpi_taskid,mpi_tasks,mpi_comm_out,istart,iend,isize)
  implicit none

  integer, intent (in) :: nx,ny,nz,mpi_comm_in
  integer, intent (out) :: istart(3),iend(3),isize(3)
  integer, intent (out) :: mpi_taskid,mpi_tasks,mpi_comm_out
  integer, intent (in) :: dims_in(2)

  integer :: i,j,k,ierr,mpicomm,mpi_comm_cart
  integer :: ipid,jpid,iproc,jproc
  integer :: numtasks,taskid
  integer :: cartid(2),dims(2)
  integer, dimension (:), allocatable :: jist,jisz,jien,kjst,kjsz,kjen
  logical :: periodic(2),remain_dims(2)

  if(nx .le. 0 .or. ny .le. 0 .or. nz .le. 0) then
     print *,'Invalid dimensions :',nx,ny,nz
     call abort()
  endif

  mpicomm = mpi_comm_in
  call MPI_COMM_SIZE (mpicomm,numtasks,ierr)
  call MPI_COMM_RANK (mpicomm,taskid,ierr)

  if(dims_in(1) .le. 0 .or. dims_in(2) .le. 0 .or.  dims_in(1)*dims_in(2) .ne. numtasks) then
     print *,'Invalid processor geometry: ',dims,' for ',numtasks, 'tasks'
     call abort()
  endif

  if(taskid .eq. 0) then 
     print *,'Using stride-1 layout'
  endif

  iproc = dims_in(1)
  jproc = dims_in(2)
  dims(1) = dims_in(2)
  dims(2) = dims_in(1)

  periodic(1) = .false.
  periodic(2) = .false.
! creating cartesian processor grid
  call MPI_Cart_create(mpicomm,2,dims,periodic,.false.,mpi_comm_cart,ierr)
! Obtaining process ids with in the cartesian grid
  call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)
! process with a linear id of 5 may have cartid of (3,1)

  ipid = cartid(2)
  jpid = cartid(1)

  allocate (jist(0:iproc-1))
  allocate (jisz(0:iproc-1))
  allocate (jien(0:iproc-1))
  allocate (kjst(0:jproc-1))
  allocate (kjsz(0:jproc-1))
  allocate (kjen(0:jproc-1))
!
!Mapping 3-D data arrays onto 2-D process grid
! (nx+2,ny,nz) => (iproc,jproc)      
! 
  call MapDataToProc(ny,iproc,jist,jien,jisz)
  call MapDataToProc(nz,jproc,kjst,kjen,kjsz)

! These are local array indices for each processor
  istart(1) = 1
  iend(1) = nx
  isize(1) = nx
  istart(2) = jist(ipid)
  iend(2) = jien(ipid)
  isize(2) = jisz(ipid)
  istart(3) = kjst(jpid)
  iend(3) = kjen(jpid)
  isize(3) = kjsz(jpid)

  deallocate(jist,jisz,jien,kjst,kjsz,kjen)

  mpi_taskid = taskid
  mpi_tasks = numtasks
  mpi_comm_out = mpi_comm_cart

end subroutine p3dfft_stub


!----------------------------------------------------------------
! calculate subdomain bounds
!----------------------------------------------------------------
subroutine MapDataToProc (data,proc,st,en,sz)
  implicit none
  integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
  integer i,size,nl,nu

  size=data/proc
  nu = data - size * proc
  nl = proc - nu
  st(0) = 1
  sz(0) = size
  en(0) = size
  do i=1,nl-1
     st(i) = st(i-1) + size
     sz(i) = size
     en(i) = en(i-1) + size
  enddo
  size = size + 1
  do i=nl,proc-1
     st(i) = en(i-1) + 1
     sz(i) = size
     en(i) = en(i-1) + size
  enddo
  en(proc-1)= data 
  sz(proc-1)= data-st(proc-1)+1
end subroutine

end module p3dfft_wrapper




