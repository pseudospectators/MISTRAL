!-------------------------------------------------------------------------------
! add a cavity mask around the domain
! note: if the obstacle extends in this zone (e.g. the flying insect touching
! the lower end of the domain), it will be overwritten.
! the cavity walls have the color "0" which helps excluding them
! for example when computing the integral forces
!-------------------------------------------------------------------------------
subroutine add_cavity(mask,mask_color,us)
    use vars
    implicit none

    real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
    real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
    integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

    integer :: ix,iy,iz
    real(kind=pr) :: ux, uy, uz, x, y, z, tmp

    ux = 0.d0
    uy = 0.d0
    uz = 0.d0

    if (iCavity=="yes") then
        ! this is the normal cavity forcing homogeneous dirichlet BC
        ux = 0.d0
        uy = 0.d0
        uz = 0.d0
    elseif (iCavity=="freestream".or. iCavity=="freestream_smooth") then
        ! cavity to force free stream inflow
        ux = uxmean
        uy = uymean
        uz = uzmean
    endif

    if (iCavity=="freestream_smooth") then
        ! this is a smooth cavity, with a smoothing layer that is not small, but extends
        ! over half the width "thick_wall" of the cavity. As for the other freestream cavity,
        ! it is not really a cavity in the sense that it is surronded by solid walls, but rather
        ! enforces the free-flow conditions "at infinity" which is just mean flow. Its particular
        ! smooth shape is for comparison with the artificial compressibility method, where we use
        ! this to minimize the reflection of pressure waves.
        if (nx>1) then
            !-----3D-----
            do iz = ra(3), rb(3)
                z = dble(iz) * zl / dble(nz)
                do iy = ra(2), rb(2)
                    y = dble(iy) * yl / dble(ny)
                    do ix = ra(1), rb(1)
                        x = dble(ix) * xl / dble(nx)
                        ! distance to borders of domain
                        tmp = minval( (/x,y,z,-(x-xl),-(y-yl),-(z-zl)/) )
                        if (tmp < thick_wall) then
                            call smoothstep(mask(ix,iy,iz), tmp, 0.5_pr*thick_wall, 0.5_pr*thick_wall)
                            us(ix,iy,iz,1:3) = (/ux,uy,uz/)
                            ! sponges have color -1 (decision made in rhs.f90 if mask or sponge)
                            mask_color(ix,iy,iz) = -1
                        endif
                    enddo
                enddo
            enddo
        else
            !-----2D-----
            do iz = ra(3), rb(3)
                z = dble(iz) * zl / dble(nz)
                do iy = ra(2), rb(2)
                    y = dble(iy) * yl / dble(ny)
                    ! distance to borders of domain
                    tmp = minval( (/z,y,-(z-zl),-(y-yl)/) )
                    if (tmp < thick_wall) then
                        call smoothstep(mask(0,iy,iz), tmp, 0.5_pr*thick_wall, 0.5_pr*thick_wall)
                        us(0,iy,iz,1:3) = (/ux,uy,uz/)
                        ! sponges have color -1 (decision made in rhs.f90 if mask or sponge)
                        mask_color(0,iy,iz) = -1
                    endif

                enddo
            enddo
        endif

    elseif (iCavity=="freestream_smooth_only_pressure") then
        ! as above, but acts ONLY on the pressure (forcing it to zero)
        do iz = ra(3), rb(3)
            z = dble(iz) * zl / dble(nz)
            do iy = ra(2), rb(2)
                y = dble(iy) * yl / dble(ny)
                do ix = ra(1), rb(1)
                    x = dble(ix) * xl / dble(nx)
                    ! distance to borders of domain
                    tmp = minval( (/x,y,z,-(x-xl),-(y-yl),-(z-zl)/) )
                    if (tmp < thick_wall) then
                        call smoothstep(mask(ix,iy,iz), tmp, 0.5_pr*thick_wall, 0.5_pr*thick_wall)
                        us(ix, iy, iz, 4) = 0.0d0
                        ! sponges have color -1 (decision made in rhs.f90 if mask or sponge)
                        mask_color(ix,iy,iz) = -1
                    endif
                enddo
            enddo
        enddo

    else
        ! real cavity, solid walls
        ! non-smooth cavity for free stream or solid walls.
        do iz = ra(3), rb(3)
            do iy = ra(2), rb(2)
                do ix = ra(1), rb(1)
                    if (nx>1) then
                        if ((ix<=cavity_size-1).or.(ix>=nx-1-cavity_size+1)) then
                            mask(ix,iy,iz) = 1.d0
                            us(ix,iy,iz,:) = (/ux,uy,uz/)
                            ! external boxes have color 0 (important for forces)
                            mask_color(ix,iy,iz) = 0
                        endif
                    endif
                    if ((iy<=cavity_size-1).or.(iy>=ny-1-cavity_size+1)) then
                        mask(ix,iy,iz) = 1.d0
                        us(ix,iy,iz,:) = (/ux,uy,uz/)
                        ! external boxes have color 0 (important for forces)
                        mask_color(ix,iy,iz) = 0
                    endif

                    if ((iz<=cavity_size-1).or.(iz>=nz-1-cavity_size+1)) then
                        mask(ix,iy,iz) = 1.d0
                        us(ix,iy,iz,:) = (/ux,uy,uz/)
                        ! external boxes have color 0 (important for forces)
                        mask_color(ix,iy,iz) = 0
                    endif
                enddo
            enddo
        enddo
    endif

end subroutine Add_Cavity
