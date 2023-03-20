submodule(OceanMod) OceanExternalRoutines
  implicit none

  contains

  subroutine init_state_ocean_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                           :: i, j, m, jm_int, ndI1, jmsI, jmvI
    real(kind=dbl),    allocatable    :: r(:)
    complex(kind=dbl), allocatable    :: velc(:), temp(:,:), spher1(:,:), torr(:,:), spher2(:,:)

    if (.not. init_through_file_ocean) then
      do i = 1, this%nd+1
        do j = 0, this%jmax
          do m = 0, j
            jm_int = jm(j,m)

            if ((j == 0) .and. (m == 0)) then
              this%sol%temp(3*(i-1)+1,jm_int)%re = (this%rad_grid%r(this%nd)/this%rad_grid%rr(i)-1)*this%rad_grid%r(1)*sqrt(4*pi)
            else if (m == 0) then
              call random_number( this%sol%temp(3*(i-1)+1, jm_int)%re )
              this%sol%temp(3*(i-1)+1, jm_int)%re = this%sol%temp(3*(i-1)+1, jm_int)%re / 1e3
            else
              call random_number( this%sol%temp(3*(i-1)+1, jm_int)%re ); call random_number( this%sol%temp(3*(i-1)+1, jm_int)%im )
              this%sol%temp(3*(i-1)+1, jm_int) = this%sol%temp(3*(i-1)+1, jm_int) / 1e3
            end if

          end do
        end do
      end do

    else
      ndI1 = nd_init_ocean+1; jmsI = jm(jmax_init_ocean,jmax_init_ocean); jmvI = jml(jmax_init_ocean,jmax_init_ocean,+1)

      allocate( r(ndI1), velc(jmvI), temp(ndI1,jmsI), spher1(ndI1,jmsI), spher2(ndI1,jmsI), torr(ndI1,jmsI) )
        spher1 = cmplx(0._dbl, 0._dbl, kind=dbl); spher2 = cmplx(0._dbl, 0._dbl, kind=dbl); torr = cmplx(0._dbl, 0._dbl, kind=dbl)

        open(unit=8, file='code/ocean/inittemp', status='old', action='read')
          do i = 1, ndI1
            read(8,*) r(i), temp(i,:)
          end do
        close(8)

        open(unit=8, file='code/ocean/initvelc', status='old', action='read')
          do i = 1, ndI1
            read(8,*) r(i), velc

            do jm_int = 2, jmsI
              spher1(i,jm_int) = velc(3*(jm_int-1)-1)
              torr(  i,jm_int) = velc(3*(jm_int-1)  )
              spher2(i,jm_int) = velc(3*(jm_int-1)+1)
            end do
          end do
        close(8)

      deallocate(velc)

      do i = 1, this%nd+1
        this%sol%temp(3*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, temp  )
        this%sol%mech(6*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, spher1)
        this%sol%torr(3*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, torr  )
        this%sol%mech(6*(i-1)+2,:) = this%rad_grid%interpolation_fn(this%jms, i, r, spher2)
      end do

      deallocate(r, spher1, spher2, torr, temp)
    end if
    
    call this%time_scheme_sub(cf=1._dbl) ; call this%vypis_ocean_sub()
    
  end subroutine init_state_ocean_sub

  subroutine set_boundary_deformation_sub(this, u_up, t_up)
    class(T_ocean),    intent(inout) :: this
    complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    integer                          :: jmsmax

    jmsmax = min(size(t_up),this%jms)

    this%sol%t_up(1:jmsmax) = t_up(1:jmsmax) / this%D_ud
    this%sol%u_up(1:jmsmax) = u_up(1:jmsmax) / this%D_ud

  end subroutine set_boundary_deformation_sub

  subroutine global_rotation_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: i, m
    real(kind=dbl)                :: coeff
    complex(kind=dbl)             :: angularMomentum

    coeff = ((1/this%r_ud-1)**5) / (1/this%r_ud**5-1)

    do m = 0, 1
      angularMomentum = 5 * this%rad_grid%intV_fn(this%rad_grid%rr * this%sol%velocity_i_fn(1,m,0)) * coeff
        
      do i = 1, this%nd+1
        this%sol%torr(3*(i-1)+1, m+2) = this%sol%torr(3*(i-1)+1, m+2) - angularMomentum * this%rad_grid%rr(i)
      end do
    end do

  end subroutine global_rotation_sub

end submodule OceanExternalRoutines