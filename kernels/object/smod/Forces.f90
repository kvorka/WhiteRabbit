submodule(PhysicalObject) Forces
  implicit none
  
  contains

  pure function buoy_rr_fn(this, i, ijm) result(buoy)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i, ijm
    complex(kind=dbl)                   :: buoy
    
    buoy = this%Ra * this%alpha_fn(i) * this%gravity%g_fn( this%rad_grid%rr(i) ) * this%sol%temp2_fn(i,ijm)
    
  end function buoy_rr_fn
  
  pure function buoy_rr_jml_fn(this, i) result(gdrho)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: gdrho(this%jmv)
    
    gdrho = ersv_fn(this%jmax, this%Ra * this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i)) * this%sol%temp_jm_fn(i))
    
  end function buoy_rr_jml_fn

  pure function coriolis_rr_jml_fn(this, i, v) result(coriolis)
    class(T_physicalObject),     intent(in) :: this
    integer,           optional, intent(in) :: i
    complex(kind=dbl), optional, intent(in) :: v(:)
    real(kind=dbl)                          :: fac
    complex(kind=dbl)                       :: coriolis(this%jmv)
  
    if ( present(i) ) coriolis = ezvv_fn( this%jmax, this%sol%velocity_jml_fn(i) )
    if ( present(v) ) coriolis = ezvv_fn( this%jmax, v )

    fac = 2 / this%Ek ; coriolis = coriolis * fac

  end function coriolis_rr_jml_fn

  subroutine coriolis_rr_jml_sub(this, v, coriolis) 
    class(T_physicalObject), intent(in)    :: this
    complex(kind=dbl),       intent(in)    :: v(:)
    complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    real(kind=dbl)                         :: fac
    
    call ezvv_sub(this%jmax, 2/this%Ek, v, coriolis)

  end subroutine coriolis_rr_jml_sub

  subroutine buoy_rr_jml_sub(this, i, T, force)
    class(T_physicalObject), intent(in)    :: this
    integer,                 intent(in)    :: i
    complex(kind=dbl),       intent(in)    :: T(:)
    complex(kind=dbl),       intent(inout) :: force(:,:)
    integer                                :: ijm, j
    real(kind=dbl)                         :: fac, fac1, fac2
      
    fac = this%Ra * this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i))
    
    do j = 1, this%jmax
      fac1 = -sqrt( (j  ) / (2*j+1._dbl) ) * fac
      fac2 = +sqrt( (j+1) / (2*j+1._dbl) ) * fac
      
      do concurrent ( ijm = j*(j+1)/2+1:j*(j+1)/2+j+1 )
        force(1,ijm) = force(1,ijm) + fac1 * T(ijm)
        force(3,ijm) = force(3,ijm) + fac2 * T(ijm)
      end do
    end do
      
  end subroutine buoy_rr_jml_sub

  subroutine global_rotation_sub(this)
    class(T_physicalObject), intent(inout) :: this
    integer                                :: i, m
    real(kind=dbl)                         :: coeff
    complex(kind=dbl)                      :: angularMomentum

    coeff = ((1/this%r_ud-1)**5) / (1/this%r_ud**5-1)

    do m = 0, 1
      angularMomentum = 5 * this%rad_grid%intV_fn(this%rad_grid%rr * this%sol%velocity_i_fn(1,m,0)) * coeff
        
      do i = 1, this%nd+1
        this%sol%torr(3*(i-1)+1, m+2) = this%sol%torr(3*(i-1)+1, m+2) - angularMomentum * this%rad_grid%rr(i)
      end do
    end do

  end subroutine global_rotation_sub
  
end submodule Forces
