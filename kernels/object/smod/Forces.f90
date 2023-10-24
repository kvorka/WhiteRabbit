submodule(PhysicalObject) Forces
  implicit none
  
  contains

  pure function buoy_rr_fn(this, ir, ijm) result(buoy)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    complex(kind=dbl)                   :: buoy
    
    buoy = this%Ra * this%alpha_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) ) * this%sol%temp_fn(ir,ijm)
    
  end function buoy_rr_fn
  
  pure function buoy_rr_jml_fn(this, ir) result(gdrho)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       allocatable :: gdrho(:)
    
    allocate( gdrho(this%jmv) ) ; gdrho = czero
    
    gdrho = ersv_fn(this%jmax, this%Ra * this%alpha_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) ) * this%sol%temp_jm_fn(ir))
    
  end function buoy_rr_jml_fn

  pure function coriolis_rr_jml_fn(this, ir, v) result(coriolis)
    class(T_physicalObject),     intent(in) :: this
    integer,           optional, intent(in) :: ir
    complex(kind=dbl), optional, intent(in) :: v(:)
    real(kind=dbl)                          :: fac
    complex(kind=dbl), allocatable          :: coriolis(:)
    
    allocate( coriolis(this%jmv) ) ; coriolis = czero
    
    if ( present(ir) ) coriolis = ezvv_fn( this%jmax, this%sol%velocity_jml_fn(ir) )
    if ( present(v)  ) coriolis = ezvv_fn( this%jmax, v )
    
    fac = 2 / this%Ek ; coriolis = coriolis * fac
    
  end function coriolis_rr_jml_fn

  subroutine coriolis_rr_jml_sub(this, v, coriolis) 
    class(T_physicalObject), intent(in)    :: this
    complex(kind=dbl),       intent(in)    :: v(:)
    complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    real(kind=dbl)                         :: fac
    
    call ezvv_sub(this%jmax, 2/this%Ek, v, coriolis)

  end subroutine coriolis_rr_jml_sub

  subroutine buoy_rr_jml_sub(this, ir, T, force)
    class(T_physicalObject), intent(in)    :: this
    integer,                 intent(in)    :: ir
    complex(kind=dbl),       intent(in)    :: T(:)
    complex(kind=dbl),       intent(inout) :: force(:,:)
    integer                                :: ijm, ij
    real(kind=dbl)                         :: fac, fac1, fac2
      
    fac = this%Ra * this%alpha_fn(ir) * this%gravity%g_fn(this%rad_grid%rr(ir))
    
    do ij = 1, this%jmax
      fac1 = -sqrt( (ij  ) / (2*ij+1._dbl) ) * fac
      fac2 = +sqrt( (ij+1) / (2*ij+1._dbl) ) * fac
      
      do concurrent ( ijm = ij*(ij+1)/2+1:ij*(ij+1)/2+ij+1 )
        force(1,ijm) = force(1,ijm) + fac1 * T(ijm)
        force(3,ijm) = force(3,ijm) + fac2 * T(ijm)
      end do
    end do
      
  end subroutine buoy_rr_jml_sub

  subroutine global_rotation_sub(this)
    class(T_physicalObject), intent(inout) :: this
    integer                                :: ir, is, ijm
    real(kind=dbl)                         :: coeff
    complex(kind=dbl)                      :: angularMomentum
    
    coeff = 5 * ((1/this%r_ud-1)**5) / (1/this%r_ud**5-1)
    
    do ijm = 2, 3
      angularMomentum = coeff * this%rad_grid%intV_fn(this%rad_grid%rr * this%sol%velocity_i_fn(0,ijm))
        
      do concurrent ( ir = 1:this%nd+1 )
        is = 3*(ir-1)+1 ; this%sol%torr(is,ijm) = this%sol%torr(is,ijm) - angularMomentum * this%rad_grid%rr(ir)
      end do
    end do
    
  end subroutine global_rotation_sub
  
end submodule Forces
