submodule (PhysicalObject) Buoyancy
  implicit none; contains
  
  module pure complex(kind=dbl) function buoy_rr_fn(this, ir, il, ijm, sgn)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, il, ijm, sgn
    real(kind=dbl)                      :: j
    complex(kind=dbl)                   :: er_buoy_jm
    
    j = i2r_fn( this%j_indx(ijm) )
    er_buoy_jm = this%Ra * this%alpha_rr_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) ) * this%temp_rr_fn(ir,ijm)
    
    select case (il)
      case( -1)
        buoy_rr_fn = +sgn * sqrt((j  )/(2*j+1)) * er_buoy_jm
        
      case ( 0)
        buoy_rr_fn = czero
        
      case (+1)
        buoy_rr_fn = -sgn * sqrt((j+1)/(2*j+1)) * er_buoy_jm
        
    end select
    
  end function buoy_rr_fn
  
  module pure subroutine er_buoy_rr_jm_sub(this, ir, force)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: force(*)
    integer                              :: ijm
    real(kind=dbl)                       :: fac
      
    fac = this%Ra * this%alpha_rr_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) )
    
    do concurrent ( ijm = 1:this%jms )
      force(ijm) = fac * this%temp_rr_fn(ir,ijm)
    end do
    
  end subroutine er_buoy_rr_jm_sub
  
  module pure subroutine buoy_rr_jml_sub(this, ir, T, force)
    class(T_physicalObject), intent(in)    :: this
    integer,                 intent(in)    :: ir
    complex(kind=dbl),       intent(in)    :: T(:)
    complex(kind=dbl),       intent(inout) :: force(:,:)
    integer                                :: ijm, ij
    real(kind=dbl)                         :: fac, fac1, fac2
      
    fac = this%Ra * this%alpha_rr_fn(ir) * this%gravity%g_fn( this%rad_grid%rr(ir) )
    
    do ij = 1, this%jmax
      fac1 = -sqrt( (ij  ) / (2*ij+one) ) * fac
      fac2 = +sqrt( (ij+1) / (2*ij+one) ) * fac
      
      do concurrent ( ijm = jm(ij,0):jm(ij,ij) )
        force(1,ijm) = force(1,ijm) + fac1 * T(ijm)
        force(3,ijm) = force(3,ijm) + fac2 * T(ijm)
      end do
    end do
      
  end subroutine buoy_rr_jml_sub
  
end submodule Buoyancy