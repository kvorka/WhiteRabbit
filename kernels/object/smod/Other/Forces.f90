submodule(PhysicalObject) Forces
  implicit none ; contains
  
  module pure subroutine coriolis_rr_jml_sub(this, v, coriolis)
    class(T_physicalObject), intent(in)    :: this
    complex(kind=dbl),       intent(in)    :: v(:)
    complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    
    select case (this%scaling)
      case ('christ')
        call ezvv_sub(this%jmax, 2._dbl, v, coriolis)
      
      case('basics')
        call ezvv_sub(this%jmax, 2/this%Ek, v, coriolis)
    end select
    
  end subroutine coriolis_rr_jml_sub
  
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
  
  module pure subroutine global_rotation_sub(this)
    class(T_physicalObject), intent(inout) :: this
    integer                                :: ir, is, ijm
    real(kind=dbl)                         :: coeff
    complex(kind=dbl)                      :: angularMomentum
    complex(kind=dbl),       allocatable   :: angularMomentum_rr(:)
    
    coeff = 5 * ((1/this%r_ud-1)**5) / (1/this%r_ud**5-1)
    
    allocate( angularMomentum_rr(this%nd+1) )
    
      do ijm = 2, 3
        do concurrent ( ir = 1:this%nd+1 )
          angularMomentum_rr(ir) = this%rad_grid%rr(ir) * this%v_rr_fn(ir,0,ijm)
        end do
        
        angularMomentum = coeff * this%rad_grid%intV_fn(angularMomentum_rr)
        
        do concurrent ( ir = 1:this%nd+1 )
          is = 3*(ir-1)+1
          this%sol%torr(is,ijm) = this%sol%torr(is,ijm) - angularMomentum * this%rad_grid%rr(ir)
        end do
      end do
      
    deallocate( angularMomentum_rr )
    
  end subroutine global_rotation_sub
  
end submodule Forces