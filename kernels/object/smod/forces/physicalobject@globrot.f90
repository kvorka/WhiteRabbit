submodule (physicalobject) globrot
  implicit none; contains
  
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
  
end submodule globrot