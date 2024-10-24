module grid_ops
  use math
  implicit none
  
  interface
    module pure subroutine grid_op_vcvv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(*)
    end subroutine grid_op_vcvv_sub
    
    module pure subroutine grid_op_vcvxv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(*)
    end subroutine grid_op_vcvxv_sub
  end interface
  
end module grid_ops