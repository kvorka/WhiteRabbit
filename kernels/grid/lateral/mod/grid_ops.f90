module grid_ops
  use math
  implicit none
  
  interface
    module pure subroutine grid_op_vcvv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcvv_sub
    
    module pure subroutine grid_op_vcvxv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcvxv_sub
    
    module pure subroutine grid_op_vcvv_vcvgv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcvv_vcvgv_sub
    
    module pure subroutine grid_op_vcvgv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcvgv_sub
    
    module pure subroutine grid_op_vcsv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcsv_sub
    
    module pure subroutine grid_op_vcst_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcst_sub
    
    module pure subroutine grid_op_vcss_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcss_sub
    
    module pure subroutine grid_op_vcss_add_vcvv_sub(step, nfour, grid)
      integer,                intent(in)    :: nfour, step
      real(kind=dbl), target, intent(inout) :: grid(step,*)
    end subroutine grid_op_vcss_add_vcvv_sub
  end interface
  
end module grid_ops