module grid_ops
  use math
  implicit none
  
  interface
    module subroutine grid_op_vcvv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcvv_sub
    
    module subroutine grid_op_vcvxv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcvxv_sub
    
    module subroutine grid_op_vcvv_vcvgv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcvv_vcvgv_sub
    
    module subroutine grid_op_vcvgv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcvgv_sub
    
    module subroutine grid_op_vcsv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcsv_sub
    
    module subroutine grid_op_vcst_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcst_sub
    
    module subroutine grid_op_vcss_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcss_sub
    
    module subroutine grid_op_vcss_add_vcvv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(16,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(16,*)
    end subroutine grid_op_vcss_add_vcvv_sub
  end interface
  
end module grid_ops