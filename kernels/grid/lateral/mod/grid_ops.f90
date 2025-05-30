module grid_ops
  use math
  implicit none
  
  interface
    module  subroutine grid_op_save_sub(nfour, tempgrid, grid, padding)
      integer,        intent(in)  :: nfour, padding
      real(kind=dbl), intent(in)  :: tempgrid(step,nfour)
      real(kind=dbl), intent(out) :: grid(*)
    end subroutine grid_op_save_sub
    
    module  subroutine grid_op_load_sub(nfour, tempgrid, grid, padding)
      integer,        intent(in)  :: nfour, padding
      real(kind=dbl), intent(out) :: tempgrid(step,nfour)
      real(kind=dbl), intent(in)  :: grid(*)
    end subroutine grid_op_load_sub
    
    module subroutine grid_op_vcvv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcvv_sub
    
    module subroutine grid_op_vcvxv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcvxv_sub
    
    module subroutine grid_op_vcvv_vcvgv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcvv_vcvgv_sub
    
    module subroutine grid_op_vcvgv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcvgv_sub
    
    module subroutine grid_op_vcsv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcsv_sub
    
    module subroutine grid_op_vcst_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcst_sub
    
    module subroutine grid_op_vcss_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcss_sub
    
    module subroutine grid_op_vcss_add_vcvv_sub(nfour, grid, tempgrid)
      integer,                intent(in)    :: nfour
      real(kind=dbl), target, intent(inout) :: grid(step,*)
      real(kind=dbl), target, intent(out)   :: tempgrid(step,*)
    end subroutine grid_op_vcss_add_vcvv_sub
  end interface
  
end module grid_ops