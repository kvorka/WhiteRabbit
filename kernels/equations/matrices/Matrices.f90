module Matrices
  use matrix
  implicit none
  
  type, public :: T_matrices
    class(T_matrix), allocatable :: temp(:), torr(:), mech(:)
    integer                      :: nd, jmax
    character(len=5)             :: grid_type
    
    contains
    
    procedure :: init_sub       => init_matrices_sub
    procedure :: deallocate_sub => deallocate_matrices_sub
    procedure :: init_mtemp_sub, init_mtorr_sub, init_mmech_sub
    
  end type T_matrices
  
  interface
    module pure subroutine init_matrices_sub(this, nd, jmax, grid_type)
      class(T_matrices), intent(inout) :: this
      integer,           intent(in)    :: nd, jmax
      character(len=*),  intent(in)    :: grid_type
    end subroutine init_matrices_sub
    
    module pure subroutine deallocate_matrices_sub(this)
      class(T_matrices), intent(inout) :: this
    end subroutine deallocate_matrices_sub

    module pure subroutine init_mtemp_sub(this)
      class(T_matrices), intent(inout) :: this
    end subroutine init_mtemp_sub
    
    module pure subroutine init_mtorr_sub(this)
      class(T_matrices), intent(inout) :: this
    end subroutine init_mtorr_sub
    
    module pure subroutine init_mmech_sub(this)
      class(T_matrices), intent(inout) :: this
    end subroutine init_mmech_sub
  end interface
  
end module Matrices