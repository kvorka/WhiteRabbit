module Matrices
  use Math
  use Matrix
  implicit none; private
  
  type, public :: T_matrices
    class(T_matrix), allocatable :: temp(:), torr(:), mech(:)
    integer                      :: nd, jmax
    character(len=5)             :: grid_type
    
    contains
    
    procedure :: init_sub       => init_matrices_sub
    procedure :: init_mtemp_sub => init_mtemp_sub
    procedure :: init_mtorr_sub => init_mtorr_sub
    procedure :: init_mmech_sub => init_mmech_sub
    procedure :: deallocate_sub => deallocate_matrices_sub
    
  end type T_matrices

  contains
  
  subroutine init_matrices_sub(this, nd, jmax, grid_type)
    class(T_matrices), intent(inout) :: this
    integer,           intent(in)    :: nd, jmax
    character(len=*),  intent(in)    :: grid_type
    
    this%nd        = nd
    this%jmax      = jmax
    this%grid_type = grid_type
    
  end subroutine init_matrices_sub
  
  subroutine init_mtemp_sub(this)
    class(T_matrices), intent(inout) :: this
    integer                          :: j
    
    allocate( this%temp(0:this%jmax) )
    
    select case (this%grid_type)
      case('chebv')
        do j = 0, this%jmax
          call this%temp(j)%init_sub(3*this%nd+1, 5, 5)
        end do
      
      case('homog')
        do j = 0, this%jmax
          call this%temp(j)%init_sub(3*this%nd+1, 3, 3)
        end do
    end select
    
  end subroutine init_mtemp_sub
  
  subroutine init_mtorr_sub(this)
    class(T_matrices), intent(inout) :: this
    integer                          :: j
    
    allocate( this%torr(1:this%jmax) )
    
    select case (this%grid_type)
      case('chebv')
        do j = 1, this%jmax
          call this%torr(j)%init_sub(3*this%nd+1, 5, 5)
        end do
        
      case('homog')
        do j = 1, this%jmax
          call this%torr(j)%init_sub(3*this%nd+1, 3, 3)
        end do
    end select
    
  end subroutine init_mtorr_sub
  
  subroutine init_mmech_sub(this)
    class(T_matrices), intent(inout) :: this
    integer                          :: j
    
    allocate( this%mech(1:this%jmax) )
    
    select case (this%grid_type)
      case('chebv')
        do j = 1, this%jmax
          call this%mech(j)%init_sub(6*this%nd+2, 11, 11)
        end do
      
      case('homog')
        do j = 1, this%jmax
          call this%mech(j)%init_sub(6*this%nd+2, 7, 7)
        end do
    end select
      
  end subroutine init_mmech_sub
  
  subroutine deallocate_matrices_sub(this)
    class(T_matrices), intent(inout) :: this
    integer                          :: j
    
    if ( allocated(this%temp) ) then
      do j = 0, this%jmax
        call this%temp(j)%deallocate_sub()
      end do
      
      deallocate( this%temp )
    end if
    
    if ( allocated(this%torr) ) then
      do j = 1, this%jmax
        call this%torr(j)%deallocate_sub()
      end do
      
      deallocate( this%torr )
    end if
    
    if ( allocated(this%mech) ) then
      do j = 1, this%jmax
        call this%mech(j)%deallocate_sub()
      end do
      
      deallocate( this%mech )
    end if
    
  end subroutine deallocate_matrices_sub
  
end module Matrices
