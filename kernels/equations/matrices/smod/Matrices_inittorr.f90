submodule(Matrices) Matrices_inittorr
  implicit none; contains
  
  module pure subroutine init_mtorr_sub(this)
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
  
end submodule Matrices_inittorr