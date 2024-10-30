submodule (matrices) inittorr
  implicit none; contains
  
  module procedure init_mtorr_sub
    integer :: j
    
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
    
  end procedure init_mtorr_sub
  
end submodule inittorr