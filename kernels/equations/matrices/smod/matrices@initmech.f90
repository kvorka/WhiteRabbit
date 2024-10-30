submodule (Matrices) initmech
  implicit none; contains
  
  module procedure init_mmech_sub
    integer :: j
    
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
      
  end procedure init_mmech_sub
  
end submodule initmech