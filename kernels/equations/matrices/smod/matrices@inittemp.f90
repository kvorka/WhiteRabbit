submodule (matrices) inittemp
  implicit none; contains
  
  module procedure init_mtemp_sub
    integer :: j
    
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
    
  end procedure init_mtemp_sub
  
end submodule inittemp