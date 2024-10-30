submodule (matrices) init
  implicit none; contains
  
  module procedure init_matrices_sub
    
    this%nd        = nd
    this%jmax      = jmax
    this%grid_type = grid_type
    
  end procedure init_matrices_sub
  
  module procedure deallocate_matrices_sub
    integer :: j
    
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
    
  end procedure deallocate_matrices_sub
  
end submodule init