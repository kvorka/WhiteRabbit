submodule(Matrices) Matrices_init
  implicit none; contains
  
  module pure subroutine init_matrices_sub(this, nd, jmax, grid_type)
    class(T_matrices), intent(inout) :: this
    integer,           intent(in)    :: nd, jmax
    character(len=*),  intent(in)    :: grid_type
    
    this%nd        = nd
    this%jmax      = jmax
    this%grid_type = grid_type
    
  end subroutine init_matrices_sub
  
  module pure subroutine deallocate_matrices_sub(this)
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
  
end submodule Matrices_init