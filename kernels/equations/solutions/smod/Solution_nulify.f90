submodule (Solution) Solution_nulify
  implicit none; contains
  
  module pure subroutine nulify_solution_sub(this)
    class(T_solution), intent(inout) :: this
    
    if ( allocated(this%temp) ) this%temp = czero
    if ( allocated(this%torr) ) this%torr = czero
    if ( allocated(this%mech) ) this%mech = czero
    
  end subroutine nulify_solution_sub
  
end submodule Solution_nulify