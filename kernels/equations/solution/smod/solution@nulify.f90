submodule (solution) nulify
  implicit none; contains
  
  module procedure nulify_solution_sub
    
    if ( allocated(this%temp) ) this%temp = czero
    if ( allocated(this%torr) ) this%torr = czero
    if ( allocated(this%mech) ) this%mech = czero
    
  end procedure nulify_solution_sub
  
end submodule nulify