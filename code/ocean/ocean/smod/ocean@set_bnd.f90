submodule (ocean) set_bnd
  implicit none; contains
  
  module procedure set_boundary_deformation_sub
    integer :: jmsmax

    jmsmax = min(size(t_up),this%jms)

    this%bnd%t_up(1:jmsmax) = t_up(1:jmsmax) / this%D_ud
    this%bnd%u_up(1:jmsmax) = u_up(1:jmsmax) / this%D_ud

  end procedure set_boundary_deformation_sub
  
end submodule set_bnd