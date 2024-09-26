submodule (OceanMod) Set_bnd_ocean
  implicit none; contains
  
  module subroutine set_boundary_deformation_sub(this, u_up, t_up)
    class(T_ocean),    intent(inout) :: this
    complex(kind=dbl), intent(in)    :: u_up(:), t_up(:)
    integer                          :: jmsmax

    jmsmax = min(size(t_up),this%jms)

    this%bnd%t_up(1:jmsmax) = t_up(1:jmsmax) / this%D_ud
    this%bnd%u_up(1:jmsmax) = u_up(1:jmsmax) / this%D_ud

  end subroutine set_boundary_deformation_sub
  
end submodule Set_bnd_ocean