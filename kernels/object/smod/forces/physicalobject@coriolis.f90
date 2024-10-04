submodule (physicalobject) coriolis
  implicit none; contains
  
  module pure subroutine coriolis_rr_jml_sub(this, v, coriolis)
    class(T_physicalObject), intent(in)    :: this
    complex(kind=dbl),       intent(in)    :: v(:)
    complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    
    select case (this%scaling)
      case ('christ')
        call ezvv_sub(this%jmax, 2._dbl, v, coriolis)
      
      case('basics')
        call ezvv_sub(this%jmax, 2/this%Ek, v, coriolis)
    end select
    
  end subroutine coriolis_rr_jml_sub
  
end submodule coriolis