submodule (sphsvt) init
  implicit none; contains
  
  module pure subroutine init_sphsvt_sub(this, jmax)
    class(T_sphsvt), intent(inout) :: this
    integer,         intent(in)    :: jmax
    
    this%jmax  = jmax
    this%jmax1 = jmax+1
    this%jmax2 = jmax+2
    this%jmax3 = jmax+3

    this%jms  =     ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jms1 =     ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    this%jms2 =     ((jmax+2)*(jmax+3)/2 + (jmax+2)) + 1

    this%jmv  = 3 * ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jmv1 = 3 * ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    
  end subroutine init_sphsvt_sub
  
end submodule init