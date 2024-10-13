submodule (lege_poly) init
  implicit none; contains
  
  module pure subroutine init_lege_sub(this, jmax)
    class(T_legep), intent(inout) :: this
    integer,        intent(in)    :: jmax
    integer                       :: jms, ij, im
    
    jms = (jmax+2)*(jmax+1)/2

    this%jmax = jmax
    
    allocate( this%amj(jms), this%bmj(jms), this%cmm(0:jmax) )
    
    do im = 0, jmax
      if ( im == 0 ) then
        this%cmm(im) = one / s4pi
      else
        this%cmm(im) = -sqrt( (2*im+one) / (2*im) )
      end if
      
      do ij = im+1, jmax
        this%amj(im*(jmax+1)-im*(im+1)/2+ij+1) = sqrt((2*ij-1)*(2*ij+one)                    /(         (ij-im)*(ij+im)))
        this%bmj(im*(jmax+1)-im*(im+1)/2+ij+1) = sqrt(         (2*ij+one)*(ij-im-1)*(ij+im-1)/((2*ij-3)*(ij-im)*(ij+im)))
      end do
    end do
    
  end subroutine init_lege_sub
  
  module pure subroutine deallocate_lege_sub(this)
    class(T_legep), intent(inout) :: this
    
    if ( allocated(this%amj) ) deallocate( this%amj )
    if ( allocated(this%bmj) ) deallocate( this%bmj )
    if ( allocated(this%cmm) ) deallocate( this%cmm )
    
  end subroutine deallocate_lege_sub
  
end submodule init