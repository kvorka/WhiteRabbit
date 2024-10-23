submodule (lege_poly) step16p
  implicit none; contains
  
  module pure subroutine mmset_16_sub(this, m, sinx, pmm, pmj2, pmj1, pmj0)
    class(T_legep), intent(in)    :: this
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(16)
    real(kind=dbl), intent(inout) :: pmm(16)
    real(kind=dbl), intent(out)   :: pmj2(16), pmj1(16), pmj0(16)
    integer                       :: i2
    
    if ( m /= 0 ) then
      do concurrent ( i2 = 1:16 )
        pmm(i2) = this%cmm(m) * sinx(i2) * pmm(i2)
      end do
    else
      do concurrent ( i2 = 1:16 )
        pmm(i2) = this%cmm(0)
      end do
    end if
    
    do concurrent ( i2 = 1:16 )
      pmj2(i2) = zero
      pmj1(i2) = zero
      pmj0(i2) = pmm(i2)
    end do
    
  end subroutine mmset_16_sub
  
  module pure subroutine recursion_16_sub(this, mj, cosx, pmj2, pmj1, pmj0)
    class(T_legep), intent(in)    :: this
    integer,        intent(in)    :: mj
    real(kind=dbl), intent(in)    :: cosx(16)
    real(kind=dbl), intent(inout) :: pmj2(16), pmj1(16), pmj0(16)
    integer                       :: i2
    
    do concurrent ( i2=1:16 )
      pmj2(i2) = this%bmj(mj) * pmj1(i2)
      pmj1(i2) = this%amj(mj) * pmj0(i2)
      pmj0(i2) = cosx(i2) * pmj1(i2) - pmj2(i2)
    end do
    
  end subroutine recursion_16_sub
  
end submodule step16p