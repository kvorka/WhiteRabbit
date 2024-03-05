module Legendre_polynomials
  use Math
  implicit none; public; contains
  
  pure subroutine pmm_setup_8_sub(pmm)
    real(kind=dbl), intent(out) :: pmm(*)
    integer                     :: i
    
    do concurrent ( i = 1:8 )
      pmm(i) = 1._dbl
    end do
    
  end subroutine pmm_setup_8_sub
  
  pure subroutine pmm_recursion_8_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
    
    do concurrent ( i = 1:8 )
      pmm(i) = fac * sinx(i) * pmm(i)
    end do
    
  end subroutine pmm_recursion_8_sub
  
  pure subroutine pmm_setup_16_sub(pmm)
    real(kind=dbl), intent(out) :: pmm(*)
    integer                     :: i
    
    do concurrent ( i = 1:16 )
      pmm(i) = 1._dbl
    end do
    
  end subroutine pmm_setup_16_sub
  
  pure subroutine pmm_recursion_16_sub(m, sinx, pmm)
    integer,        intent(in)    :: m
    real(kind=dbl), intent(in)    :: sinx(*)
    real(kind=dbl), intent(inout) :: pmm(*)
    integer                       :: i
    real(kind=dbl)                :: fac
    
    fac = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) )
    
    do concurrent ( i = 1:16 )
      pmm(i) = fac * sinx(i) * pmm(i)
    end do
    
  end subroutine pmm_recursion_16_sub
  
end module Legendre_polynomials