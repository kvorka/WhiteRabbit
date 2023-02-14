module LU
  use Math
  implicit none

  public :: ludecomposition_sub
  public :: lusolution_sub
  public :: multipl1_fn
  public :: multipl2_fn

  contains

  subroutine ludecomposition_sub(n, ld, lu, a, al, indx)
    integer,        intent(in)    :: n, ld, lu
    real(kind=dbl), intent(inout) :: a(:,:), al(:,:)
    integer,        intent(out)   :: indx(:)
    integer                       :: i, k, l
    real(kind=dbl)                :: dum
    real(kind=dbl), allocatable   :: pom(:)

    l = ld
      do i = 1, ld
        do k = (ld+2-i), ld+1+lu
          a(k-l,i) = a(k,i)
        end do
    
        l = l-1; a(ld+1+lu-l:ld+1+lu,i) = 0._dbl
      end do
    
    do k = 1, n
      l = min(ld+k,n)
        i = maxloc(abs(a(1,k:l)),1)+k-1; dum = a(1,i); indx(k) = i
      
        if (i /= k) then
          allocate(pom(ld+1+lu))
            pom    = a(:,k)
            a(:,k) = a(:,i)
            a(:,i) = pom
          deallocate(pom)
        end if
      
      do i = k+1, l
        dum = a(1,i)/a(1,k); al(i-k,k) = dum
          a(1:ld+lu,i) = a(2:ld+1+lu,i) - dum*a(2:ld+1+lu,k)
          a(ld+1+lu,i) = 0._dbl
      end do
    end do

  end subroutine ludecomposition_sub

  subroutine lusolution_sub(n, ld, lu, a, al, indx, b)
    integer,           intent(in)    :: n, ld, lu
    real(kind=dbl),    intent(in)    :: a(:,:), al(:,:)
    integer,           intent(in)    :: indx(:)
    complex(kind=dbl), intent(inout) :: b(:)
    integer                          :: i, k, l
    complex(kind=dbl)                :: dum

    l = ld
      do k = 1, n
        i = indx(k)
          if (i /= k) then
            dum  = b(k)
            b(k) = b(i)
            b(i) = dum
          end if

        if (l < n) l = l + 1
          b(k+1:l) = b(k+1:l) - al(1:l-k,k)*b(k)
      end do

    l = 1
      do i = n, 1, -1
        b(i) = (b(i) - sum(a(2:l,i)*b(i+1:i+l-1)))/a(1,i)
        if (l < (ld+1+lu)) l = l + 1
      end do

  end subroutine lusolution_sub
  
  pure complex(kind=dbl) function multipl1_fn(i, n, ld, lu, a, x)
    integer,           intent(in) :: i, n, ld, lu
    real(kind=dbl),    intent(in) :: a(:,:)
    complex(kind=dbl), intent(in) :: x(:)

    multipl1_fn = sum( a(max(1,ld+2-i):min(ld+1+lu,ld+1+n-i),i) * x(max(1,i-ld):min(i+lu,n)) )

  end function multipl1_fn
  
  pure function multipl2_fn(i, n, ld, lu, a, x) result(multipl)
    integer,           intent(in) :: i, n, ld, lu
    real(kind=dbl),    intent(in) :: a(:,:)
    complex(kind=dbl), intent(in) :: x(:,:)
    integer                       :: j, indx1, indx2, indm1, indm2
    real(kind=dbl),   allocatable :: mat(:) 
    complex(kind=dbl)             :: multipl(size(x,2))

    indx1 = max(1,i-ld) ; indm1 = max(1,ld+2-i)
    indx2 = min(i+lu,n) ; indm2 = min(ld+1+lu,ld+1+n-i)
    
    allocate( mat(indm2-indm1+1) ) ; mat = a(indm1:indm2,i)
      do j = 1, size(x,2)
        multipl(j) = sum( mat(:) * x(indx1:indx2,j) )
      end do
    deallocate( mat )

  end function multipl2_fn

end module LU