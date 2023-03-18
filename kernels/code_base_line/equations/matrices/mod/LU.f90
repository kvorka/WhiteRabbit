module LU
  use Math
  implicit none; private

  public :: ludecomposition_sub
  public :: lusolution_sub
  public :: multipl1_fn
  public :: multipl2_fn

  contains

  subroutine ludecomposition_sub(n, ld, lu, Upper, Lower, Indx)
    integer,        intent(in)    :: n, ld, lu
    real(kind=dbl), intent(inout) :: Upper(:,:)
    real(kind=dbl), intent(out)   :: Lower(:,:)
    integer,        intent(out)   :: Indx(:)
    integer                       :: i, j, k, ldu
    real(kind=dbl), allocatable   :: pom(:)

    k = ld; ldu = ld+1+lu
    
    do i = 1, ld
      Upper(ld+2-i-k:ldu-k,i) = Upper(ld+2-i:ldu,i)

      k = k-1
        Upper(ldu-k:ldu,i) = 0._dbl
    end do
    
    allocate(pom(1:ldu))
      do j = 1, n
        k = min(ld+j,n); i = maxloc(abs(Upper(1,j:k)),1)+j-1; Indx(j) = i
      
        if (i /= j) then
          pom        = Upper(:,j)
          Upper(:,j) = Upper(:,i)
          Upper(:,i) = pom
        end if
      
        do i = j+1, k
          Lower(i-j,j)     = Upper(1,i) / Upper(1,j)
          Upper(1:ldu-1,i) = Upper(2:ldu,i) - Lower(i-j,j) * Upper(2:ldu,j)
          Upper(ldu,i) = 0._dbl
        end do
      end do
    deallocate(pom)

  end subroutine ludecomposition_sub

  subroutine lusolution_sub(n, ld, lu, Upper, Lower, Indx, b)
    integer,           intent(in)    :: n, ld, lu
    real(kind=dbl),    intent(in)    :: Upper(:,:), Lower(:,:)
    integer,           intent(in)    :: Indx(:)
    complex(kind=dbl), intent(inout) :: b(:)
    integer                          :: i, j, k, ldu
    complex(kind=dbl)                :: dum

    ldu = ld + 1 + lu
    
    do j = 1, n
      i = Indx(j)
        if (i /= j) then
          dum  = b(i)
          b(i) = b(j)
          b(j) = dum
        end if
        
      k = min(n,ld+j) ; b(j+1:k) = b(j+1:k) - Lower(1:k-j,j) * b(j)
    end do
    
    do i = n, 1, -1
      k = min(ldu,n-i+1) ; b(i) = ( b(i) - sum( Upper(2:k,i)*b(i+1:i+k-1) ) ) / Upper(1,i)
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