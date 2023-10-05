submodule(RadialGrid) Interpolation
  implicit none
  
  contains
  
  pure function interpolation_fn(this, dimOut, i, rr1, field) result(resField)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, dimOut
    real(kind=dbl),      intent(in) :: rr1(:)
    complex(kind=dbl),   intent(in) :: field(:,:)
    complex(kind=dbl), allocatable  :: resField(:)
    integer                         :: ii, dimjms, dimnd
    
    dimnd  = size(field,1)
    dimjms = min(dimOut, size(field,2))
    
    allocate( resField(dimOut) ) ; resField = czero
    
    if (i == 1) then
      resField( 1:dimjms ) = field( 1, 1:dimjms )
    
    else if (i == this%nd+1) then
      resField( 1:dimjms ) = field( dimnd, 1:dimjms )
    
    else
      do ii = 1, dimnd-1
        if ( (this%rr(i) >= rr1(ii)) .and. (this%rr(i) <= rr1(ii+1)) ) then
          resField(1:dimjms) = ( ( this%rr(i) - rr1(ii)    ) * field(ii+1, 1:dimjms) + &
                             &   ( rr1(ii+1)  - this%rr(i) ) * field(ii  , 1:dimjms)   ) / ( rr1(ii+1) - rr1(ii) )
          exit
        end if
      end do
    end if
    
  end function interpolation_fn
  
end submodule Interpolation
