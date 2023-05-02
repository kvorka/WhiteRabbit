submodule(RadialGrid) Interpolation
  implicit none

  contains

  pure function interpolation_fn(this, dimOut, i, rr1, field) result(resField)
    class(T_radialGrid), intent(in) :: this
    integer,             intent(in) :: i, dimOut
    real(kind=dbl),      intent(in) :: rr1(:)
    complex(kind=dbl),   intent(in) :: field(:,:)
    complex(kind=dbl)               :: resField(dimOut)
    integer                         :: ii, dimensions(2)
    
    dimensions = shape(field) ; resField = czero
    
    if (i == 1) then
      if (dimensions(2) >= dimOut) then
        resField(1:dimOut) = field(1,1:dimOut)
      else
        resField(1:dimensions(2)) = field(1,:)
      end if
    
    else if (i == this%nd+1) then
      if (dimensions(2) >= dimOut) then
        resField(1:dimOut) = field(dimensions(1), 1:dimOut)
      else
        resField(1:dimensions(2)) = field(dimensions(1),:)
      end if
    
    else
      do ii = 1, (dimensions(1)-1)
        if ( (this%rr(i) >= rr1(ii)) .and. (this%rr(i) <= rr1(ii+1))) then
          if (dimensions(2) >= dimOut) then
            resField(1:dimOut) = ( ( this%rr(i) - rr1(ii)    ) * field(ii+1, 1:dimOut) + &
                               &   (  rr1(ii+1) - this%rr(i) ) * field(ii  , 1:dimOut)   ) / ( rr1(ii+1) - rr1(ii) )
          else
            resField(1:dimensions(2)) = ( ( this%rr(i) - rr1(ii)    ) * field(ii+1,:) + &
                                      &   ( rr1(ii+1)  - this%rr(i) ) * field(ii  ,:)   ) / ( rr1(ii+1) - rr1(ii) )
          end if
          
          exit
          
        end if
      end do
    
    end if
    
  end function interpolation_fn

end submodule Interpolation
