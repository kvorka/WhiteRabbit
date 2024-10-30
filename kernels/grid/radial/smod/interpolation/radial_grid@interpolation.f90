submodule (radial_grid) interpolation
  implicit none ; contains
  
  module procedure interpolation_fn
    integer :: ii, dimjms, dimnd
    
    dimnd  = size(field,1) ; dimjms = min(dimOut, size(field,2))
    
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
    
  end procedure interpolation_fn
  
end submodule interpolation
