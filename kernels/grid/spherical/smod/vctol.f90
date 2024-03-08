submodule (SphericalHarmonics) vctol
  implicit none ; contains
  
  pure subroutine grid_op_vctol_sub(nfour, nstep, grid)
    integer,                intent(in)    :: nfour, nstep
    real(kind=dbl), target, intent(inout) :: grid(*)
  end subroutine grid_op_vctol_sub
  
  module subroutine vctol_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    integer                             :: j, m, mj
    real(kind=dbl)                      :: diff, maxdiff, xrand
    complex(kind=dbl),    allocatable   :: cc(:), cr(:)
    
    allocate( cc(this%jms2), cr(this%jms2) )
      
      do m = 0, this%jmax2
        do j = m, this%jmax2
          mj = m*this%jmax3-m*(m+1)/2+j+1
          
          if ( m /= 0 ) then
            call random_number( xrand ) ; cc(mj)%re = xrand
            call random_number( xrand ) ; cc(mj)%im = xrand
          else
            call random_number( xrand ) ; cc(mj)%re = xrand
            cc(mj)%im = czero
          end if
        end do
      end do
      
    do
      cr = czero
      
      !Transform
      call this%lege_transform_sub( 1, 1, cc(1), cr(1), grid_op_vctol_sub )
        
      maxdiff = 0._dbl
      
      do mj = 1, this%jms2
        diff = abs( abs( cc(mj) / cr(mj) ) / this%scale / sqrt(4*pi) - 1 )
          
        if ( diff > maxdiff ) maxdiff = diff
      end do
      
      if ( this%tolm < 1.0d-100) then
        this%tolm = 1.0d-90
        exit
      else if ( maxdiff <= 1.0d-4 ) then
        exit
      else
        this%tolm = this%tolm / 10
      end if
    end do
    
  end subroutine vctol_sub
  
end submodule vctol