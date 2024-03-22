submodule (SphericalHarmonics) vctol
  implicit none ; contains
  
  pure subroutine grid_op_vctol_sub(nfour, nstep, grid)
    integer,                intent(in)    :: nfour, nstep
    real(kind=dbl), target, intent(inout) :: grid(*)
  end subroutine grid_op_vctol_sub
  
  module subroutine vctol_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    complex(kind=dbl),    allocatable   :: cc(:), cr(:)
    
    allocate( cc(this%jms2), cr(this%jms2) )
    
    call random_number( cc(1:this%jmax2+1)%re )
    call random_number( cc(this%jmax2+2:)%re )
    call random_number( cc(this%jmax2+2:)%im )
    
    do
      !Set the test array to zeros
      call zero_carray_sub( this%jms2, cr(1) )
      
      !Transform
      call this%lege_transform_sub( 1, 1, cc(1), cr(1), grid_op_vctol_sub )
      
      !Presision check
      if ( maxval( abs( abs(cc/cr) - 1 ) ) <= 1.0d-4 ) then
        exit
      else if ( this%tolm < 1.0d-100) then
        this%tolm = 1.0d-90
        exit
      else
        write(*,*) maxval( abs( abs(cc/cr) - 1 ) )
        this%tolm = this%tolm / 10
      end if
    end do
    
    deallocate( cc, cr )
    
  end subroutine vctol_sub
  
end submodule vctol