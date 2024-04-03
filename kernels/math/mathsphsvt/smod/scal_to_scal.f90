submodule (Sphsvt) scal_to_scal
  implicit none; contains
  
  module pure subroutine scal2scal_jm_to_mj_sub(this, cjm, cab, ncab, cabpadding)
    class(T_sphsvt),   intent(in)    :: this
    integer,           intent(in)    :: ncab, cabpadding
    complex(kind=dbl), intent(in)    :: cjm(*)
    complex(kind=dbl), intent(inout) :: cab(ncab,*)
    integer                          :: m, j
    
    do m = 0, this%jmax
      do j = m, this%jmax
        cab(cabpadding, m*this%jmax3-m*(m+1)/2+j+1) = cjm(j*(j+1)/2+m+1)
      end do
    end do
    
  end subroutine scal2scal_jm_to_mj_sub
  
  module pure subroutine scal2scal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
    class(T_sphsvt),   intent(in)    :: this
    integer,           intent(in)    :: ncr, ncjm, crpadding, cjmpadding
    complex(kind=dbl), intent(in)    :: cr(ncr,*)
    complex(kind=dbl), intent(inout) :: cjm(ncjm,*)
    integer                          :: j, m, ijm, imj
    
    do j = 0, this%jmax
      m = 0
        ijm = j*(j+1)/2+m+1
        imj = m*this%jmax3-m*(m+1)/2+j+1
          cjm(cjmpadding,ijm)%re = cr(crpadding,imj)%re
          cjm(cjmpadding,ijm)%im = zero
      
      do m = 1, j
        ijm = j*(j+1)/2+m+1
        imj = m*this%jmax3-m*(m+1)/2+j+1
          cjm(cjmpadding,ijm) = cr(crpadding,imj)
      end do
    end do
    
  end subroutine scal2scal_mj_to_jm_sub
  
end submodule scal_to_scal