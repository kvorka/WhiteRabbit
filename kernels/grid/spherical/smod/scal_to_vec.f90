submodule (SphericalHarmonics) scal_to_vec
  implicit none; contains
  
  module pure subroutine scal2vecscal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: ncr, crpadding, ncjm, cjmpadding
    complex(kind=dbl),    intent(inout) :: cr(ncr,*)
    complex(kind=dbl),    intent(inout) :: cjm(ncjm,*)
    integer                             :: i, j, m, mj, mj1, mj2, ijm
    complex(kind=dbl)                   :: cr12
    
    do concurrent ( mj = 1:this%jms2 )
      cr12               = ( +cr(crpadding,mj) + cr(crpadding+1,mj) * cunit ) * sq2_1
      cr(crpadding+1,mj) = ( -cr(crpadding,mj) + cr(crpadding+1,mj) * cunit ) * sq2_1
      cr(crpadding  ,mj) = cr12
    end do
      
    j = 0
      m = 0
        ijm = 1
        mj  = m*this%jmax3-m*(m+1)/2+j+1
        mj2 = mj + this%jmax + 3 - m
        
        cjm(cjmpadding+2,ijm) =      cr(crpadding  ,mj2 )   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                            &        cr(crpadding+2,mj+1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                            & conjg( cr(crpadding  ,mj2 ) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(cjmpadding+2,ijm)%im = zero
        
    do j = 1, this%jmax
      m = 0
        ijm = ijm+1
        mj  = m*this%jmax3-m*(m+1)/2+j+1
        mj2 = mj + this%jmax + 2 - m
        
        cjm(cjmpadding  ,ijm) =        cr(crpadding  ,mj2-1)   * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                              &        cr(crpadding+2,mj -1)   * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                              & conjg( cr(crpadding  ,mj2-1) ) * cleb1_fn(j-1,m-1,1,+1,j,m) ; cjm(cjmpadding,ijm)%im = zero
        cjm(cjmpadding+1,ijm) =        cr(crpadding  ,mj2  )   * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                              &        cr(crpadding+2,mj   )   * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                              & conjg( cr(crpadding  ,mj2  ) ) * cleb1_fn(j  ,m-1,1,+1,j,m) ; cjm(cjmpadding+1,ijm)%re = zero
        cjm(cjmpadding+2,ijm) =        cr(crpadding  ,mj2+1)   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                              &        cr(crpadding+2,mj +1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                              & conjg( cr(crpadding  ,mj2+1) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(cjmpadding+2,ijm)%im = zero
      
      do m = 1, j
        ijm = ijm+1
        mj  = m*this%jmax3-m*(m+1)/2+j+1
        mj1 = mj - this%jmax - 3 + m
        mj2 = mj + this%jmax + 2 - m
        
        cjm(cjmpadding  ,ijm) = cr(crpadding  ,mj2-1) * cleb1_fn(j-1,m+1,1,-1,j,m) + &
                              & cr(crpadding+2,mj -1) * cleb1_fn(j-1,m+0,1, 0,j,m) + &
                              & cr(crpadding+1,mj1-1) * cleb1_fn(j-1,m-1,1,+1,j,m)
        cjm(cjmpadding+1,ijm) = cr(crpadding  ,mj2  ) * cleb1_fn(j  ,m+1,1,-1,j,m) + &
                              & cr(crpadding+2,mj   ) * cleb1_fn(j  ,m+0,1, 0,j,m) + &
                              & cr(crpadding+1,mj1  ) * cleb1_fn(j  ,m-1,1,+1,j,m)
        cjm(cjmpadding+2,ijm) = cr(crpadding  ,mj2+1) * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                              & cr(crpadding+2,mj +1) * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                              & cr(crpadding+1,mj1+1) * cleb1_fn(j+1,m-1,1,+1,j,m)
      end do
    end do
    
  end subroutine scal2vecscal_mj_to_jm_sub
  
end submodule scal_to_vec