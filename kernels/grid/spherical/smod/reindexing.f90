submodule (SphericalHarmonics) reindexing
  implicit none; contains
  
  module pure subroutine vec2vec_jml_to_jml_sub(this, cjml, cab, ncab, cabpadding)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: ncab, cabpadding
    complex(kind=dbl),    intent(in)    :: cjml(*)
    complex(kind=dbl),    intent(inout) :: cab(ncab,*)
    integer                             :: ijm
    
    do concurrent( ijm = 1:this%jmv )
      cab(cabpadding,ijm) = cjml(ijm)
    end do
    
  end subroutine vec2vec_jml_to_jml_sub
  
  module pure subroutine scal2scal_jm_to_mj_sub(this, cjm, cab, ncab, cabpadding)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: ncab, cabpadding
    complex(kind=dbl),    intent(in)    :: cjm(*)
    complex(kind=dbl),    intent(inout) :: cab(ncab,*)
    integer                             :: m, j
    
    do m = 0, this%jmax
      do j = m, this%jmax
        cab(cabpadding, m*(this%jmax+3)-m*(m+1)/2+j+1) = cjm(j*(j+1)/2+m+1)
      end do
    end do
    
  end subroutine scal2scal_jm_to_mj_sub
  
  module pure subroutine vec2scal_jml_to_mj_sub(this, cab, ncab, cc)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: ncab
    complex(kind=dbl),    intent(in)  :: cab(ncab,*)
    complex(kind=dbl),    intent(out) :: cc(3,ncab,*)
    integer                           :: j, m, l, mj, lmj, i1
    real(kind=dbl)                    :: cleb
    complex(kind=dbl), allocatable    :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(ncab), sum2(ncab), sum3(ncab) )
    
    do m = 0, this%maxj
      do j = m, this%maxj
        do concurrent ( i1 = 1:ncab )
          sum1(i1) = czero
          sum2(i1) = czero
          sum3(i1) = czero
        end do
        
        if (m == 0) then
          do l = abs(j-1), min(this%jmax+1, j+1)
            lmj = 3*(l*(l+1)/2+m+1)+j-l
            
            cleb = cleb1_fn(j,m,1,0,l,m)
              do concurrent ( i1 = 1:ncab )
                sum3(i1) = sum3(i1) + cab(i1,lmj-3) * cleb
              end do
            
            cleb = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
              do concurrent ( i1 = 1:ncab )
                sum1(i1) = sum1(i1) + conjg( cab(i1,lmj) ) * cleb
              end do
            
            cleb = cleb1_fn(j,m,1,+1,l,m+1)
              do concurrent ( i1 = 1:ncab )
                sum2(i1) = sum2(i1) + cab(i1,lmj) * cleb
              end do
          end do
        else
          do l = abs(j-1), min(this%jmax+1, j+1)
            lmj = 3*(l*(l+1)/2+m-1)+j-l
            
            if (l > m) then
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                do concurrent ( i1 = 1:ncab )
                  sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
                end do
              
              cleb = cleb1_fn(j,m,1,0,l,m)
                do concurrent ( i1 = 1:ncab )
                  sum3(i1) = sum3(i1) + cab(i1,lmj+3) * cleb
                end do
              
              cleb = cleb1_fn(j,m,1,+1,l,m+1)
                do concurrent ( i1 = 1:ncab )
                  sum2(i1) = sum2(i1) + cab(i1,lmj+6) * cleb
                end do
                
            else if (l > m-1) then
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                do concurrent ( i1 = 1:ncab )
                  sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
                end do
              
              cleb = cleb1_fn(j,m,1,0,l,m)
                do concurrent ( i1 = 1:ncab )
                  sum3(i1) = sum3(i1) + cab(i1,lmj+3) * cleb
                end do
            
            else
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                do concurrent (i1 = 1:ncab)
                  sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
                end do
            end if
          end do
        end if
        
        mj = m*(this%jmax+3)-m*(m+1)/2+j+1
          do concurrent ( i1 = 1:ncab )
            cc(1,i1,mj) =         ( +sum1(i1) - sum2(i1) ) * sq2_1
            cc(2,i1,mj) = cunit * ( -sum1(i1) - sum2(i1) ) * sq2_1
            cc(3,i1,mj) =           +sum3(i1)
          end do
      end do
    end do
    
    deallocate( sum1, sum2, sum3 )
    
  end subroutine vec2scal_jml_to_mj_sub
  
  module pure subroutine gradvec2vec_jmlk_to_jml_sub(this, ri, v, dv_r, cab, ncab, cabpadding)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: cabpadding, ncab
    real(kind=dbl),       intent(in)    :: ri
    complex(kind=dbl),    intent(in)    :: v(*), dv_r(*)
    complex(kind=dbl),    intent(inout) :: cab(ncab,*)
    integer                             :: j, m, l, lmj, ijml, i1
    real(kind=dbl)                      :: cleb, fac1, fac2
    complex(kind=dbl), allocatable      :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(2), sum2(2), sum3(2) )
    
    do j = 0, this%jmax+1
      fac1 = sqrt((j  )/(2*j+1._dbl))
      fac2 = sqrt((j+1)/(2*j+1._dbl))
      
      do m = 0, j
        do concurrent ( i1 = 1:2 )
          sum1(i1) = czero
          sum2(i1) = czero
          sum3(i1) = czero
        end do
        
        if (m == 0) then
          do l = abs(j-1), min(this%jmax, j+1)
            lmj = 3*(l*(l+1)/2+m+1)+j-l
            
            cleb = cleb1_fn(j,m,1, 0,l,m+0)
              sum3(1) = sum3(1) + ( dv_r(lmj-3) + (j+1) * v(lmj-3) / ri ) * cleb
              sum3(2) = sum3(2) - ( dv_r(lmj-3) - (j  ) * v(lmj-3) / ri ) * cleb
            
            cleb = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
              sum1(1) = sum1(1) + conjg( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
              sum1(2) = sum1(2) - conjg( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
            
            cleb = cleb1_fn(j,m,1,+1,l,m+1)
              sum2(1) = sum2(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
              sum2(2) = sum2(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
          end do
        
        else
          do l = max(abs(m-1), j-1), min(this%jmax, j+1)
            lmj = 3*(l*(l+1)/2+m-1)+j-l
            
            if (l > m) then
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
              
              cleb = cleb1_fn(j,m,1, 0,l,m+0)
                sum3(1) = sum3(1) + ( dv_r(lmj+3) + (j+1) * v(lmj+3) / ri ) * cleb
                sum3(2) = sum3(2) - ( dv_r(lmj+3) - (j  ) * v(lmj+3) / ri ) * cleb
              
              cleb = cleb1_fn(j,m,1,+1,l,m+1)
                sum2(1) = sum2(1) + ( dv_r(lmj+6) + (j+1) * v(lmj+6) / ri ) * cleb
                sum2(2) = sum2(2) - ( dv_r(lmj+6) - (j  ) * v(lmj+6) / ri ) * cleb
            
            else if (l > m-1) then
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
              
              cleb = cleb1_fn(j,m,1, 0,l,m+0)
                sum3(1) = sum3(1) + ( dv_r(lmj+3) + (j+1) * v(lmj+3) / ri ) * cleb
                sum3(2) = sum3(2) - ( dv_r(lmj+3) - (j  ) * v(lmj+3) / ri ) * cleb
              
            else
              cleb = cleb1_fn(j,m,1,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
            end if
          end do
        end if
        
        ijml = 3*(j*(j+1)/2+m)+abs(j-1)-j
          cab(cabpadding  ,ijml) =         ( +sum1(1) - sum2(1) ) * sq2_1 * fac1
          cab(cabpadding+1,ijml) = cunit * ( -sum1(1) - sum2(1) ) * sq2_1 * fac1
          cab(cabpadding+2,ijml) =         ( +sum3(1)           )         * fac1
        
        ijml = 3*(j*(j+1)/2+m)+(j+1)-j
          cab(cabpadding  ,ijml) =         ( +sum1(2) - sum2(2) ) * sq2_1 * fac2
          cab(cabpadding+1,ijml) = cunit * ( -sum1(2) - sum2(2) ) * sq2_1 * fac2
          cab(cabpadding+2,ijml) =         ( +sum3(2)           )         * fac2
      end do
    end do
    
    deallocate( sum1, sum2, sum3 )
    
  end subroutine gradvec2vec_jmlk_to_jml_sub
  
  module pure subroutine devtens2scal_jml2_to_mj_sub(this, ctjml2, cr, ncr, crpadding)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: ncr, crpadding
    complex(kind=dbl),    intent(in)    :: ctjml2(*)
    complex(kind=dbl),    intent(inout) :: cr(ncr,*)
    integer                             :: j, m, l, k, lm
    complex(kind=dbl)                   :: csum, cpom
    
    do k = -2, 2
      do l = 0, this%maxj
        do m = 0, l
          csum = czero
          
          do j = max(abs(m+k),abs(l-2)), min(this%jmax,l+2)
            if ( m >= -k ) then
              csum = csum + cleb2_fn(l,m,2,k,j,m+k) * ctjml2(5*(j*(j+1)/2+abs(m+k))+l-j-1)
            else
              csum = csum + cleb2_fn(l,m,2,k,j,m+k) * (-1)**(j+m+k+l) * conjg( ctjml2(5*(j*(j+1)/2+abs(m+k))+l-j-1) )
            end if
          end do
          
          lm = m*(this%jmax+3)-m*(m+1)/2+l+1
          
          if ( k <= 0 ) then
            cr(crpadding+k+2,lm) = csum
            
          else if ( k == 1 ) then
            cpom               = cr(crpadding+1,lm)
            cr(crpadding+1,lm) =           cpom - csum
            cr(crpadding+3,lm) = cunit * ( cpom + csum )
            
          else
            cpom               = cr(crpadding,lm)
            cr(crpadding  ,lm) =            cpom + csum 
            cr(crpadding+4,lm) = cunit * ( -cpom + csum )
          end if
          
        end do
      end do
    end do
    
  end subroutine devtens2scal_jml2_to_mj_sub
  
  module pure subroutine scal2scal_mj_to_jm_sub(this, cr, ncr, crpadding, cjm, ncjm, cjmpadding)
    class(T_lateralGrid), intent(in)    :: this
    integer,              intent(in)    :: ncr, ncjm, crpadding, cjmpadding
    complex(kind=dbl),    intent(in)    :: cr(ncr,*)
    complex(kind=dbl),    intent(inout) :: cjm(ncjm,*)
    integer                             :: j, m, ijm, imj
    
    do j = 0, this%jmax
      m = 0
        ijm = j*(j+1)/2+m+1
        imj = m*(this%jmax+3)-m*(m+1)/2+j+1
          cjm(cjmpadding,ijm)%re = cr(crpadding,imj)%re
          cjm(cjmpadding,ijm)%im = zero
      
      do m = 1, j
        ijm = j*(j+1)/2+m+1
        imj = m*(this%jmax+3)-m*(m+1)/2+j+1
          cjm(cjmpadding,ijm) = cr(crpadding,imj)
      end do
    end do
    
  end subroutine scal2scal_mj_to_jm_sub
  
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
        mj  = m*(this%jmax+3)-m*(m+1)/2+j+1
        mj2 = mj + this%jmax + 3 - m
        
        cjm(cjmpadding+2,ijm) =      cr(crpadding  ,mj2 )   * cleb1_fn(j+1,m+1,1,-1,j,m) + &
                            &        cr(crpadding+2,mj+1)   * cleb1_fn(j+1,m+0,1, 0,j,m) + &
                            & conjg( cr(crpadding  ,mj2 ) ) * cleb1_fn(j+1,m-1,1,+1,j,m) ; cjm(cjmpadding+2,ijm)%im = zero
        
    do j = 1, this%jmax
      m = 0
        ijm = ijm+1
        mj  = m*(this%jmax+3)-m*(m+1)/2+j+1
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
        mj  = m*(this%jmax+3)-m*(m+1)/2+j+1
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
  
  module pure subroutine scal2devtens_mj_to_jml2_sub(this, cr, ncr, crpadding, ctjml2)
    class(T_lateralGrid), intent(in) :: this
    integer,              intent(in)  :: ncr, crpadding
    complex(kind=dbl),    intent(in)  :: cr(ncr,*)
    complex(kind=dbl),    intent(out) :: ctjml2(*)
    integer                           :: j, m, l, k, i1, lm
    complex(kind=dbl)                 :: csum
    complex(kind=dbl), allocatable    :: caux(:)
    
    allocate( caux(0:4) )
    
    do j = 0, this%jmax
      do m = 0, j
        do l = abs(j-2), j+2
          csum = czero
          
          do k = -2, 2
            lm = abs(m-k)*(this%jmax+3)-abs(m-k)*(abs(m-k)+1)/2+l+1
            
            if (m < k) then
              do concurrent ( i1 = 0:4 )
                caux(i1) = (-1)**(k-m) * conjg( cr(crpadding+i1,lm) )
              end do
            else
              do concurrent ( i1 = 0:4 )
                caux(i1) = cr(crpadding+i1,lm)
              end do
            end if
            
            select case (k)
              case (-2)
                csum = csum + ( caux(0) + cunit * caux(4) ) * cleb2_fn(l,m+2,2,-2,j,m)
              case (-1)
                csum = csum + ( caux(1) - cunit * caux(3) ) * cleb2_fn(l,m+1,2,-1,j,m)
              case ( 0)
                csum = csum + ( caux(2) +         caux(2) ) * cleb2_fn(l,m  ,2, 0,j,m)
              case (+1)
                csum = csum - ( caux(1) + cunit * caux(3) ) * cleb2_fn(l,m-1,2,+1,j,m)
              case (+2)
                csum = csum + ( caux(0) - cunit * caux(4) ) * cleb2_fn(l,m-2,2,+2,j,m)
            end select
          end do
          
          ctjml2(5*(j*(j+1)/2+m)+l-j-1) = csum / 2
        end do
      end do
    end do
    
    deallocate( caux )
    
  end subroutine scal2devtens_mj_to_jml2_sub
  
end submodule reindexing