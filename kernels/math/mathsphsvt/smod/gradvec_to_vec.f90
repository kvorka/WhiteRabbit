submodule (Sphsvt) gradvec_to_vec
  implicit none; contains
  
  module pure subroutine gradvec2vec_jmlk_to_jml_sub(this, ri, v, dv_r, cab, ncab, cabpadding)
    class(T_sphsvt),   intent(in)    :: this
    integer,           intent(in)    :: cabpadding, ncab
    real(kind=dbl),    intent(in)    :: ri
    complex(kind=dbl), intent(in)    :: v(*), dv_r(*)
    complex(kind=dbl), intent(inout) :: cab(ncab,*)
    integer                          :: j, m, l, lmj, ijml, i1
    real(kind=dbl)                   :: cleb, fac1, fac2, fac3, fac4
    complex(kind=dbl), allocatable   :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(2), sum2(2), sum3(2) )
    
    do j = 0, this%jmax1
      fac1 = sqrt((j  )/(2*j+one))
      fac2 = sqrt((j+1)/(2*j+one))
      fac3 = (j+1)/ri
      fac4 = (j  )/ri
      
      m = 0
        sum1(1) = czero; sum1(2) = czero
        sum2(1) = czero; sum2(2) = czero
        sum3(1) = czero; sum3(2) = czero
        
        do l = abs(j-1), min(this%jmax, j+1)
          lmj = 3*(l*(l+1)/2+m+1)+j-l
          
          cleb = cleb1_fn(j,m,1, 0,l,m+0)
            sum3(1) = sum3(1) + ( dv_r(lmj-3) + fac3 * v(lmj-3) ) * cleb
            sum3(2) = sum3(2) - ( dv_r(lmj-3) - fac4 * v(lmj-3) ) * cleb
          
          cleb = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
            sum1(1) = sum1(1) + conjg( dv_r(lmj) + fac3 * v(lmj) ) * cleb
            sum1(2) = sum1(2) - conjg( dv_r(lmj) - fac4 * v(lmj) ) * cleb
          
          cleb = cleb1_fn(j,m,1,+1,l,m+1)
            sum2(1) = sum2(1) + ( dv_r(lmj) + fac3 * v(lmj) ) * cleb
            sum2(2) = sum2(2) - ( dv_r(lmj) - fac4 * v(lmj) ) * cleb
        end do
        
        ijml = 3*(j*(j+1)/2+m)+abs(j-1)-j
          cab(cabpadding  ,ijml) =         ( +sum1(1) - sum2(1) ) * sq2_1 * fac1
          cab(cabpadding+1,ijml) = cunit * ( -sum1(1) - sum2(1) ) * sq2_1 * fac1
          cab(cabpadding+2,ijml) =         ( +sum3(1)           )         * fac1
        
        ijml = 3*(j*(j+1)/2+m)+(j+1)-j
          cab(cabpadding  ,ijml) =         ( +sum1(2) - sum2(2) ) * sq2_1 * fac2
          cab(cabpadding+1,ijml) = cunit * ( -sum1(2) - sum2(2) ) * sq2_1 * fac2
          cab(cabpadding+2,ijml) =         ( +sum3(2)           )         * fac2
      
      do m = 1, j
        sum1(1) = czero; sum1(2) = czero
        sum2(1) = czero; sum2(2) = czero
        sum3(1) = czero; sum3(2) = czero
        
        do l = j-1, min(this%jmax, j+1)
          lmj = 3*(l*(l+1)/2+m-1)+j-l
          
          !every time
            cleb = cleb1_fn(j,m,1,-1,l,m-1)
              sum1(1) = sum1(1) + ( dv_r(lmj) + fac3 * v(lmj) ) * cleb
              sum1(2) = sum1(2) - ( dv_r(lmj) - fac4 * v(lmj) ) * cleb
            
          if ( l > m-1 ) then
            cleb = cleb1_fn(j,m,1, 0,l,m+0)
              sum3(1) = sum3(1) + ( dv_r(lmj+3) + fac3 * v(lmj+3) ) * cleb
              sum3(2) = sum3(2) - ( dv_r(lmj+3) - fac4 * v(lmj+3) ) * cleb
          end if
            
          if ( l > m ) then
            cleb = cleb1_fn(j,m,1,+1,l,m+1)
              sum2(1) = sum2(1) + ( dv_r(lmj+6) + fac3 * v(lmj+6) ) * cleb
              sum2(2) = sum2(2) - ( dv_r(lmj+6) - fac4 * v(lmj+6) ) * cleb
          end if
        end do
        
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
  
end submodule gradvec_to_vec