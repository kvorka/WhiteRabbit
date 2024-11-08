submodule (sphsvt) gradvec_to_vec
  implicit none; contains
  
  module procedure gradvec2vec_jmlk_to_jml_sub
    integer                        :: j, m, l, lmj, ijml, i1
    real(kind=dbl)                 :: cg, fac1, fac2, fac3, fac4
    complex(kind=dbl), allocatable :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(2), sum2(2), sum3(2) )
    
    do j = 0, this%jmax1
      fac1 = sqrt((j  )/(2*j+one))
      fac2 = sqrt((j+1)/(2*j+one))
      fac3 = (j+1)/ri
      fac4 = (j  )/ri
      
      m = 0
      call zero_carray_sub( 2, sum1(1) )
      call zero_carray_sub( 2, sum2(1) )
      call zero_carray_sub( 2, sum3(1) )
        
        do l = abs(j-1), min(this%jmax, j+1)
          lmj = 3*(l*(l+1)/2+m+1)+j-l
          
          cg = cleb1_fn(j,m,1, 0,l,m+0)
            sum3(1) = sum3(1) + ( dv_r(lmj-3) + fac3 * v(lmj-3) ) * cg
            sum3(2) = sum3(2) - ( dv_r(lmj-3) - fac4 * v(lmj-3) ) * cg
          
          cg = cleb1_fn(j,m,1,-1,l,m-1) * (-1)**(j+l)
            sum1(1) = sum1(1) + conjg( dv_r(lmj) + fac3 * v(lmj) ) * cg
            sum1(2) = sum1(2) - conjg( dv_r(lmj) - fac4 * v(lmj) ) * cg
          
          cg = cleb1_fn(j,m,1,+1,l,m+1)
            sum2(1) = sum2(1) + ( dv_r(lmj) + fac3 * v(lmj) ) * cg
            sum2(2) = sum2(2) - ( dv_r(lmj) - fac4 * v(lmj) ) * cg
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
        call zero_carray_sub( 2, sum1(1) )
        call zero_carray_sub( 2, sum2(1) )
        call zero_carray_sub( 2, sum3(1) )
        
        do l = j-1, min(this%jmax, j+1)
          lmj = 3*(l*(l+1)/2+m-1)+j-l
          
          !every time
            cg = cleb1_fn(j,m,1,-1,l,m-1)
              sum1(1) = sum1(1) + ( dv_r(lmj) + fac3 * v(lmj) ) * cg
              sum1(2) = sum1(2) - ( dv_r(lmj) - fac4 * v(lmj) ) * cg
            
          if ( l > m-1 ) then
            cg = cleb1_fn(j,m,1, 0,l,m+0)
              sum3(1) = sum3(1) + ( dv_r(lmj+3) + fac3 * v(lmj+3) ) * cg
              sum3(2) = sum3(2) - ( dv_r(lmj+3) - fac4 * v(lmj+3) ) * cg
          end if
            
          if ( l > m ) then
            cg = cleb1_fn(j,m,1,+1,l,m+1)
              sum2(1) = sum2(1) + ( dv_r(lmj+6) + fac3 * v(lmj+6) ) * cg
              sum2(2) = sum2(2) - ( dv_r(lmj+6) - fac4 * v(lmj+6) ) * cg
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
    
  end procedure gradvec2vec_jmlk_to_jml_sub
  
end submodule gradvec_to_vec