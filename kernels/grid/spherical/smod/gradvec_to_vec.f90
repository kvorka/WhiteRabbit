submodule (SphericalHarmonics) gradvec_to_vec
  implicit none; contains
  
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
    
    do j = 0, this%jmax1
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
  
end submodule gradvec_to_vec