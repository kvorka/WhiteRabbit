submodule (sphsvt) vec_to_scal
  implicit none; contains
  
  module procedure vec2scal_jml_to_mj_sub
    integer                        :: j, m, l, mj, lmj, i1, indx
    real(kind=dbl)                 :: cleb
    complex(kind=dbl), allocatable :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(ncab), sum2(ncab), sum3(ncab) )
    
    m = 0
      do j = m, this%jmax2
        do concurrent ( i1 = 1:ncab )
          sum1(i1) = czero
          sum2(i1) = czero
          sum3(i1) = czero
        end do
        
        do l = abs(j-1), min(this%jmax1, j+1)
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
        
        mj = m*this%jmax3-m*(m+1)/2+j+1
          do concurrent ( i1 = 1:ncab )
            indx = 3*(i1-1)+ccpadding
              cc(indx  ,mj) =         ( +sum1(i1) - sum2(i1) ) * sq2_1
              cc(indx+1,mj) = cunit * ( -sum1(i1) - sum2(i1) ) * sq2_1
              cc(indx+2,mj) =           +sum3(i1)
          end do
      end do
    
    do m = 1, this%jmax2
      do j = m, this%jmax2
        do concurrent ( i1 = 1:ncab )
          sum1(i1) = czero
          sum2(i1) = czero
          sum3(i1) = czero
        end do
        
        do l = j-1, min(this%jmax1, j+1)
          lmj = 3*(l*(l+1)/2+m-1)+j-l
          
          !every time
            cleb = cleb1_fn(j,m,1,-1,l,m-1)
              do concurrent ( i1 = 1:ncab )
                sum1(i1) = sum1(i1) + cab(i1,lmj) * cleb
              end do
          
          if ( l > m-1 ) then
            cleb = cleb1_fn(j,m,1,0,l,m)
              do concurrent ( i1 = 1:ncab )
                sum3(i1) = sum3(i1) + cab(i1,lmj+3) * cleb
              end do
          end if
          
          if ( l > m ) then
            cleb = cleb1_fn(j,m,1,+1,l,m+1)
              do concurrent ( i1 = 1:ncab )
                sum2(i1) = sum2(i1) + cab(i1,lmj+6) * cleb
              end do
          end if
        end do
        
        mj = m*this%jmax3-m*(m+1)/2+j+1
          do concurrent ( i1 = 1:ncab )
            indx = 3*(i1-1)+ccpadding
              cc(indx  ,mj) =         ( +sum1(i1) - sum2(i1) ) * sq2_1
              cc(indx+1,mj) = cunit * ( -sum1(i1) - sum2(i1) ) * sq2_1
              cc(indx+2,mj) =           +sum3(i1)
          end do
      end do
    end do
    
    deallocate( sum1, sum2, sum3 )
    
  end procedure vec2scal_jml_to_mj_sub
  
end submodule vec_to_scal