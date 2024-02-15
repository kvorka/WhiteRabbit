module sph_vectors
  use Clebsch_Gordan
  implicit none; public; contains
  
  pure subroutine gradvec_to_vectors_sub(jmax, padding, nvec, ri, v, dv_r, cab)
    integer,           intent(in)    :: jmax, padding, nvec
    real(kind=dbl),    intent(in)    :: ri
    complex(kind=dbl), intent(in)    :: v(*), dv_r(*)
    complex(kind=dbl), intent(inout) :: cab(*)
    integer                          :: j, m, l, lmj, ijml, i1
    real(kind=dbl)                   :: cleb, fac1, fac2
    complex(kind=dbl), allocatable   :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(2), sum2(2), sum3(2) )
    
    do j = 0, jmax+1
      fac1 = sqrt((j  )/(2*j+1._dbl))
      fac2 = sqrt((j+1)/(2*j+1._dbl))
      
      do m = 0, j
        do concurrent ( i1 = 1:2 )
          sum1(i1) = czero
          sum2(i1) = czero
          sum3(i1) = czero
        end do
        
        if (m == 0) then
          do l = abs(j-1), min(jmax, j+1)
            lmj = 3*(l*(l+1)/2+m+1)+j-l
            
            cleb = cleb_fn(j-l,0,l,m)
              sum3(1) = sum3(1) + ( dv_r(lmj-3) + (j+1) * v(lmj-3) / ri ) * cleb
              sum3(2) = sum3(2) - ( dv_r(lmj-3) - (j  ) * v(lmj-3) / ri ) * cleb
            
            cleb = cleb_fn(j-l,-1,l,m-1) * (-1)**(j+l)
              sum1(1) = sum1(1) + conjg( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
              sum1(2) = sum1(2) - conjg( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
            
            cleb = cleb_fn(j-l,+1,l,m+1)
              sum2(1) = sum2(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
              sum2(2) = sum2(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
          end do
        
        else
          do l = max(abs(m-1), j-1), min(jmax, j+1)
            lmj = 3*(l*(l+1)/2+m-1)+j-l
            
            if (l > m) then
              cleb = cleb_fn(j-l,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
              
              cleb = cleb_fn(j-l,0,l,m)
                sum3(1) = sum3(1) + ( dv_r(lmj+3) + (j+1) * v(lmj+3) / ri ) * cleb
                sum3(2) = sum3(2) - ( dv_r(lmj+3) - (j  ) * v(lmj+3) / ri ) * cleb
              
              cleb = cleb_fn(j-l,+1,l,m+1)
                sum2(1) = sum2(1) + ( dv_r(lmj+6) + (j+1) * v(lmj+6) / ri ) * cleb
                sum2(2) = sum2(2) - ( dv_r(lmj+6) - (j  ) * v(lmj+6) / ri ) * cleb
            
            else if (l > m-1) then
              cleb = cleb_fn(j-l,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
              
              cleb = cleb_fn(j-l,0,l,m)
                sum3(1) = sum3(1) + ( dv_r(lmj+3) + (j+1) * v(lmj+3) / ri ) * cleb
                sum3(2) = sum3(2) - ( dv_r(lmj+3) - (j  ) * v(lmj+3) / ri ) * cleb
              
            else
              cleb = cleb_fn(j-l,-1,l,m-1)
                sum1(1) = sum1(1) + ( dv_r(lmj) + (j+1) * v(lmj) / ri ) * cleb
                sum1(2) = sum1(2) - ( dv_r(lmj) - (j  ) * v(lmj) / ri ) * cleb
            end if
          end do
        end if
        
        ijml = 3*(j*(j+1)/2+m)+abs(j-1)-j-1
          cab(padding  +nvec*ijml) =         ( +sum1(1) - sum2(1) ) * sq2_1 * fac1
          cab(padding+1+nvec*ijml) = cunit * ( -sum1(1) - sum2(1) ) * sq2_1 * fac1
          cab(padding+2+nvec*ijml) =         ( +sum3(1)           )         * fac1
        
        ijml = 3*(j*(j+1)/2+m)+(j+1)-j-1
          cab(padding  +nvec*ijml) =         ( +sum1(2) - sum2(2) ) * sq2_1 * fac2
          cab(padding+1+nvec*ijml) = cunit * ( -sum1(2) - sum2(2) ) * sq2_1 * fac2
          cab(padding+2+nvec*ijml) =         ( +sum3(2)           )         * fac2
      end do
    end do
    
    deallocate( sum1, sum2, sum3 )
    
  end subroutine gradvec_to_vectors_sub
  
  pure subroutine vectors_to_scalars_sub(jmax, nvec, cab, cc)
    integer,           intent(in)  :: jmax, nvec
    complex(kind=dbl), intent(in)  :: cab(*)
    complex(kind=dbl), intent(out) :: cc(*)
    integer                        :: j, m, l, mj, lmj, i1, i2
    real(kind=dbl)                 :: cleb
    complex(kind=dbl), allocatable :: sum1(:), sum2(:), sum3(:)
    
    allocate( sum1(nvec), sum2(nvec), sum3(nvec) )
    
    do m = 0, jmax+2
      do j = m, jmax+2
        do concurrent ( i1 = 1:nvec )
          sum1(i1) = czero
          sum2(i1) = czero
          sum3(i1) = czero
        end do
        
        if (m == 0) then
          do l = abs(j-1), min(jmax+1, j+1)
            lmj = nvec*(3*(l*(l+1)/2+m+1)+j-l-1)
            
            cleb = cleb_fn(j-l,0,l,m)
              do concurrent ( i1 = 1:nvec )
                sum3(i1) = sum3(i1) + cab(i1-3*nvec+lmj) * cleb
              end do
            
            cleb = cleb_fn(j-l,-1,l,m-1) * (-1)**(j+l)
              do concurrent ( i1 = 1:nvec )
                sum1(i1) = sum1(i1) + conjg( cab(i1+lmj) ) * cleb
              end do
            
            cleb = cleb_fn(j-l,+1,l,m+1)
              do concurrent ( i1 = 1:nvec )
                sum2(i1) = sum2(i1) + cab(i1+lmj) * cleb
              end do
          end do
        else
          do l = abs(j-1), min(jmax+1, j+1)
            lmj = nvec*(3*(l*(l+1)/2+m-1)+j-l-1)
            
            if (l > m) then
              cleb = cleb_fn(j-l,-1,l,m-1)
                do concurrent ( i1 = 1:nvec )
                  sum1(i1) = sum1(i1) + cab(i1+lmj) * cleb
                end do
              
              cleb = cleb_fn(j-l,0,l,m)
                do concurrent ( i1 = 1:nvec )
                  sum3(i1) = sum3(i1) + cab(i1+3*nvec+lmj) * cleb
                end do
              
              cleb = cleb_fn(j-l,+1,l,m+1)
                do concurrent ( i1 = 1:nvec )
                  sum2(i1) = sum2(i1) + cab(i1+6*nvec+lmj) * cleb
                end do
                
            else if (l > m-1) then
              cleb = cleb_fn(j-l,-1,l,m-1)
                do concurrent ( i1 = 1:nvec )
                  sum1(i1) = sum1(i1) + cab(i1+lmj) * cleb
                end do
              
              cleb = cleb_fn(j-l,0,l,m)
                do concurrent ( i1 = 1:nvec )
                  sum3(i1) = sum3(i1) + cab(i1+3*nvec+lmj) * cleb
                end do
            
            else
              cleb = cleb_fn(j-l,-1,l,m-1)
                do concurrent (i1 = 1:nvec)
                  sum1(i1) = sum1(i1) + cab(i1+lmj) * cleb
                end do
            end if
          end do
        end if
        
        mj = m*(jmax+3)-m*(m+1)/2+j
          do concurrent ( i1 = 1:nvec )
            i2 = 3*(i1-1)+1
            
            cc(i2  +3*nvec*mj) =         ( +sum1(i1) - sum2(i1) ) * sq2_1
            cc(i2+1+3*nvec*mj) = cunit * ( -sum1(i1) - sum2(i1) ) * sq2_1
            cc(i2+2+3*nvec*mj) =           +sum3(i1)
          end do
      end do
    end do
    
    deallocate( sum1, sum2, sum3 )
    
  end subroutine vectors_to_scalars_sub
  
  pure subroutine cartesian_to_cyclic_sub(padding, nvec, length, cr)
    integer,           intent(in)    :: padding, nvec, length
    complex(kind=dbl), intent(inout) :: cr(*)
    integer                          :: mj
    complex(kind=dbl)                :: cr12
    
    do concurrent ( mj = 0:length-1 )
      cr12                  = ( +cr(padding+nvec*mj) + cr(padding+1+nvec*mj) * cunit ) * sq2_1
      cr(padding+1+nvec*mj) = ( -cr(padding+nvec*mj) + cr(padding+1+nvec*mj) * cunit ) * sq2_1
      cr(padding+  nvec*mj) = cr12
    end do
    
  end subroutine cartesian_to_cyclic_sub
  
end module sph_vectors
