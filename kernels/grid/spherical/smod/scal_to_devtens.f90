submodule (SphericalHarmonics) scal_to_devtens
  implicit none; contains
  
  module pure subroutine scal2devtens_mj_to_jml2_sub(this, cr, ncr, crpadding, ctjml2)
    class(T_lateralGrid), intent(in) :: this
    integer,              intent(in)  :: ncr, crpadding
    complex(kind=dbl),    intent(in)  :: cr(ncr,*)
    complex(kind=dbl),    intent(out) :: ctjml2(*)
    integer                           :: j, m, l, k, i1, lm, ijml2
    complex(kind=dbl)                 :: csum
    complex(kind=dbl), allocatable    :: caux(:)
    
    allocate( caux(0:4) )
    
    do j = 0, this%jmax
      m = 0
        do l = abs(j-2), j+2
          csum = czero
          
          do k = -2, 2
            lm = abs(k)*this%jmax3-abs(k)*(abs(k)+1)/2+l+1
            
            if (m < k) then
              do concurrent ( i1 = 0:4 )
                caux(i1) = (-1)**(k) * conjg( cr(crpadding+i1,lm) )
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
          
          ijml2 = 5*(j*(j+1)/2+m)+l-j-1
          select case (l-j)
            case (-2)
              ctjml2(ijml2)%re = csum%re / 2
              ctjml2(ijml2)%im = zero
            case (-1)
              ctjml2(ijml2)%re = zero
              ctjml2(ijml2)%im = csum%im / 2
            case (0)
              ctjml2(ijml2)%re = csum%re / 2
              ctjml2(ijml2)%im = zero
            case (+1)
              ctjml2(ijml2)%re = zero
              ctjml2(ijml2)%im = csum%im / 2
            case (+2)
              ctjml2(ijml2)%re = csum%re / 2
              ctjml2(ijml2)%im = zero
          end select
        end do
        
      do m = 1, j
        do l = abs(j-2), j+2
          csum = czero
          
          do k = -2, 2
            lm = abs(m-k)*this%jmax3-abs(m-k)*(abs(m-k)+1)/2+l+1
            
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
  
end submodule scal_to_devtens