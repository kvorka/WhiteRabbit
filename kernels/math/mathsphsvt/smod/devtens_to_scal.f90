submodule (Sphsvt) devtens_to_scal
  implicit none; contains
  
  module pure subroutine devtens2scal_jml2_to_mj_sub(this, ctjml2, cr, ncr, crpadding)
    class(T_sphsvt),   intent(in)    :: this
    integer,           intent(in)    :: ncr, crpadding
    complex(kind=dbl), intent(in)    :: ctjml2(*)
    complex(kind=dbl), intent(inout) :: cr(ncr,*)
    integer                          :: j, m, l, k, lm
    complex(kind=dbl)                :: csum, cpom
    
    do k = -2, 2
      do l = 0, this%jmax2
        do m = 0, l
          csum = czero
          
          do j = max(abs(m+k),abs(l-2)), min(this%jmax,l+2)
            if ( m >= -k ) then
              csum = csum + cleb2_fn(l,m,2,k,j,m+k) * ctjml2(5*(j*(j+1)/2+abs(m+k))+l-j-1)
            else
              csum = csum + cleb2_fn(l,m,2,k,j,m+k) * (-1)**(j+m+k+l) * conjg( ctjml2(5*(j*(j+1)/2+abs(m+k))+l-j-1) )
            end if
          end do
          
          lm = m*this%jmax3-m*(m+1)/2+l+1
          
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
  
end submodule devtens_to_scal