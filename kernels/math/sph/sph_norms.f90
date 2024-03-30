module sph_norms
  use Conversions
  implicit none; public; contains
  
  pure function scalproduct_fn(np, cajm, cbjm) result(sp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(*), cbjm(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: sp
    
    !j == 0
    sp = cajm(1)%re * cbjm(1)%re
    
    !higher degrees
    do j = 1, np
      !m == 0
      indx = j*(j+1)/2+1
      sp   = sp + c2r_fn( cajm(indx) * conjg(cbjm(indx)) )
      
      do m = 1, j
        indx = j*(j+1)/2+m+1
        sp   = sp + 2 * c2r_fn( cajm(indx) * conjg(cbjm(indx)) )
      end do
    end do
    
  end function scalproduct_fn
  
  pure function dotproduct_fn(np, cajml, cbjml) result(vp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(*), cbjml(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: vp
    
    !j == 0
    vp = cajml(1)%re * cbjml(1)%re
    
    !higher degrees
    do j = 1, np
      !m == 0
      indx = 3*(j*(j+1)/2)-1
      vp   = vp + c2r_fn( dot_product( cajml(indx:indx+2), cbjml(indx:indx+2) ) )
      
      do m = 1, j
        indx = 3*(j*(j+1)/2+m)-1
        vp   = vp + 2 * c2r_fn( dot_product( cajml(indx:indx+2), cbjml(indx:indx+2) ) )
      end do
    end do
    
  end function dotproduct_fn
  
  pure function tensproduct_fn(np, cajml2, cbjml2) result(tp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(*), cbjml2(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: tp
    
    !j == 0
    tp = cajml2(1)%re * cbjml2(1)%re
    
    !j == 1
    tp = tp +     c2r_fn( dot_product( cajml2(4: 6), cbjml2(4: 6) ) ) + &
            & 2 * c2r_fn( dot_product( cajml2(9:11), cbjml2(9:11) ) )
    
    !higher degrees
    do j = 2, np
      !m == 0
      indx = 5*(j*(j+1)/2)-3
      tp   = tp + c2r_fn( dot_product( cajml2(indx:indx+4), cbjml2(indx:indx+4) ) )
      
      do m = 1, j
        indx = 5*(j*(j+1)/2+m)-3
        tp   = tp + 2 * c2r_fn( dot_product( cajml2(indx:indx+4), cbjml2(indx:indx+4) ) )
      end do
    end do
    
  end function tensproduct_fn
  
end module sph_norms