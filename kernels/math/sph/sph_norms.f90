module sph_norms
  use Conversions
  implicit none; public; contains
  
  pure function scalproduct_fn(np, cajm, cbjm) result(sp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(*), cbjm(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: sp
    
    !j == 0
    sp = c2r_fn( cajm(1) * conjg( cbjm(1) ) )
    
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
    vp = c2r_fn( cajml(1) * conjg(cbjml(1)) )
    
    !higher degrees
    do j = 1, np
      !m == 0
      indx = 3*(j*(j+1)/2)-1
      vp   = vp + sum( c2r_fn( cajml(indx:indx+2) * conjg(cbjml(indx:indx+2)) ) )
      
      do m = 1, j
        indx = 3*(j*(j+1)/2+m)-1
        vp   = vp + 2 * sum( c2r_fn( cajml(indx:indx+2) * conjg(cbjml(indx:indx+2)) ) )
      end do
    end do
    
  end function dotproduct_fn
  
  pure function tensproduct_fn(np, cajml2, cbjml2) result(tp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(*), cbjml2(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: tp
    
    !j == 0
    tp = c2r_fn( cajml2(1) * conjg(cbjml2(1)) ) 
    
    !j == 1
    tp = tp +     sum( c2r_fn( cajml2(4: 6) * conjg(cbjml2(4: 6)) ) ) + &
            & 2 * sum( c2r_fn( cajml2(9:11) * conjg(cbjml2(9:11)) ) )
    
    !higher degrees
    do j = 2, np
      !m == 0
      indx = 5*(j*(j+1)/2)-3
      tp   = tp + sum( c2r_fn( cajml2(indx:indx+4) * conjg(cbjml2(indx:indx+4)) ) )
      
      do m = 1, j
        indx = 5*(j*(j+1)/2+m)-3
        tp   = tp + 2 * sum( c2r_fn( cajml2(indx:indx+4) * conjg(cbjml2(indx:indx+4)) ) )
      end do
    end do
    
  end function tensproduct_fn
  
  pure function scalnorm2_fn(np, cajm) result(sp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: sp
    
    !j == 0
    sp = c2r_fn( cajm(1) * conjg( cajm(1) ) )
    
    !higher degrees
    do j = 1, np
      !m == 0
      indx = j*(j+1)/2+1
      sp   = sp + c2r_fn( cajm(indx) * conjg(cajm(indx)) )
      
      do m = 1, j
        indx = j*(j+1)/2+m+1
        sp   = sp + 2 * c2r_fn( cajm(indx) * conjg(cajm(indx)) )
      end do
    end do
    
  end function scalnorm2_fn
  
  pure function vectnorm2_fn(np, cajml) result(vp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: vp
    
    !j == 0
    vp = c2r_fn( cajml(1) * conjg(cajml(1)) )
    
    !higher degrees
    do j = 1, np
      !m == 0
      indx = 3*(j*(j+1)/2)-1
      vp   = vp + sum( c2r_fn( cajml(indx:indx+2) * conjg(cajml(indx:indx+2)) ) )
      
      do m = 1, j
        indx = 3*(j*(j+1)/2+m)-1
        vp   = vp + 2 * sum( c2r_fn( cajml(indx:indx+2) * conjg(cajml(indx:indx+2)) ) )
      end do
    end do
    
  end function vectnorm2_fn
  
  pure function tensnorm2_fn(np, cajml2) result(tp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(*)
    integer                       :: j, m, indx
    real(kind=dbl)                :: tp
    
    !j == 0
    tp = c2r_fn( cajml2(1) * conjg(cajml2(1)) ) 
    
    !j == 1
    tp = tp +     sum( c2r_fn( cajml2(4: 6) * conjg(cajml2(4: 6)) ) ) + &
            & 2 * sum( c2r_fn( cajml2(9:11) * conjg(cajml2(9:11)) ) )
    
    !higher degrees
    do j = 2, np
      !m == 0
      indx = 5*(j*(j+1)/2)-3
      tp   = tp + sum( c2r_fn( cajml2(indx:indx+4) * conjg(cajml2(indx:indx+4)) ) )
      
      do m = 1, j
        indx = 5*(j*(j+1)/2+m)-3
        tp   = tp + 2 * sum( c2r_fn( cajml2(indx:indx+4) * conjg(cajml2(indx:indx+4)) ) )
      end do
    end do
    
  end function tensnorm2_fn
  
end module sph_norms