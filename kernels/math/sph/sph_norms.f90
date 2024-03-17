module sph_norms
  use Conversions
  implicit none; public; contains
  
  pure function scalproduct_fn(np, cajm, cbjm) result(sp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(*), cbjm(*)
    integer                       :: j, j0
    real(kind=dbl)                :: sp
    
    !j == 0
    sp = c2r_fn( cajm(1) * conjg( cbjm(1) ) )
    
    !higher degrees
    do j = 1, np
      j0 = j*(j+1)/2+1
      
      sp = sp + c2r_fn( cajm(j0) * conjg(cbjm(j0)) ) + 2 * sum( c2r_fn( cajm(j0+1:j0+j) * conjg(cbjm(j0+1:j0+j)) ) )
    end do
    
  end function scalproduct_fn
  
  pure function dotproduct_fn(np, cajml, cbjml) result(vp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(*), cbjml(*)
    integer                       :: j, j0
    real(kind=dbl)                :: vp
    
    !j == 0
    vp = c2r_fn( cajml(1) * conjg(cbjml(1)) )
    
    !higher degrees
    do j = 1, np
      j0 = 3*(j*(j+1)/2)
      
      vp = vp +     sum( c2r_fn( cajml(j0-1:j0    +1) * conjg(cbjml(j0-1:j0    +1)) ) ) + &
              & 2 * sum( c2r_fn( cajml(j0+2:j0+3*j+1) * conjg(cbjml(j0+3:j0+3*j+1)) ) )
    end do
    
  end function dotproduct_fn
  
  pure function tensproduct_fn(np, cajml2, cbjml2) result(tp)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(:), cbjml2(:)
    integer                       :: j, j0
    real(kind=dbl)                :: tp
    
    !j == 0 and j == 1
    tp = c2r_fn( cajml2(1) * conjg(cbjml2(1)) ) +     sum( c2r_fn( cajml2(4: 6) * conjg(cbjml2(4: 6)) ) ) + &
                                                & 2 * sum( c2r_fn( cajml2(9:11) * conjg(cbjml2(9:11)) ) )
    
    !higher degrees
    do j = 2, np
      j0 = 5*(j*(j+1)/2)-3
      
      tp = tp +      sum( c2r_fn( cajml2(j0  :j0    +4) * conjg(cbjml2(j0  :j0    +4)) ) ) + &
               & 2 * sum( c2r_fn( cajml2(j0+5:j0+5*j+4) * conjg(cbjml2(j0+5:j0+5*j+4)) ) )
    end do
  
  end function tensproduct_fn
  
  pure real(kind=dbl) function snorm_fn(np, cajm)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(:)
    
    snorm_fn = sqrt( scalproduct_fn(np, cajm, cajm) )
    
  end function snorm_fn
  
  pure real(kind=dbl) function vnorm_fn(np, cajml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:)
    
    vnorm_fn = sqrt( dotproduct_fn(np, cajml, cajml) )
    
  end function vnorm_fn
  
  pure real(kind=dbl) function tnorm_fn(np, cajml2)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(:)
    
    tnorm_fn = sqrt( tensproduct_fn(np, cajml2, cajml2) )
    
  end function tnorm_fn
  
end module sph_norms