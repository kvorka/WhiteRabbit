module sph_norms
  use sph_indexing
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
      
      !sum over orders is split to m == 0 and higher orders (these are multiplied by factor 2)
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
      
      !sum over orders is split to m == 0 and higher orders (these are multiplied by factor 2)
      vp = vp +     sum( c2r_fn( cajml(j0-1:j0    +1) * conjg(cbjml(j0-1:j0    +1)) ) ) + &
              & 2 * sum( c2r_fn( cajml(j0+2:j0+3*j+1) * conjg(cbjml(j0+3:j0+3*j+1)) ) )
    end do
    
  end function dotproduct_fn
  
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
    integer                       :: j
    
    tnorm_fn = abs( cajml2(1) )**2 + sum( abs( cajml2(4:6) )**2 ) + 2 * sum( abs( cajml2(9:11) )**2 )
    
    do j = 2, np
      tnorm_fn = tnorm_fn +     sum( abs( cajml2(jml2(j,0,-2):jml2(j,0,+2)) )**2 ) + &
               &            2 * sum( abs( cajml2(jml2(j,1,-2):jml2(j,j,+2)) )**2 )
    end do
    
    tnorm_fn = sqrt( tnorm_fn )
    
  end function tnorm_fn
  
end module sph_norms