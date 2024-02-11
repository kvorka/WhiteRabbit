module sph_norms
  use sph_indexing
  use Conversions
  implicit none; public; contains
  
  pure real(kind=dbl) function snorm_fn(np, cajm)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(:)
    integer                       :: ij
    
    ij = 0
      snorm_fn = abs( cajm(1) )**2
    
    do ij = 1, np
      snorm_fn = snorm_fn + abs( cajm(jm(ij,0)) )**2 + 2 * sum( abs( cajm(jm(ij,1):jm(ij,ij)) )**2 )
    end do
    
    snorm_fn = sqrt( snorm_fn )
    
  end function snorm_fn
  
  pure real(kind=dbl) function vnorm_fn(np, cajml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:)
    integer                       :: ij
    
    vnorm_fn = abs( cajml(1) )**2
    
    do ij = 1, np
      vnorm_fn = vnorm_fn +     sum( abs( cajml(jml(ij,0,-1):jml(ij, 0,+1)) )**2 ) + &
               &            2 * sum( abs( cajml(jml(ij,1,-1):jml(ij,ij,+1)) )**2 )
    end do
    
    vnorm_fn = sqrt( vnorm_fn )
    
  end function vnorm_fn
  
  pure real(kind=dbl) function tnorm_fn(np, cajml2)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(:)
    integer                       :: ij
    
    tnorm_fn = abs( cajml2(1) )**2 + sum( abs( cajml2(4:6) )**2 ) + 2 * sum( abs( cajml2(9:11) )**2 )
    
    do ij = 2, np
      tnorm_fn = tnorm_fn +     sum( abs( cajml2(jml2(ij,0,-2):jml2(ij, 0,+2)) )**2 ) + &
               &            2 * sum( abs( cajml2(jml2(ij,1,-2):jml2(ij,ij,+2)) )**2 )
    end do
    
    tnorm_fn = sqrt( tnorm_fn )
    
  end function tnorm_fn
  
  pure real(kind=dbl) function scalproduct_fn(np, cajm, cbjm)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(:), cbjm(:)
    integer                       :: ij
    
    scalproduct_fn = c2r_fn( cajm(1) * conjg( cbjm(1) ) )
    
    do ij = 1, np
      scalproduct_fn = scalproduct_fn +     c2r_fn(      cajm(jm(ij,0))           * conjg( cbjm(jm(ij,0)) )             ) + &
                     &                  2 * c2r_fn( sum( cajm(jm(ij,1):jm(ij,ij)) * conjg( cbjm(jm(ij,1):jm(ij,ij)) ) ) ) 
    end do
    
  end function scalproduct_fn
  
  pure real(kind=dbl) function dotproduct_fn(np, cajml, cbjml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:), cbjml(:)
    integer                       :: ij
    
    dotproduct_fn = c2r_fn( cajml(1) * conjg(cbjml(1)) )
    
    do ij = 1, np
      dotproduct_fn = dotproduct_fn +   c2r_fn(sum(cajml(jml(ij,0,-1):jml(ij, 0,+1))*conjg(cbjml(jml(ij,0,-1):jml(ij, 0,+1))))) + &
                                    2 * c2r_fn(sum(cajml(jml(ij,1,-1):jml(ij,ij,+1))*conjg(cbjml(jml(ij,1,-1):jml(ij,ij,+1)))))
    end do
    
  end function dotproduct_fn
  
end module sph_norms