module Spherical_func
  use Math
  implicit none
  
  public :: jm, jml, jml2, ersv_fn, ervs_fn, ervv_fn, ezvv_fn, ezvv_sub, snorm_fn, vnorm_fn, tnorm_fn, dotproduct_fn, scalproduct_fn
  
  contains
  
  pure integer function jm(ij, im)
    integer, intent(in) :: ij, im
    
    jm = ij*(ij+1)/2+im+1
    
  end function jm
  
  pure integer function jml(ij, im, il)
    integer, intent(in) :: ij, im, il
    
    jml = 3*(ij*(ij+1)/2+im)+il
    
  end function jml
  
  pure integer function jml2(ij, im, il)
    integer, intent(in) :: ij, im, il
    
    jml2 = 5*(ij*(ij+1)/2+im)+il-1
    
  end function jml2
  
  pure function ervs_fn(np, cajml) result(cjm)
    integer,           intent(in)  :: np
    complex(kind=dbl), intent(in)  :: cajml(:)
    complex(kind=dbl), allocatable :: cjm(:)
    integer                        :: ij, ijm
    real(kind=dbl)                 :: fac1, fac2
    
    allocate( cjm(jm(np,np)) ) ; cjm = czero
    
    ij = 0 ; cjm(1) = -cajml(1)
    
    do ij = 1, np
      fac1 = +sqrt( ( ij   ) / ( 2*ij + 1._dbl ) )
      fac2 = -sqrt( ( ij+1 ) / ( 2*ij + 1._dbl ) )
      
      do concurrent ( ijm = jm(ij,0) : jm(ij,ij) )
        cjm(ijm) = fac1 * cajml(3*ijm-4) + fac2 * cajml(3*ijm-2)
      end do
    end do
    
  end function ervs_fn
  
  pure function ersv_fn(np, cajm) result(cjml)
    integer,           intent(in)  :: np
    complex(kind=dbl), intent(in)  :: cajm(:)
    complex(kind=dbl), allocatable :: cjml(:)
    integer                        :: ij, ijm
    real(kind=dbl)                 :: fac1, fac2
    
    allocate( cjml(jml(np,np,+1)) ) ; cjml = czero
    
    ij = 0 ; cjml(1) = -cajm(1)
    
    do ij = 1, np
      fac1 = +sqrt( (ij  ) / ( 2*ij+1._dbl ) )
      fac2 = -sqrt( (ij+1) / ( 2*ij+1._dbl ) )
      
      do concurrent ( ijm = jm(ij,0) : jm(ij,ij) )
        cjml(3*ijm-4) = fac1 * cajm(ijm)
        cjml(3*ijm-2) = fac2 * cajm(ijm)
      end do
    end do
    
  end function ersv_fn
  
  pure function ervv_fn(np, cajml) result(cjml)
    integer,           intent(in)  :: np
    complex(kind=dbl), intent(in)  :: cajml(:)
    complex(kind=dbl), allocatable :: cjml(:)
    integer                        :: ij, ijm
    real(kind=dbl)                 :: fac1, fac2
    
    allocate( cjml(jml(np,np,+1)) ) ; cjml = czero
    
    ij = 0 ; cjml(1) = czero
    
    do ij = 1, np
      fac1 = sqrt( ( ij   ) / ( 2*ij+1._dbl ) )
      fac2 = sqrt( ( ij+1 ) / ( 2*ij+1._dbl ) )
      
      do concurrent ( ijm = jm(ij,0) : jm(ij,ij) )
        cjml(3*ijm-4) = cunit *   fac2 * cajml(3*ijm-3)
        cjml(3*ijm-3) = cunit * ( fac2 * cajml(3*ijm-4) + fac1 * cajml(3*ijm-2) )
        cjml(3*ijm-2) = cunit *   fac1 * cajml(3*ijm-3)
      end do
    end do
    
  end function ervv_fn
  
  pure function ezvv_fn(np, cajml) result(cjml)
    integer,           intent(in)  :: np
    complex(kind=dbl), intent(in)  :: cajml(:)
    complex(kind=dbl), allocatable :: cjml(:)
    integer                        :: ij, im, jm0, jm1, jm2
    
    allocate( cjml(jml(np,np,+1)) ) ; cjml = czero
    
    ij = 0 ; cjml(1) = sqrt(2._dbl/3._dbl) * cajml(3)
    
    ij = 1
      do im = 0, ij
        jm1 = 3 * ( (ij-1)*(ij  ) / 2 + im )
        jm0 = 3 * ( (ij  )*(ij+1) / 2 + im )
        jm2 = 3 * ( (ij+1)*(ij+2) / 2 + im )
        
        cjml(jm0-1) =                                                                      - im * cajml(jm0-1) / (ij       )
        cjml(jm0 )  = sqrt((ij+1)*((ij  )**2-im**2)/(2*ij+1._dbl)) * cajml(jm1+1) / (ij  ) - im * cajml(jm0  ) / (ij*(ij+1)) + &
                    & sqrt((ij  )*((ij+1)**2-im**2)/(2*ij+1._dbl)) * cajml(jm2-1) / (ij+1)
        cjml(jm0+1) =                                                                        im * cajml(jm0+1) / (    ij+1 ) + &
                    & sqrt((ij+2)*((ij+1)**2-im**2)/(2*ij+3._dbl)) * cajml(jm2  ) / (ij+1)
      end do
    
    do ij = 2, np-1
      do im = 0, ij
        jm1 = 3 * ( (ij-1)*(ij  ) / 2 + im )
        jm0 = 3 * ( (ij  )*(ij+1) / 2 + im )
        jm2 = 3 * ( (ij+1)*(ij+2) / 2 + im )
        
        cjml(jm0-1) = sqrt((ij-1)*((ij  )**2-im**2)/(2*ij-1._dbl)) * cajml(jm1  ) / (ij  ) - im * cajml(jm0-1) / (ij       )
        cjml(jm0  ) = sqrt((ij+1)*((ij  )**2-im**2)/(2*ij+1._dbl)) * cajml(jm1+1) / (ij  ) - im * cajml(jm0  ) / (ij*(ij+1)) + &
                    & sqrt((ij  )*((ij+1)**2-im**2)/(2*ij+1._dbl)) * cajml(jm2-1) / (ij+1)
        cjml(jm0+1) =                                                                        im * cajml(jm0+1) / (    ij+1 ) + &
                    & sqrt((ij+2)*((ij+1)**2-im**2)/(2*ij+3._dbl)) * cajml(jm2  ) / (ij+1)
      end do
    end do
    
    ij = np
      do im = 0, ij
        jm1 = 3 * ( (ij-1)*(ij  ) / 2 + im )
        jm0 = 3 * ( (ij  )*(ij+1) / 2 + im )
        jm2 = 3 * ( (ij+1)*(ij+2) / 2 + im )
        
        cjml(jm0-1) = sqrt(((ij-1)*(ij**2-im**2))/(2*ij-1._dbl)) * cajml(jm1  ) / ij - im * cajml(jm0-1) / (ij       )
        cjml(jm0  ) = sqrt(((ij+1)*(ij**2-im**2))/(2*ij+1._dbl)) * cajml(jm1+1) / ij - im * cajml(jm0  ) / (ij*(ij+1))
        cjml(jm0+1) =                                                                  im * cajml(jm0+1) / (    ij+1 )
      end do
    
    cjml = cunit * cjml
    
  end function ezvv_fn
  
  subroutine ezvv_sub(np, fac, cajml, cjml)
    integer,           intent(in)    :: np
    real(kind=dbl),    intent(in)    :: fac
    complex(kind=dbl), intent(in)    :: cajml(:)
    complex(kind=dbl), intent(inout) :: cjml(:,:)
    integer                          :: ij, im, ijm, jm0, jm1, jm2
    complex(kind=dbl)                :: cfac
    
    cfac = cunit * fac
    
    ij = 0
      im = 0
        ijm = 1
        
        cjml(1,ijm) = czero
        cjml(2,ijm) = czero
        cjml(3,ijm) = cjml(3,ijm) + sqrt(2._dbl/3._dbl) * cajml(3) * cfac
    
    ij = 1
      do im = 0, ij
        ijm = ijm+1
        
        jm1 = 3 * ( (ij-1)*(ij  ) / 2 + im )
        jm0 = 3 * ( (ij  )*(ij+1) / 2 + im )
        jm2 = 3 * ( (ij+1)*(ij+2) / 2 + im )
        
        cjml(1,ijm) = cjml(1,ijm) + (                                          -im * cajml(jm0-1) /  ij           ) * cfac
        cjml(2,ijm) = cjml(2,ijm) + ( sqrt((ij+1)*((ij  )**2-im**2)/(2*ij+1._dbl)) * cajml(jm1+1) /  ij         - &
                                    &                                           im * cajml(jm0  ) / (ij*(ij+1)) + &
                                    & sqrt((ij  )*((ij+1)**2-im**2)/(2*ij+1._dbl)) * cajml(jm2-1) /     (ij+1)    ) * cfac
        cjml(3,ijm) = cjml(3,ijm) + (                                           im * cajml(jm0+1) /     (ij+1)  + &
                                    & sqrt((ij+2)*((ij+1)**2-im**2)/(2*ij+3._dbl)) * cajml(jm2  ) /     (ij+1)    ) * cfac
      end do
    
    do ij = 2, np-1
      do im = 0, ij
        ijm = ijm+1
        
        jm1 = 3 * ( (ij-1)*(ij  ) / 2 + im )
        jm0 = 3 * ( (ij  )*(ij+1) / 2 + im )
        jm2 = 3 * ( (ij+1)*(ij+2) / 2 + im )
        
        cjml(1,ijm) = cjml(1,ijm) + ( sqrt((ij-1)*((ij  )**2-im**2)/(2*ij-1._dbl)) * cajml(jm1  ) /  ij         - &
                                    &                                           im * cajml(jm0-1) /  ij           ) * cfac
        cjml(2,ijm) = cjml(2,ijm) + ( sqrt((ij+1)*((ij  )**2-im**2)/(2*ij+1._dbl)) * cajml(jm1+1) /  ij         - &
                                    &                                           im * cajml(jm0  ) / (ij*(ij+1)) + &
                                    & sqrt((ij  )*((ij+1)**2-im**2)/(2*ij+1._dbl)) * cajml(jm2-1) /     (ij+1)    ) * cfac
        cjml(3,ijm) = cjml(3,ijm) + (                                           im * cajml(jm0+1) /     (ij+1)  + &
                                    & sqrt((ij+2)*((ij+1)**2-im**2)/(2*ij+3._dbl)) * cajml(jm2  ) /     (ij+1)    ) * cfac
      end do
    end do
    
    ij = np
      do im = 0, ij
        ijm = ijm+1
        
        jm1 = 3 * ( (ij-1)*(ij  ) / 2 + im )
        jm0 = 3 * ( (ij  )*(ij+1) / 2 + im )
        jm2 = 3 * ( (ij+1)*(ij+2) / 2 + im )
        
        cjml(1,ijm) = cjml(1,ijm) + ( sqrt((ij-1)*((ij  )**2-im**2)/(2*ij-1._dbl)) * cajml(jm1  ) /  ij       - &
                                    &                                           im * cajml(jm0-1) /  ij         ) * cfac
        cjml(2,ijm) = cjml(2,ijm) + ( sqrt((ij+1)*((ij  )**2-im**2)/(2*ij+1._dbl)) * cajml(jm1+1) /  ij       - &
                                    &                                           im * cajml(jm0  ) / (ij*(ij+1)) ) * cfac
        cjml(3,ijm) = cjml(3,ijm) + (                                           im * cajml(jm0+1) /     (ij+1)  ) * cfac
      end do
      
  end subroutine ezvv_sub
  
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
  
end module Spherical_func