module sph_unitvec_op
  use Math
  use sph_indexing
  implicit none; public; contains
  
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
  
  pure subroutine ezvv_sub(np, fac, cajml, cjml)
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
  
end module sph_unitvec_op