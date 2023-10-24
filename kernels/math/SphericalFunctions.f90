module Spherical_func
  use Math
  implicit none
  
  public :: jm, jml, jml2
  public :: ersv_fn, ervs_fn, ervv_fn, ezvv_fn, ezvv_sub
  public :: snorm_fn, vnorm_fn, tnorm_fn, dotproduct_fn, scalproduct_fn
  
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
    
    ij = 0
      cjm(1) = -cajml(1)
    
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
    
    ij = 0
      cjml(1) = -cajm(1)
    
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
    
    ij = 0
      cjml(1) = czero
    
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
    integer,          intent(in) :: np
    complex(kind=dbl),intent(in) :: cajml(:)
    complex(kind=dbl)            :: cjml(jml(np,np,+1))
    integer                      :: j, m, jm0, jm1, jm2

    cjml(1) = sqrt(2._dbl/3._dbl)*cajml(3)

    j = 1
      do m = 0, j
        jm1 = 3*((j-1)*(j  )/2+m)
        jm0 = 3*((j  )*(j+1)/2+m)
        jm2 = 3*((j+1)*(j+2)/2+m)

        cjml(jm0-1) =                                                                                - m*cajml(jm0-1)/(j      )
        cjml(jm0 )  = sqrt(((j+1._dbl)*((j       )*(j       )-m*m))/(2*j+1._dbl))*cajml(jm1+1)/(j  ) - m*cajml(jm0  )/(j*(j+1)) + &
                    & sqrt(((j       )*((j+1._dbl)*(j+1._dbl)-m*m))/(2*j+1._dbl))*cajml(jm2-1)/(j+1)
        cjml(jm0+1) =                                                                                  m*cajml(jm0+1)/(   j+1 ) + &
                    & sqrt(((j+2._dbl)*((j+1._dbl)*(j+1._dbl)-m*m))/(2*j+3._dbl))*cajml(jm2  )/(j+1)
      end do

    do j = 2, np-1
      do m = 0, j
        jm1 = 3*((j-1)*(j  )/2+m)
        jm0 = 3*((j  )*(j+1)/2+m)
        jm2 = 3*((j+1)*(j+2)/2+m)

        cjml(jm0-1) = sqrt(((j-1._dbl)*((j       )*(j       )-m*m))/(2*j-1._dbl))*cajml(jm1  )/(j  ) - m*cajml(jm0-1)/(j      )
        cjml(jm0  ) = sqrt(((j+1._dbl)*((j       )*(j       )-m*m))/(2*j+1._dbl))*cajml(jm1+1)/(j  ) - m*cajml(jm0  )/(j*(j+1)) + &
                    & sqrt(((j       )*((j+1._dbl)*(j+1._dbl)-m*m))/(2*j+1._dbl))*cajml(jm2-1)/(j+1)
        cjml(jm0+1) =                                                                                  m*cajml(jm0+1)/(   j+1 ) + &
                    & sqrt(((j+2._dbl)*((j+1._dbl)*(j+1._dbl)-m*m))/(2*j+3._dbl))*cajml(jm2  )/(j+1)
      end do
    end do

    j = np
      do m = 0, j
        jm1 = 3*((j-1)*(j  )/2+m)
        jm0 = 3*((j  )*(j+1)/2+m)
        jm2 = 3*((j+1)*(j+2)/2+m)

        cjml(jm0-1) = sqrt(((j-1._dbl)*((j  )*(j  )-m*m))/(2*j-1._dbl))*cajml(jm1  )/(j  ) - m*cajml(jm0-1)/(j      )
        cjml(jm0  ) = sqrt(((j+1._dbl)*((j  )*(j  )-m*m))/(2*j+1._dbl))*cajml(jm1+1)/(j  ) - m*cajml(jm0  )/(j*(j+1))
        cjml(jm0+1) =                                                                        m*cajml(jm0+1)/(   j+1 )
      end do

    cjml = cunit * cjml

  end function ezvv_fn
  
  subroutine ezvv_sub(np, fac, cajml, cjml)
    integer,           intent(in)    :: np
    real(kind=dbl),    intent(in)    :: fac
    complex(kind=dbl), intent(in)    :: cajml(:)
    complex(kind=dbl), intent(inout) :: cjml(:,:)
    integer                          :: j, m, ijm, jm0, jm1, jm2
    complex(kind=dbl)                :: cfac

    cfac = cunit * fac

    j = 0
      m = 0
        ijm = 1

        cjml(1,ijm) = czero
        cjml(2,ijm) = czero
        cjml(3,ijm) = cjml(3,ijm) + sqrt(2._dbl/3._dbl) * cajml(3) * cfac

    j = 1
      do m = 0, j
        ijm = ijm+1

        jm1 = 3*((j-1)*(j  )/2+m)
        jm0 = 3*((j  )*(j+1)/2+m)
        jm2 = 3*((j+1)*(j+2)/2+m)

        cjml(1,ijm) = cjml(1,ijm) + (                                           -m * cajml(jm0-1) /  j          ) * cfac
        cjml(2,ijm) = cjml(2,ijm) + ( sqrt(((j+1)*((j  )*(j  )-m*m))/(2*j+1._dbl)) * cajml(jm1+1) /  j -        &
                                    &                                            m * cajml(jm0  ) / (j*(j+1)) + &
                                    & sqrt(((j  )*((j+1)*(j+1)-m*m))/(2*j+1._dbl)) * cajml(jm2-1) /    (j+1)    ) * cfac
        cjml(3,ijm) = cjml(3,ijm) + (                                            m * cajml(jm0+1) /    (j+1)  + &
                                    & sqrt(((j+2)*((j+1)*(j+1)-m*m))/(2*j+3._dbl)) * cajml(jm2  ) /    (j+1)    ) * cfac
      end do

    do j = 2, np-1
      do m = 0, j
        ijm = ijm+1

        jm1 = 3*((j-1)*(j  )/2+m)
        jm0 = 3*((j  )*(j+1)/2+m)
        jm2 = 3*((j+1)*(j+2)/2+m)

        cjml(1,ijm) = cjml(1,ijm) + ( sqrt(((j-1)*((j  )*(j  )-m*m))/(2*j-1._dbl)) * cajml(jm1  ) /  j -        &
                                    &                                            m * cajml(jm0-1) /  j          ) * cfac
        cjml(2,ijm) = cjml(2,ijm) + ( sqrt(((j+1)*((j  )*(j  )-m*m))/(2*j+1._dbl)) * cajml(jm1+1) /  j -        &
                                    &                                            m * cajml(jm0  ) / (j*(j+1)) + &
                                    & sqrt(((j  )*((j+1)*(j+1)-m*m))/(2*j+1._dbl)) * cajml(jm2-1) /    (j+1)    ) * cfac
        cjml(3,ijm) = cjml(3,ijm) + (                                            m * cajml(jm0+1) /    (j+1) +  &
                                    & sqrt(((j+2)*((j+1)*(j+1)-m*m))/(2*j+3._dbl)) * cajml(jm2  ) /    (j+1)    ) * cfac
      end do
    end do

    j = np
      do m = 0, j
        ijm = ijm+1

        jm1 = 3*((j-1)*(j  )/2+m)
        jm0 = 3*((j  )*(j+1)/2+m)
        jm2 = 3*((j+1)*(j+2)/2+m)

        cjml(1,ijm) = cjml(1,ijm) + ( sqrt(((j-1)*((j  )*(j  )-m*m))/(2*j-1._dbl)) * cajml(jm1  ) /  j -      &
                                    &                                            m * cajml(jm0-1) /  j        ) * cfac
        cjml(2,ijm) = cjml(2,ijm) + ( sqrt(((j+1)*((j  )*(j  )-m*m))/(2*j+1._dbl)) * cajml(jm1+1) /  j -      &
                                    &                                            m * cajml(jm0  ) / (j*(j+1)) ) * cfac
        cjml(3,ijm) = cjml(3,ijm) + (                                            m * cajml(jm0+1) /    (j+1)  ) * cfac
      end do

  end subroutine ezvv_sub
  
  pure real(kind=dbl) function snorm_fn(np, cajm)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(:)
    integer                       :: j, m

    snorm_fn = 0._dbl

    do j = 0, np
      m = 0
        snorm_fn = snorm_fn + abs(cajm(jm(j,m)))**2

      do m = 1, j
        snorm_fn = snorm_fn + 2*abs(cajm(jm(j,m)))**2
      end do
    end do

    snorm_fn = sqrt(snorm_fn)

  end function snorm_fn
  
  pure real(kind=dbl) function vnorm_fn(np, cajml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:)
    integer                       :: j, m, l

    vnorm_fn = 0._dbl

    do j = 0, np
      m = 0
        do l = abs(j-1)-j, +1
          vnorm_fn = vnorm_fn + abs(cajml(jml(j,m,l)))**2
        end do

      do m = 1, j
        do l = abs(j-1)-j, +1
          vnorm_fn = vnorm_fn + 2*abs(cajml(jml(j,m,l)))**2
        end do
      end do
    end do

    vnorm_fn = sqrt(vnorm_fn)

  end function vnorm_fn
  
  pure real(kind=dbl) function tnorm_fn(np, cajml2)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml2(:)
    integer                       :: j, m, l

    tnorm_fn = 0._dbl

    do j = 0, np
      m = 0
        do l = abs(j-2)-j, +2
          tnorm_fn = tnorm_fn + abs(cajml2(jml2(j,m,l)))**2
        end do

      do m = 1, j
        do l = abs(j-2)-j, +2
          tnorm_fn = tnorm_fn + 2*abs(cajml2(jml2(j,m,l)))**2
        end do
      end do
    end do

    tnorm_fn = sqrt(tnorm_fn)

  end function tnorm_fn
  
  pure real(kind=dbl) function scalproduct_fn(np, cajm, cbjm)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(:), cbjm(:)
    integer                       :: j, m
    
    scalproduct_fn = 0._dbl
    
    do j = 0, np
      m = 0
        scalproduct_fn = scalproduct_fn + c2r_fn( cajm(jm(j,m)) * conjg(cbjm(jm(j,m))) )
      
      do m = 1, j
        scalproduct_fn = scalproduct_fn +  2  * c2r_fn( cajm(jm(j,m)) * conjg(cbjm(jm(j,m))) )
      end do
    end do
    
  end function scalproduct_fn
  
  pure real(kind=dbl) function dotproduct_fn(np, cajml, cbjml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:), cbjml(:)
    integer                       :: j, m, l
    
    dotproduct_fn = 0._dbl
    
    do j = 0, np
      m = 0
        do l = abs(j-1)-j, +1
          dotproduct_fn = dotproduct_fn + c2r_fn( cajml(jml(j,m,l)) * conjg(cbjml(jml(j,m,l))) )
        end do
      
      do m = 1, j
        do l = abs(j-1)-j, +1
          dotproduct_fn = dotproduct_fn + 2 * c2r_fn( cajml(jml(j,m,l)) * conjg(cbjml(jml(j,m,l))) )
        end do
      end do
    end do
    
  end function dotproduct_fn
  
end module Spherical_func