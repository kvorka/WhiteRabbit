module Spherical_func
  use Math
  implicit none

  public :: ersv_fn     !er*skalar1          -> vektor2
  public :: ervs_fn     !er*vektor1          -> skalar2
  public :: ervv_fn     !er*vektor1          -> vektor2
  public :: ezvv_fn     !ez*vektor1          -> vektor2

  public :: ezvv_sub
  
  public :: snorm_fn       !L2 norma skalarneho  radu
  public :: vnorm_fn       !L2 norma vektoroveho radu
  public :: tnorm_fn       !L2 norma tenzoroveho radu (k == 2)
  public :: dotproduct_fn  !L2 norma skalarneho sucinu vektorov
  public :: scalproduct_fn !L2 norma sucinu skalarov
  
  public :: jm    !(j,m)     -> jm   = j*(j+1)/2+m+1
  public :: jml   !(j,m,l-j) -> jml  = 3*(j*(j+1)/2+m)+l-j
  public :: jml2  !(j,m,l-j) -> jml2 = 5*(j*(j+1)/2+m)+l-(j+2)+1

  contains

  pure integer function jm(j,m)
    integer, intent(in) :: j, m

    jm = j*(j+1)/2+m+1

  end function jm

  pure integer function jml(j,m,l)
    integer, intent(in) :: j, m, l

    jml = 3*(j*(j+1)/2+m)+l

  end function jml

  pure integer function jml2(j, m, l)
    integer, intent(in) :: j, m, l

    jml2 = 5*(j*(j+1)/2+m)+l-1

  end function jml2

  pure function ervs_fn(np, cajml) result(cjm)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:)
    complex(kind=dbl)             :: cjm(jm(np,np))
    integer                       :: j, m

    cjm(1) = -cajml(1)

    do j = 1, np
      do m = 0, j
        cjm(jm(j,m)) = sqrt(j/(2*j+1._dbl)) * cajml(jml(j,m,-1)) - sqrt((j+1._dbl)/(2*j+1._dbl)) * cajml(jml(j,m,+1))
      end do
    end do

  end function ervs_fn

  pure function ersv_fn(np, cajm) result(cjml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajm(:)
    complex(kind=dbl)             :: cjml(jml(np,np,+1))
    integer                       :: j, m

    cjml(1) = -cajm(1)

    do j = 1, np
      do m = 0, j
        cjml(jml(j,m,-1)) = +sqrt((j       )/(2*j+1._dbl))*cajm(jm(j,m))
        cjml(jml(j,m, 0)) = cmplx(0._dbl, 0._dbl, kind=dbl)
        cjml(jml(j,m,+1)) = -sqrt((j+1._dbl)/(2*j+1._dbl))*cajm(jm(j,m))
      end do
    end do

  end function ersv_fn

  pure function ervv_fn(np, cajml) result(cjml)
    integer,           intent(in) :: np
    complex(kind=dbl), intent(in) :: cajml(:)
    complex(kind=dbl)             :: cjml(jml(np,np,+1))
    integer                       :: j, m

    cjml(1) = cmplx(0._dbl, 0._dbl, kind=dbl)

    do j = 1, np
      do m = 0, j
        cjml(jml(j,m,-1)) = sqrt((j+1._dbl)/(2*j+1._dbl))*cajml(jml(j,m, 0))
        cjml(jml(j,m, 0)) = sqrt((j+1._dbl)/(2*j+1._dbl))*cajml(jml(j,m,-1)) + sqrt(j/(2*j+1._dbl))*cajml(jml(j,m,+1))
        cjml(jml(j,m,+1)) =                                                    sqrt(j/(2*j+1._dbl))*cajml(jml(j,m, 0))
      end do
    end do

    cjml = cunit*cjml

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

    cjml = cunit*cjml

  end function ezvv_fn

  subroutine ezvv_sub(np, cajml, cjml)
    integer,           intent(in)  :: np
    complex(kind=dbl), intent(in)  :: cajml(:)
    complex(kind=dbl), intent(out) :: cjml(:)
    integer                        :: j, m, jm0, jm1, jm2

    cjml(1) = sqrt(2._dbl/3._dbl) * cajml(3)

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
        scalproduct_fn = scalproduct_fn + real(cajm(jm(j,m))*conjg(cbjm(jm(j,m))), kind=dbl)

      do m = 1, j
        scalproduct_fn = scalproduct_fn +  2  * real(cajm(jm(j,m))*conjg(cbjm(jm(j,m))), kind=dbl)
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
          dotproduct_fn = dotproduct_fn + real(cajml(jml(j,m,l))*conjg(cbjml(jml(j,m,l))), kind=dbl)
        end do

      do m = 1, j
        do l = abs(j-1)-j, +1
          dotproduct_fn = dotproduct_fn + 2*real(cajml(jml(j,m,l))*conjg(cbjml(jml(j,m,l))), kind=dbl)
        end do
      end do
    end do

  end function dotproduct_fn

end module Spherical_func