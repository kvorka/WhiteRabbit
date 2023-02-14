module Vector_analysis
  use Math
  implicit none
  
  public :: vec2scals_sub
  public :: vecxyz2rtp_sub
  private :: cleb1_fn
  
  contains
  
  subroutine vec2scals_sub(jmax_in, vec, x, y, z)
    integer,           intent(in)  :: jmax_in
    complex(kind=dbl), intent(in)  :: vec(:)
    complex(kind=dbl), intent(out) :: x(:), y(:), z(:)
    integer                        :: j, m, l, jm, s
    complex(kind=dbl)              :: sum1, sum2, sum3
  
    x = cmplx(0._dbl, 0._dbl, kind=dbl)
    y = cmplx(0._dbl, 0._dbl, kind=dbl)
    z = cmplx(0._dbl, 0._dbl, kind=dbl)
  
    m = 0
      do j = m, jmax_in+1
        sum1 = cmplx(0._dbl, 0._dbl, kind=dbl)
        sum2 = cmplx(0._dbl, 0._dbl, kind=dbl)
        sum3 = cmplx(0._dbl, 0._dbl, kind=dbl)
  
        s = +1
        do l = abs(j-1), min(jmax_in, j+1)
          s = -s
          sum1 = sum1 + s*conjg( vec(3*(l*(l+1)/2+m+1)+j-l) ) * cleb1_fn(j,m,1,-1,l,m-1)
          sum3 = sum3 +          vec(3*(l*(l+1)/2+m  )+j-l)   * cleb1_fn(j,m,1, 0,l,m  )
          sum2 = sum2 +          vec(3*(l*(l+1)/2+m+1)+j-l)   * cleb1_fn(j,m,1,+1,l,m+1)
        end do
  
        jm = j*(j+1)/2+m+1
          x(jm) =         (sum1-sum2) / sqrt(2._dbl)
          y(jm) = cunit * (sum1+sum2) / sqrt(2._dbl)
          z(jm) =          sum3
      end do
  
    do m = 1, jmax_in+1
      do j = m, jmax_in+1
        sum1 = cmplx(0._dbl, 0._dbl, kind=dbl)
        sum2 = cmplx(0._dbl, 0._dbl, kind=dbl)
        sum3 = cmplx(0._dbl, 0._dbl, kind=dbl)
  
        do l = abs(j-1), min(jmax_in, j+1)
                       sum1 = sum1 + vec(3*(l*(l+1)/2+m-1)+j-l) * cleb1_fn(j,m,1,-1,l,m-1)
          if (l > m-1) sum3 = sum3 + vec(3*(l*(l+1)/2+m  )+j-l) * cleb1_fn(j,m,1, 0,l,m  )
          if (l > m+0) sum2 = sum2 + vec(3*(l*(l+1)/2+m+1)+j-l) * cleb1_fn(j,m,1,+1,l,m+1)
        end do
  
        jm = j*(j+1)/2+m+1
          x(jm) =         (sum1-sum2) / sqrt(2._dbl)
          y(jm) = cunit * (sum1+sum2) / sqrt(2._dbl)
          z(jm) =          sum3
      end do
    end do
  
  end subroutine vec2scals_sub

  subroutine vecxyz2rtp_sub(theta, phi, vx, vy, vz, vr, vth, vph)
    real(kind=dbl), intent(in)  :: theta, phi, vx, vy, vz
    real(kind=dbl), intent(out) :: vr, vth, vph
    real(kind=dbl)              :: cph, sph, ct, st

    cph = cos(phi   * pi / 180)
    sph = sin(phi   * pi / 180)
    ct  = cos(theta * pi / 180)
    st  = sin(theta * pi / 180)
    
    vr  =  vx * cph * st + vy * sph * st + vz * ct
    vth =  vx * cph * ct + vy * sph * ct - vz * st
    vph = -vx * sph      + vy * cph
    
  end subroutine vecxyz2rtp_sub

    function cleb1_fn(j1, m1, j2, m2, j, m) result(cleb1)
      integer,       intent(in) :: j1, m1, j2, m2, j, m
      real(kind=dbl)            :: cleb1
    
      if ((j2 /= 1) .or. (abs(j1-j) > 1) .or. ((j1+j) == 0) .or. (abs(m2) > 1) .or. abs(m1) > j1) then
        cleb1 = 0._dbl

      else if (m2 == -1) then
        if (j1 == j-1) cleb1 = +sqrt((j-m-1._dbl)*(j-m       )/((2*j-1._dbl)*(  j       ))/2)
        if (j1 == j  ) cleb1 = +sqrt((j+m+1._dbl)*(j-m       )/((  j+1._dbl)*(  j       ))/2)
        if (j1 == j+1) cleb1 = +sqrt((j+m+2._dbl)*(j+m+1._dbl)/((  j+1._dbl)*(2*j+3._dbl))/2)

      else if (m2 == 0) then
        if (j1 == j-1) cleb1 = +sqrt((j+m       )*(j-m       )/ (2*j-1._dbl)/(  j       ))
        if (j1 == j  ) cleb1 = +sqrt((  m       )*(  m       )/ (  j+1._dbl)/(  j       ))
        if (j1 == j+1) cleb1 = -sqrt((j+m+1._dbl)*(j-m+1._dbl)/((  j+1._dbl)*(2*j+3._dbl)))
    
      else
        if (j1 == j-1) cleb1 = +sqrt((j+m-1._dbl)*(j+m       )/((2*j-1._dbl)*(  j       ))/2)
        if (j1 == j  ) cleb1 = -sqrt((j+m       )*(j-m+1._dbl)/((  j+1._dbl)*(  j       ))/2)
        if (j1 == j+1) cleb1 = +sqrt((j-m+1._dbl)*(j-m+2._dbl)/((  j+1._dbl)*(2*j+3._dbl))/2)
      end if
    
    end function cleb1_fn
     
  end module Vector_analysis