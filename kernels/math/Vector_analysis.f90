module Vector_analysis
  use Math
  use Clebsch_Legendre
  use Spherical_func
  implicit none
  
  public :: vec2scals_sub, vecxyz2rtp_sub
  
  contains
  
  module pure subroutine vec2scals_sub(jmax_in, vec, x, y, z)
    integer,           intent(in)  :: jmax_in
    complex(kind=dbl), intent(in)  :: vec(:)
    complex(kind=dbl), intent(out) :: x(:), y(:), z(:)
    integer                        :: ij, im, il, ijm, ijml, zn
    complex(kind=dbl)              :: sum1, sum2, sum3
    
    x = czero ; y = czero ; z = czero
    
    im = 0
      do ij = im, jmax_in+1
        sum1 = czero ; sum2 = czero ; sum3 = czero
        
        zn = +1
        do il = abs(ij-1), min(jmax_in, ij+1)
          zn = -zn ; ijml = jml(il,im,ij-il)
          
          sum1 = sum1 + zn * conjg( vec(ijml+3) ) * cleb1_fn(ij,im,1,-1,il,im-1)
          sum3 = sum3 +             vec(ijml  )   * cleb1_fn(ij,im,1, 0,il,im  )
          sum2 = sum2 +             vec(ijml+3)   * cleb1_fn(ij,im,1,+1,il,im+1)
        end do
        
        ijm = ij*(ij+1)/2+im+1
          x(ijm) =         ( sum1 - sum2 ) / sqrt(2._dbl)
          y(ijm) = cunit * ( sum1 + sum2 ) / sqrt(2._dbl)
          z(ijm) =           sum3 
      end do
    
    do im = 1, jmax_in+1
      do ij = im, jmax_in+1
        sum1 = czero ; sum2 = czero ; sum3 = czero
        
        do il = abs(ij-1), min(jmax_in, ij+1)
          ijml = jml(il,im,ij-il)
          
          if ( il > im ) then
              sum1 = sum1 + vec(ijml  ) * cleb1_fn(ij,im,1,-1,il,im-1)
              sum3 = sum3 + vec(ijml+3) * cleb1_fn(ij,im,1, 0,il,im  )
              sum2 = sum2 + vec(ijml+6) * cleb1_fn(ij,im,1,+1,il,im+1)
          else if ( il > im-1 ) then
              sum1 = sum1 + vec(ijml  ) * cleb1_fn(ij,im,1,-1,il,im-1)
              sum3 = sum3 + vec(ijml+3) * cleb1_fn(ij,im,1, 0,il,im  )
          else
              sum1 = sum1 + vec(ijml  ) * cleb1_fn(ij,im,1,-1,il,im-1)
          end if
        end do
        
        ijm = ij*(ij+1)/2+im+1
          x(ijm) =         ( sum1 - sum2 ) / sqrt(2._dbl)
          y(ijm) = cunit * ( sum1 + sum2 ) / sqrt(2._dbl)
          z(ijm) =           sum3 
      end do
    end do
    
  end subroutine vec2scals_sub
  
  module pure subroutine vecxyz2rtp_sub(theta, phi, vx, vy, vz, vr, vth, vph)
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
  
end module Vector_analysis