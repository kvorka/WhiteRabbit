module FFT
  use Math
  implicit none
  
  integer, private :: npoints, psym, factor(21), sym(21), unsym(21)

  public  :: init_fft_sub, fft_sub, fft_inv_sub
  private :: ecwc06, ecyc06, ecqc06, ecrc06, ecsc06, ectc06, ecuc06, ecvc06, ecxc06

  contains

  subroutine init_fft_sub(Nfour)
    integer, intent(in) :: Nfour
    integer             :: f, j, n, p, ptwo, q, r, pmax, twogrp, pp(10), qq(20)
    
    npoints = Nfour
    pmax    = 19
    twogrp  = 8

    n = npoints; psym = 1; f = 2; p = 0; q = 0
    
    do
      if (n <= 1) exit
      
      do j = f, pmax
         if (n == (n/j)*j) exit
      end do

      f = j
      n = n/f

      if (n /= (n/f)*f) then
        q = q + 1
        qq(q) = f
        
      else
        n     = n/f
        p     = p + 1
        pp(p) = f
        psym  = psym*f
        
      end if
    end do

    r = 1; if (q == 0) r = 0

    if (p >= 1) then
      sym(1:p)    = pp(p:1:-1)
      factor(1:p) = pp(p:1:-1)
      
      factor(p+q+1:2*p+q) = pp(1:p)
      sym(p+r+1:2*p+r)    = pp(1:p)
    end if
    
    if (q >= 1) then
      unsym(1:q)      = qq(1:q)
      factor(p+1:p+q) = qq(1:q)
      
      sym(p+1) = npoints/psym**2
    end if

    factor(2*p+q+1) = 0
    ptwo            = 1

    j = 0

    do
      j = j + 1; if (factor(j) == 0) exit
      
      if (factor(j) == 2) then
        ptwo      = 2 * ptwo
        factor(j) = 1
      
        if ((ptwo >= twogrp) .or. (factor(j+1) /= 2)) then
          factor(j) = ptwo
          ptwo      = 1
          
        end if
      end if
    end do
    
    if (p == 0) r = 0; sym(2*p+r+1) = 0
    if (q <= 1) q = 0; unsym(q+1)   = 0
    
  end subroutine init_fft_sub

  subroutine fft_inv_sub(x, y)
    real(kind=dbl), intent(inout) :: x(:), y(:)
    
    y = -y
    call fft_sub(x,y)
    y = -y
    
  end subroutine fft_inv_sub

  subroutine fft_sub(x, y)
    real(kind=dbl), intent(inout) :: x(:), y(:)

    call ecwc06(x, y)
    call ecyc06(x, y)

    x = x / sqrt(real(npoints, kind=dbl))
    y = y / sqrt(real(npoints, kind=dbl))

  end subroutine fft_sub

    subroutine ecwc06(x, y)
      real(kind=dbl), intent(inout) :: x(:), y(:)
      integer                       :: f, p, m
  
      f = 0; m = npoints
  
      do
        f = f + 1; p = factor(f)
        
        if (p == 0) then 
          exit
          
        else if (p /= 1) then
          m = m/p
          
          select case (p)
            case(2)     ; call ecvc06(x, y, m)
            case(3)     ; call ecuc06(x, y, m)
            case(4)     ; call ectc06(x, y, m)
            case(5)     ; call ecsc06(x, y, m)
            case(8)     ; call ecrc06(x, y, m)
            case default; call ecqc06(x, y, m, p)
          end select
        end if
      end do
      
    end subroutine ecwc06
    
      subroutine ecqc06(x, y, m, p)
        integer,        intent(in)    :: m, p
        real(kind=dbl), intent(inout) :: x(:), y(:)
        integer                       :: j, k, u, v
        real(kind=dbl)                :: is, iu, rs, ru, t, xt, yt, angle
        real(kind=dbl)                :: a(18), b(18), c(18), s(18), ia(9), ib(9), ra(9), rb(9), aa(9,9), bb(9,9)
      
        do u = 1, p/2
          angle = 2 * pi * u / p
          
          a(u)   = cos(angle)
          b(u)   = sin(angle)
          
          a(p-u) = a(u)
          b(p-u) = -b(u)
        end do
      
        do u = 1, p/2
          do v = 1, p/2
            aa(v,u) = a(u*v - ((u*v)/p)*p)
            bb(v,u) = b(u*v - ((u*v)/p)*p)
          end do
        end do
      
        do j = 0, m/2
          angle = 2 * pi * j / m / p
          
          c(1) = cos(angle)
          s(1) = sin(angle)
            do u = 2, p-1
              c(u) = c(u-1)*c(1) - s(u-1)*s(1)
              s(u) = s(u-1)*c(1) + c(u-1)*s(1)
            end do
      
          do k = j+1, npoints, m*p
            xt  = x(k)
            yt  = y(k)
            rs  = x(m+k) + x((p-1)*m+k)
            is  = y(m+k) + y((p-1)*m+k)
            ru  = x(m+k) - x((p-1)*m+k)
            iu  = y(m+k) - y((p-1)*m+k)
      
            do u = 1, p/2
              ra(u) = xt + rs*aa(u,1)
              ia(u) = yt + is*aa(u,1)
              rb(u) = ru*bb(u,1)
              ib(u) = iu*bb(u,1)
            end do
      
            xt = xt + rs
            yt = yt + is
      
            do u = 2, p/2
              rs  = x(u*m+k) + x((p-u)*m+k)
              is  = y(u*m+k) + y((p-u)*m+k)
              ru  = x(u*m+k) - x((p-u)*m+k)
              iu  = y(u*m+k) - y((p-u)*m+k)
              
              xt  = xt + rs
              yt  = yt + is
      
              do v = 1, p/2
                ra(v) = ra(v) + rs*aa(v,u)
                ia(v) = ia(v) + is*aa(v,u)
                rb(v) = rb(v) + ru*bb(v,u)
                ib(v) = ib(v) + iu*bb(v,u)
              end do
            end do
      
            x(k) = xt
            y(k) = yt
      
            do u = 1, p/2
              xt = ra(u) + ib(u)
              yt = ia(u) - rb(u)
                x(u*m+k) = xt*c(u) + yt*s(u)
                y(u*m+k) = yt*c(u) - xt*s(u)
                
              xt = ra(u) - ib(u)
              yt = ia(u) + rb(u)
                x((p-u)*m+k) = xt*c(p-u) + yt*s(p-u)
                y((p-u)*m+k) = yt*c(p-u) - xt*s(p-u)
            end do
          end do
      
          if ((j > 0) .and. (2*j < m)) then
            do u = 1, p-1
              t    = +c(u)*a(u) + s(u)*b(u)
              s(u) = -s(u)*a(u) + c(u)*b(u)
              c(u) = t
            end do
            
            do k = m+1-j, npoints, m*p
              xt  = x(k)
              yt  = y(k)
              rs  = x(m+k) + x((p-1)*m+k)
              is  = y(m+k) + y((p-1)*m+k)
              ru  = x(m+k) - x((p-1)*m+k)
              iu  = y(m+k) - y((p-1)*m+k)
              
              do u = 1, p/2
                ra(u) = xt + rs*aa(u,1)
                ia(u) = yt + is*aa(u,1)
                rb(u) = ru*bb(u,1)
                ib(u) = iu*bb(u,1)
              end do
              
              xt = xt + rs
              yt = yt + is
              
              do u = 2, p/2
                rs  = x(u*m+k) + x((p-u)*m+k)
                is  = y(u*m+k) + y((p-u)*m+k)
                ru  = x(u*m+k) - x((p-u)*m+k)
                iu  = y(u*m+k) - y((p-u)*m+k)
                
                xt  = xt + rs
                yt  = yt + is
              
                do v = 1, p/2
                  ra(v) = ra(v) + rs*aa(v,u)
                  ia(v) = ia(v) + is*aa(v,u)
                  rb(v) = rb(v) + ru*bb(v,u)
                  ib(v) = ib(v) + iu*bb(v,u)
                end do
              end do
              
              x(k) = xt
              y(k) = yt
              
              do u = 1, p/2
                xt = ra(u) + ib(u)
                yt = ia(u) - rb(u)
                  x(u*m+k) = xt*c(u) + yt*s(u)
                  y(u*m+k) = yt*c(u) - xt*s(u)
                  
                xt = ra(u) - ib(u)
                yt = ia(u) + rb(u)
                  x((p-u)*m+k) = xt*c(p-u) + yt*s(p-u)
                  y((p-u)*m+k) = yt*c(p-u) - xt*s(p-u)
              end do
            end do
          end if
        end do
      
      end subroutine ecqc06
      
      subroutine ecrc06(x, y, m)
        integer,        intent(in)    :: m
        real(kind=dbl), intent(inout) :: x(:), y(:)
        integer                       :: j, k
        real(kind=dbl)                :: e, t, c1, c2, c3, c4, c5, c6, c7, s1, s2, s3, s4, s5, s6, s7, angle,                 &
        & i1, i2, i3, i4, i5, i6, i7, is0, is1, is2, is3, iss0, iss1, isu0, isu1, iu0, iu1, iu2, iu3, ius0, ius1, iuu0, iuu1, &
        & r1, r2, r3, r4, r5, r6, r7, rs0, rs1, rs2, rs3, rss0, rss1, rsu0, rsu1, ru0, ru1, ru2, ru3, rus0, rus1, ruu0, ruu1
      
        e = cos(pi/4)
      
        do j = 0, m/2
          angle = pi * j / m / 4
            c1 = cos(angle)
            s1 = sin(angle)
              c2 = c1*c1 - s1*s1
              s2 = s1*c1 + c1*s1
              c3 = c2*c1 - s2*s1
              s3 = s2*c1 + c2*s1
              c4 = c2*c2 - s2*s2
              s4 = s2*c2 + c2*s2
              c5 = c4*c1 - s4*s1
              s5 = s4*c1 + c4*s1
              c6 = c4*c2 - s4*s2
              s6 = s4*c2 + c4*s2
              c7 = c4*c3 - s4*s3
              s7 = s4*c3 + c4*s3
          
          do k = j+1, npoints, 8*m
            rs0 = x(k    ) + x(k+4*m)
            is0 = y(k    ) + y(k+4*m)
            ru0 = x(k    ) - x(k+4*m)
            iu0 = y(k    ) - y(k+4*m)
            rs1 = x(k+  m) + x(k+5*m)
            is1 = y(k+  m) + y(k+5*m)
            ru1 = x(k+  m) - x(k+5*m)
            iu1 = y(k+  m) - y(k+5*m)
            rs2 = x(k+2*m) + x(k+6*m)
            is2 = y(k+2*m) + y(k+6*m)
            ru2 = x(k+2*m) - x(k+6*m)
            iu2 = y(k+2*m) - y(k+6*m)
            rs3 = x(k+3*m) + x(k+7*m)
            is3 = y(k+3*m) + y(k+7*m)
            ru3 = x(k+3*m) - x(k+7*m)
            iu3 = y(k+3*m) - y(k+7*m)
            
            rss0 = rs0 + rs2
            iss0 = is0 + is2
            rsu0 = rs0 - rs2
            isu0 = is0 - is2
            rss1 = rs1 + rs3
            iss1 = is1 + is3
            rsu1 = rs1 - rs3
            isu1 = is1 - is3
            rus0 = ru0 - iu2
            ius0 = iu0 + ru2
            ruu0 = ru0 + iu2
            iuu0 = iu0 - ru2
            rus1 = ru1 - iu3
            ius1 = iu1 + ru3
            ruu1 = ru1 + iu3
            iuu1 = iu1 - ru3
            t = (rus1+ius1)*e
            ius1 = (ius1-rus1)*e
            rus1 = t
            t = (ruu1+iuu1)*e
            iuu1 = (iuu1-ruu1)*e
            ruu1 = t
            
            r1 = ruu0 + ruu1
            i1 = iuu0 + iuu1
            r2 = rsu0 + isu1
            i2 = isu0 - rsu1
            r3 = rus0 + ius1
            i3 = ius0 - rus1
            r4 = rss0 - rss1
            i4 = iss0 - iss1
            r5 = ruu0 - ruu1
            i5 = iuu0 - iuu1
            r6 = rsu0 - isu1
            i6 = isu0 + rsu1
            r7 = rus0 - ius1
            i7 = ius0 + rus1
            
            x(k    ) = rss0 + rss1
            y(k    ) = iss0 + iss1
            x(k+  m) = r4*c4 + i4*s4
            y(k+  m) = i4*c4 - r4*s4
            x(k+2*m) = r2*c2 + i2*s2
            y(k+2*m) = i2*c2 - r2*s2
            x(k+3*m) = r6*c6 + i6*s6
            y(k+3*m) = i6*c6 - r6*s6
            x(k+4*m) = r1*c1 + i1*s1
            y(k+4*m) = i1*c1 - r1*s1
            x(k+5*m) = r5*c5 + i5*s5
            y(k+5*m) = i5*c5 - r5*s5
            x(k+6*m) = r3*c3 + i3*s3
            y(k+6*m) = i3*c3 - r3*s3
            x(k+7*m) = r7*c7 + i7*s7
            y(k+7*m) = i7*c7 - r7*s7
          end do
      
          if ((j > 0) .and. (2*j < m)) then
            t    = (c1+s1)*e
            s1   = (c1-s1)*e
            c1   = t
      
            t    = s2
            s2   = c2
            c2   = t
      
            t    = (-c3+s3)*e
            s3   = (c3+s3)*e
            c3   = t
      
            c4   = -c4
      
            t    = -(c5+s5)*e
            s5   = (-c5+s5)*e
            c5   = t
      
            t    = -s6
            s6   = -c6
            c6   = t
      
            t    = (c7-s7)*e
            s7   = -(c7+s7)*e
            c7   = t
            
            do k = m+1-j, npoints, 8*m
              rs0 = x(k    ) + x(k+4*m)
              is0 = y(k    ) + y(k+4*m)
              ru0 = x(k    ) - x(k+4*m)
              iu0 = y(k    ) - y(k+4*m)
              rs1 = x(k+  m) + x(k+5*m)
              is1 = y(k+  m) + y(k+5*m)
              ru1 = x(k+  m) - x(k+5*m)
              iu1 = y(k+  m) - y(k+5*m)
              rs2 = x(k+2*m) + x(k+6*m)
              is2 = y(k+2*m) + y(k+6*m)
              ru2 = x(k+2*m) - x(k+6*m)
              iu2 = y(k+2*m) - y(k+6*m)
              rs3 = x(k+3*m) + x(k+7*m)
              is3 = y(k+3*m) + y(k+7*m)
              ru3 = x(k+3*m) - x(k+7*m)
              iu3 = y(k+3*m) - y(k+7*m)
              
              rss0 = rs0 + rs2
              iss0 = is0 + is2
              rsu0 = rs0 - rs2
              isu0 = is0 - is2
              rss1 = rs1 + rs3
              iss1 = is1 + is3
              rsu1 = rs1 - rs3
              isu1 = is1 - is3
              rus0 = ru0 - iu2
              ius0 = iu0 + ru2
              ruu0 = ru0 + iu2
              iuu0 = iu0 - ru2
              rus1 = ru1 - iu3
              ius1 = iu1 + ru3
              ruu1 = ru1 + iu3
              iuu1 = iu1 - ru3
              t = (rus1+ius1)*e
              ius1 = (ius1-rus1)*e
              rus1 = t
              t = (ruu1+iuu1)*e
              iuu1 = (iuu1-ruu1)*e
              ruu1 = t
              
              r1 = ruu0 + ruu1
              i1 = iuu0 + iuu1
              r2 = rsu0 + isu1
              i2 = isu0 - rsu1
              r3 = rus0 + ius1
              i3 = ius0 - rus1
              r4 = rss0 - rss1
              i4 = iss0 - iss1
              r5 = ruu0 - ruu1
              i5 = iuu0 - iuu1
              r6 = rsu0 - isu1
              i6 = isu0 + rsu1
              r7 = rus0 - ius1
              i7 = ius0 + rus1
              
              x(k    ) = rss0 + rss1
              y(k    ) = iss0 + iss1
              x(k+  m) = r4*c4 + i4*s4
              y(k+  m) = i4*c4 - r4*s4
              x(k+2*m) = r2*c2 + i2*s2
              y(k+2*m) = i2*c2 - r2*s2
              x(k+3*m) = r6*c6 + i6*s6
              y(k+3*m) = i6*c6 - r6*s6
              x(k+4*m) = r1*c1 + i1*s1
              y(k+4*m) = i1*c1 - r1*s1
              x(k+5*m) = r5*c5 + i5*s5
              y(k+5*m) = i5*c5 - r5*s5
              x(k+6*m) = r3*c3 + i3*s3
              y(k+6*m) = i3*c3 - r3*s3
              x(k+7*m) = r7*c7 + i7*s7
              y(k+7*m) = i7*c7 - r7*s7
            end do
          end if
        end do
      
      end subroutine ecrc06
      
      subroutine ecsc06(x, y, m)
        integer,        intent(in)    :: m
        real(kind=dbl), intent(inout) :: x(:), y(:)
        integer                       :: j, k
        real(kind=dbl)                :: t, a1, a2, as, au, b1, b2, c1, c2, c3, c4, s1, s2, s3, s4, angle
        real(kind=dbl)                :: i0, i1, i2, i3, i4, ia1, ia2, ias, iau, ib1, ib2, is1, is2, iss, iu1, iu2
        real(kind=dbl)                :: r0, r1, r2, r3, r4, ra1, ra2, ras, rau, rb1, rb2, rs1, rs2, rss, ru1, ru2
      
        a1 = cos(2*pi/5)
        b1 = sin(2*pi/5)
        a2 = cos(4*pi/5)
        b2 = sin(4*pi/5)
        as = -0.25_dbl
        au = sqrt(5._dbl)/4
      
        do j = 0,  m/2
          angle = 2 * pi * j / m / 5
            c1 = cos(angle)
            s1 = sin(angle)
              c2 = c1*c1 - s1*s1
              s2 = s1*c1 + c1*s1
              c3 = c2*c1 - s2*s1
              s3 = s2*c1 + c2*s1
              c4 = c2*c2 - s2*s2
              s4 = s2*c2 + c2*s2
          
          do k = j+1, npoints, 5*m
            r0 = x(k)
            i0 = y(k)
            rs1 = x(k+m) + x(k+4*m)
            is1 = y(k+m) + y(k+4*m)
            ru1 = x(k+m) - x(k+4*m)
            iu1 = y(k+m) - y(k+4*m)
            rs2 = x(k+2*m) + x(k+3*m)
            is2 = y(k+2*m) + y(k+3*m)
            ru2 = x(k+2*m) - x(k+3*m)
            iu2 = y(k+2*m) - y(k+3*m)
            
            rss = rs1 + rs2
            iss = is1 + is2
            ras = r0 + rss*as
            ias = i0 + iss*as
            rau = (rs1-rs2)*au
            iau = (is1-is2)*au
            ra1 = ras + rau
            ia1 = ias + iau
            ra2 = ras - rau
            ia2 = ias - iau
            rb1 = ru1*b1 + ru2*b2
            ib1 = iu1*b1 + iu2*b2
            rb2 = ru1*b2 - ru2*b1
            ib2 = iu1*b2 - iu2*b1
            
            r1 = ra1 + ib1
            i1 = ia1 - rb1
            r2 = ra2 + ib2
            i2 = ia2 - rb2
            r3 = ra2 - ib2
            i3 = ia2 + rb2
            r4 = ra1 - ib1
            i4 = ia1 + rb1
            
            x(k    ) = r0 + rss
            y(k    ) = i0 + iss
            x(k+  m) = r1*c1 + i1*s1
            y(k+  m) = i1*c1 - r1*s1
            x(k+2*m) = r2*c2 + i2*s2
            y(k+2*m) = i2*c2 - r2*s2
            x(k+3*m) = r3*c3 + i3*s3
            y(k+3*m) = i3*c3 - r3*s3
            x(k+4*m) = r4*c4 + i4*s4
            y(k+4*m) = i4*c4 - r4*s4
          end do
      
          if ((j > 0) .and. (2*j < m)) then
            t  = c1*a1 + s1*b1
            s1 = c1*b1 - s1*a1
            c1 = t
      
            t  = c2*a2 + s2*b2
            s2 = c2*b2 - s2*a2
            c2 = t
      
            t  = c3*a2 - s3*b2
            s3 = -c3*b2 - s3*a2
            c3 = t
      
            t  = c4*a1 - s4*b1
            s4 = -c4*b1 - s4*a1
            c4 = t
      
            do k = m+1-j, npoints, 5*m
              r0 = x(k)
              i0 = y(k)
              rs1 = x(k+m) + x(k+4*m)
              is1 = y(k+m) + y(k+4*m)
              ru1 = x(k+m) - x(k+4*m)
              iu1 = y(k+m) - y(k+4*m)
              rs2 = x(k+2*m) + x(k+3*m)
              is2 = y(k+2*m) + y(k+3*m)
              ru2 = x(k+2*m) - x(k+3*m)
              iu2 = y(k+2*m) - y(k+3*m)
              
              rss = rs1 + rs2
              iss = is1 + is2
              ras = r0 + rss*as
              ias = i0 + iss*as
              rau = (rs1-rs2)*au
              iau = (is1-is2)*au
              ra1 = ras + rau
              ia1 = ias + iau
              ra2 = ras - rau
              ia2 = ias - iau
              rb1 = ru1*b1 + ru2*b2
              ib1 = iu1*b1 + iu2*b2
              rb2 = ru1*b2 - ru2*b1
              ib2 = iu1*b2 - iu2*b1
              
              r1 = ra1 + ib1
              i1 = ia1 - rb1
              r2 = ra2 + ib2
              i2 = ia2 - rb2
              r3 = ra2 - ib2
              i3 = ia2 + rb2
              r4 = ra1 - ib1
              i4 = ia1 + rb1
              
              x(k    ) = r0 + rss
              y(k    ) = i0 + iss
              x(k+  m) = r1*c1 + i1*s1
              y(k+  m) = i1*c1 - r1*s1
              x(k+2*m) = r2*c2 + i2*s2
              y(k+2*m) = i2*c2 - r2*s2
              x(k+3*m) = r3*c3 + i3*s3
              y(k+3*m) = i3*c3 - r3*s3
              x(k+4*m) = r4*c4 + i4*s4
              y(k+4*m) = i4*c4 - r4*s4
            end do
          end if
        end do
      
      end subroutine ecsc06
      
      subroutine ectc06(x, y, m)
        integer,        intent(in)    :: m
        real(kind=dbl), intent(inout) :: x(:), y(:)
        integer                       :: j, k
        real(kind=dbl)                :: t, i1, i2, i3, is0, is1, iu0, iu1, c1, c2, c3,     &
                                      &     r1, r2, r3, rs0, rs1, ru0, ru1, s1, s2, s3, ang
      
        do j = 0, m/2
          ang = pi * j / m / 2
            c1 = cos(ang)
            s1 = sin(ang)
              c2 = c1*c1 - s1*s1
              s2 = s1*c1 + c1*s1
              c3 = c2*c1 - s2*s1
              s3 = s2*c1 + c2*s1
      
          do k = j+1, npoints, 4*m
            rs0 = x(k) + x(k+2*m)
            is0 = y(k) + y(k+2*m)
            ru0 = x(k) - x(k+2*m)
            iu0 = y(k) - y(k+2*m)
            rs1 = x(k+m) + x(k+3*m)
            is1 = y(k+m) + y(k+3*m)
            ru1 = x(k+m) - x(k+3*m)
            iu1 = y(k+m) - y(k+3*m)
            
            r1 = ru0 + iu1
            i1 = iu0 - ru1
            r2 = rs0 - rs1
            i2 = is0 - is1
            r3 = ru0 - iu1
            i3 = iu0 + ru1
            
            x(k    ) = rs0 + rs1
            y(k    ) = is0 + is1
            x(k+  m) = r2*c2 + i2*s2
            y(k+  m) = i2*c2 - r2*s2
            x(k+2*m) = r1*c1 + i1*s1
            y(k+2*m) = i1*c1 - r1*s1
            x(k+3*m) = r3*c3 + i3*s3
            y(k+3*m) = i3*c3 - r3*s3
          end do
      
          if ((j > 0) .and. (2*j < m)) then
            t = c1
            c1 = s1
            s1 = t
      
            c2 = -c2
      
            t = c3
            c3 = -s3
            s3 = -t
            
            do k = m+1-j, npoints, 4*m
              rs0 = x(k) + x(k+2*m)
              is0 = y(k) + y(k+2*m)
              ru0 = x(k) - x(k+2*m)
              iu0 = y(k) - y(k+2*m)
              rs1 = x(k+m) + x(k+3*m)
              is1 = y(k+m) + y(k+3*m)
              ru1 = x(k+m) - x(k+3*m)
              iu1 = y(k+m) - y(k+3*m)
              
              r1 = ru0 + iu1
              i1 = iu0 - ru1
              r2 = rs0 - rs1
              i2 = is0 - is1
              r3 = ru0 - iu1
              i3 = iu0 + ru1
              
              x(k    ) = rs0 + rs1
              y(k    ) = is0 + is1
              x(k+  m) = r2*c2 + i2*s2
              y(k+  m) = i2*c2 - r2*s2
              x(k+2*m) = r1*c1 + i1*s1
              y(k+2*m) = i1*c1 - r1*s1
              x(k+3*m) = r3*c3 + i3*s3
              y(k+3*m) = i3*c3 - r3*s3
            end do
          end if
        end do
      
      end subroutine ectc06
      
      subroutine ecuc06(x, y, m)
        integer,        intent(in)    :: m
        real(kind=dbl), intent(inout) :: x(:), y(:)
        integer                       :: j, k
        real(kind=dbl)                :: a, b, t, c1, c2, s1, s2, i0, i1, i2, ia, ib, is, r0, r1, r2, ra, rb, rs, angle
      
        a = -0.5_dbl; b = sqrt(0.75_dbl)
      
        do j = 0, m/2
          angle = 2 * pi * j / m / 3
            c1 = cos(angle)
            s1 = sin(angle)
          
          c2 = c1*c1 - s1*s1
          s2 = s1*c1 + c1*s1
          
          do k = j+1, npoints, 3*m
            r0 = x(k)
            i0 = y(k)
            rs = x(k+m) + x(k+2*m)
            is = y(k+m) + y(k+2*m)
      
            x(k) = r0 + rs
            y(k) = i0 + is
      
            ra = r0 + rs*a
            ia = i0 + is*a
            rb = (x(k+m) - x(k+2*m))*b
            ib = (y(k+m) - y(k+2*m))*b
      
            r1 = ra + ib
            i1 = ia - rb
            r2 = ra - ib
            i2 = ia + rb
      
            x(k+m) = r1*c1 + i1*s1
            y(k+m) = i1*c1 - r1*s1
            x(k+2*m) = r2*c2 + i2*s2
            y(k+2*m) = i2*c2 - r2*s2
          end do
      
          if ((j > 0) .and. (2*j < m)) then
            t  = c1*a + s1*b
            s1 = c1*b - s1*a
            c1 = t
      
            t  = +c2*a - s2*b
            s2 = -c2*b - s2*a
            c2 = t
            
            do k = m+1-j, npoints, 3*m
              r0 = x(k)
              i0 = y(k)
              rs = x(k+m) + x(k+2*m)
              is = y(k+m) + y(k+2*m)
      
              x(k) = r0 + rs
              y(k) = i0 + is
      
              ra = r0 + rs*a
              ia = i0 + is*a
              rb = (x(k+m) - x(k+2*m))*b
              ib = (y(k+m) - y(k+2*m))*b
      
              r1 = ra + ib
              i1 = ia - rb
              r2 = ra - ib
              i2 = ia + rb
      
              x(k+m) = r1*c1 + i1*s1
              y(k+m) = i1*c1 - r1*s1
              x(k+2*m) = r2*c2 + i2*s2
              y(k+2*m) = i2*c2 - r2*s2
            end do
          end if
        end do
      
      end subroutine ecuc06
      
      subroutine ecvc06(x, y, m)
        integer,        intent(in)    :: m
        real(kind=dbl), intent(inout) :: x(:), y(:)
        integer                       :: j, k
        real(kind=dbl)                :: c, s, is, iu, rs, ru, angle
      
        do j = 0, m/2
          angle = pi * j / m
            c = cos(angle)
            s = sin(angle)
      
          do k = j+1, npoints, 2*m
            rs = x(k) + x(k+m)
            is = y(k) + y(k+m)
            ru = x(k) - x(k+m)
            iu = y(k) - y(k+m)
      
            x(k) = rs
            y(k) = is
      
            x(k+m) = ru*c + iu*s
            y(k+m) = iu*c - ru*s
          end do
      
          if ((j > 0) .and. (2*j < m)) then
            c = -c
            
            do k = m+1-j, npoints, 2*m
              rs = x(k) + x(k+m)
              is = y(k) + y(k+m)
              ru = x(k) - x(k+m)
              iu = y(k) - y(k+m)
      
              x(k) = rs
              y(k) = is
      
              x(k+m) = ru*c + iu*s
              y(k+m) = iu*c - ru*s
            end do
          end if
        end do
      
      end subroutine ecvc06
      
    subroutine ecyc06(x, y)
      real(kind=dbl), intent(inout) :: x(:), y(:)
      integer                       :: dk, i, ii, il, j, jj, jl, k, kk, ks, lk, mods, mult, punsym, test, varmodulo(20)
      real(kind=dbl)                :: t
      
      call ecxc06(x, y)
  
      if (unsym(1) /= 0) then
        punsym = npoints/psym**2
        mult   = punsym/unsym(1)
        test   = (unsym(1)*unsym(2)-1)*mult*psym
        
        lk = mult
        dk = mult
        
        do k = 2, 20
          if (unsym(k) == 0) exit
        
          lk = lk * unsym(k-1)
          dk = dk / unsym(k)
        
          varmodulo(k) = (lk-dk)*psym; mods = k
        end do
        
        if (mods >= 3) then
          k = (mods+3)/2
          
          do j = 3, k
            kk                  = varmodulo(j)
            varmodulo(j)        = varmodulo(mods+3-j)
            varmodulo(mods+3-j) = kk
          end do
        end if
        
        jl = (punsym-3) * psym
        ks = (punsym  ) * psym
        
        do j = psym, jl, psym
          jj = j
        
          do
            jj = jj*mult
            
            if (mods >= 3) then
              do i = 3, mods
                jj = jj - (jj/varmodulo(i))*varmodulo(i)
              end do
            end if
        
            if (jj < test) then
              jj = jj - (jj/varmodulo(2))*varmodulo(2)
            else
              jj = jj - (jj/varmodulo(2))*varmodulo(2) + varmodulo(2)
            end if
        
            if (jj >= j) exit
          end do
        
          if (jj /= j) then
            lk = jj - j
            ii = j + 1
            il = j + psym
        
            do i = ii, il
              do k = i, npoints, ks
                kk = k + lk
        
                t     = x(k)
                x(k)  = x(kk)
                x(kk) = t
        
                t     = y(k)
                y(k)  = y(kk)
                y(kk) = t
              end do
            end do
          end if
        end do
      end if
      
    end subroutine ecyc06
  
      subroutine ecxc06(x, y)
        real(kind=dbl), intent(inout) :: x(:), y(:)
        real(kind=dbl)                :: t
        integer                       :: j, kk, level, i(20), k(20), l(20)
        integer                       :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, &
                                      & k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, &
                                      & l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
    
        if (sym(1) /= 0) then
          l = 1; i = 1
          
          l(1) = npoints
          i(1) = l(1) / sym(1)
          
          do j = 2, 20
            if (sym(j) == 0) exit
          
            l(j) = l(j-1) / sym(j-1)
            i(j) = l(j)   / sym(j)
          end do
          
          kk = 0; level = 20; k(level) = 1
          
          l1 = l(1); l2 = l(2); l3 = l(3); l4 = l(4); l5 = l(5); l6 = l(6); l7 = l(7); l8 = l(8); l9 = l(9); l10 = l(10)
          i1 = i(1); i2 = i(2); i3 = i(3); i4 = i(4); i5 = i(5); i6 = i(6); i7 = i(7); i8 = i(8); i9 = i(9); i10 = i(10)
          
          do
            level = level - 1
            
            do j = 10, level
              k(level+10-j) = k(level+11-j)
            end do
            
            do k10 = k(11), l10, i10
              k(10) = k10
              do k9 = k10, l9, i9
                k(9) = k9
                do k8 = k9, l8, i8
                  k(8) = k8
                  do k7 = k8, l7, i7
                    k(7) = k7
                    do k6 = k7, l6, i6
                      k(6) = k6
                      do k5 = k6, l5, i5
                        k(5) = k5
                        do k4 = k5, l4, i4
                          k(4) = k4
                          do k3 = k4, l3, i3
                            k(3) = k3
                            do k2 = k3, l2, i2
                              k(2) = k2
                              do k1 = k2, l1, i1
                                k(1) = k1
        
                                kk = kk + 1
                                if (kk < k1) then
                                  t     = x(kk)
                                  x(kk) = x(k1)
                                  x(k1) = t
        
                                  t     = y(kk)
                                  y(kk) = y(k1)
                                  y(k1) = t
                                end if
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
            
            level = 10
            
            do
              if (level >= 20) return
              
              level = level + 1
              k(level) = k(level) + i(level)
              if (k(level) <= l(level)) exit
            end do
            
          end do
        end if
        
      end subroutine ecxc06
      
end module FFT
