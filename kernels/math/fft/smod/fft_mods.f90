submodule(Fourier_transform) fft_mods
  implicit none
  
  real(kind=dbl), parameter :: C31 = -0.5_dbl
  real(kind=dbl), parameter :: C32 = +0.86602540378443864676_dbl
  real(kind=dbl), parameter :: C51 = +0.25_dbl
  real(kind=dbl), parameter :: C52 = +0.5590169943749474241_dbl
  real(kind=dbl), parameter :: C53 = +0.6180339887498948482_dbl
  real(kind=dbl), parameter :: C54 = -0.9510565162951535721_dbl
  
  contains
  
  module pure subroutine fxzini(n, it, t)
    integer,        intent(in)  :: n
    integer,        intent(out) :: it(n)
    real(kind=dbl), intent(out) :: t(2,0:n-1)
    integer                     :: j, k, l, i, is, ic, i4, it1, isj, ipc, ip, ipd
    integer,        allocatable :: ipn(:), itc(:), itw(:), iw(:)
    
    allocate( ipn(4), itc(2:5), itw(0:n-1) , iw(n-2) ) ; ipn = [5,4,3,2] ; t = 0._dbl ; iw = 0 ; itw = 0
    
    itc = 0 ; j = n
      do i = 5, 2, -1
        itc(i) = 0 ; k = 0
        
        do while ( k == 0 )
          k = mod( j, i )
          
          if (k == 0) then
            itc(i) = itc(i) + 1
            j = j / i
          end if
        end do
      end do
    
    itw = 0 ; k = n ; l = 1
      do ipc = 4, 1, -1
        ip = ipn(ipc)
        
        do ic = 1, itc(ip)
          do concurrent ( i = 0:l-1 , j = 0:k/ip-1 , ipd = 1:ip-1 )
            itw(i+l*(ipd+j*ip)) = itw(i+l*(ipd+j*ip)) + ipd * k/ip
          end do
          
          k = k / ip
          l = l * ip
        end do
      end do
    
    k = n ; l = 1 ; is = 0
      do ipc = 1, 4
        ip = ipn(ipc)
        
        do ic = 1, itc(ip)
          do i = 0, l-1
            j = itw(k*i)
            
            do ipd = 1, ip-1
              t(1,is+(ip-1)*i+ipd-1) = cos( 2 * pi * j * ipd / ( ip * l ) )
              t(2,is+(ip-1)*i+ipd-1) = sin( 2 * pi * j * ipd / ( ip * l ) )
            end do
          end do
          
          is = is + l * (ip-1)        
          k  = k / ip
          l  = l * ip
        end do
      end do
    
    itw = 0 ; k = n ; l = 1
      do ipc = 1, 4
        ip = ipn(ipc)
        
        do ic = 1, itc(ip)
          do concurrent ( i = 0:l-1 , j = 0:k/ip-1 , ipd = 1:ip-1 )
            itw(i+l*(ipd+j*ip)) = itw(i+l*(ipd+j*ip)) + ipd * k/ip
          end do
          
          k = k / ip
          l = l * ip
        end do
      end do
    
    it1 = 0 ; i4 = 1
      do ipc = 1, 4
        ip = ipn(ipc)
        
        do ic = 1, itc(ip)
          it1 = it1 + (ip-2) * i4 ; i4 = i4 * 4
        end do
      end do
    
    it(n-1) = it1 ; it(n) = sum( itc ) ; ic = 0
      do j = 1, n-2
        if ( iw(j) /= 1 ) then
          ic = ic + 1 ; it(ic) = j ; iw(j) = 1 ; isj = itw(j)
            do
              if ( isj == j ) then
                it(ic) = it(ic) + imm ; exit
              end if
              
              ic      = ic + 1
              it(ic)  = isj
              iw(isj) = 1             
              isj     = itw(isj)
            end do
        end if
      end do
      
    deallocate( ipn, itc, itw , iw )
    
  end subroutine fxzini
  
  module pure recursive subroutine fxztal(m, k, l, x, t, ic, itsum, is, j1, it1)
    integer,        intent(in)    :: m, k, l, ic, itsum, is, j1, it1
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,0:*)
    integer                       :: ip, j, isd, kd, ld, k1, jd, icd, icdd, it1d
    
    kd = k ; ld = l ; k1 = 1 ; j = j1 ; it1d = it1 ; ip = mod(it1d,4)+2
    
    if ( ic /= 1 ) then
      select case (ip)
        case (4)
          call fxzm4a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
        case (2)
          call fxzm2a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
        case (3)
          call fxzm3a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
        case (5)
          call fxzm5a(m, k1, ld, x(1,1,j*ld), t(1,is+j*(ip-1)))
      end select
    end if
    
    k1 = k1 * ip ; ld = ld / ip ; isd = is + kd * (ip-1) ; kd = kd * ip ; j = j * ip ; it1d = it1d / 4 ; icd = ic + 1
    
    if ( (m * ip * ld > 2048) .and. (icd <= itsum) ) then
      do jd = 0, ip-1
         call fxztal(m, kd, ld, x, t, icd, itsum, isd, j+jd, it1d)
      end do
    else
      do icdd = icd, itsum
        ip = mod(it1d,4)+2
        
        select case (ip)
          case (4)
            call fxzm4a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
          case (2)
            call fxzm2a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
          case (3)
            call fxzm3a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
          case (5)
            call fxzm5a(m, k1, ld, x(1,1,j*ld), t(1,isd+j*(ip-1)))
        end select
        
        k1 = k1 * ip ; ld = ld/ip ; isd = isd + kd * (ip-1) ; kd = kd * ip ; j = j * ip ; it1d = it1d / 4
      end do
    end if
    
  end subroutine fxztal
  
  module pure subroutine fxzm2a(m, k, l, x, t)
    integer,        intent(in)    :: m, k, l
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,l/2,0:1,0:k-1)
    integer                       :: i, j, iv
    real(kind=dbl)                :: x1re, x1im
  
    do j = 0, k-1
       do concurrent ( i = 1:L/2 , iv = 1:m )
         x1re = x(iv,1,i,0,j) - t(1,j) * x(iv,1,i,1,j) ; x1im = x(iv,2,i,0,j) - t(2,j) * x(iv,1,i,1,j)
         
         x(iv,1,i,1,j) =     x1re          + t(2,j) * x(iv,2,i,1,j) ; x(iv,2,i,1,j) =     x1im          - t(1,j) * x(iv,2,i,1,j)
         x(iv,1,i,0,j) = 2 * x(iv,1,i,0,j) -          x(iv,1,i,1,j) ; x(iv,2,i,0,j) = 2 * x(iv,2,i,0,j) -          x(iv,2,i,1,j)
       end do
    end do
  
  end subroutine fxzm2a
  
  module pure subroutine fxzm2b(m, l, x)
    integer,        intent(in)    :: m, l
    real(kind=dbl), intent(inout) :: x(m,2,l/2,0:1)
    integer                       :: i, iv
    
    do concurrent ( i = 1:l/2 , iv = 1:m )
      x(iv,1,i,1) =     x(iv,1,i,0) - x(iv,1,i,1) ; x(iv,2,i,1) =     x(iv,2,i,0) - x(iv,2,i,1)
      x(iv,1,i,0) = 2 * x(iv,1,i,0) - x(iv,1,i,1) ; x(iv,2,i,0) = 2 * x(iv,2,i,0) - x(iv,2,i,1)
    end do
    
  end subroutine fxzm2b
  
  module pure subroutine fxzm3a(m, k, l, x, t)
    integer,        intent(in)    :: m, k, l
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,l/3,0:2,0:k-1)
    integer                       :: i, j, ij, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im
    real(kind=dbl)                :: t1re, t1im, t2re, t2im
    
    ij = 0
    
    do j = 0, k-1
       t1re = t(1,ij  ) ; t1im = t(2,ij  )            
       t2re = t(1,ij+1) ; t2im = t(2,ij+1)
                   
       do concurrent ( i = 1:L/3 , iv = 1:m )
         x0re =        t1re * x(iv,1,i,1,j) - t1im * x(iv,2,i,1,j) ; x0im =        t1re * x(iv,2,i,1,j) + t1im * x(iv,1,i,1,j)
         x1re = x0re - t2re * x(iv,1,i,2,j) + t2im * x(iv,2,i,2,j) ; x1im = x0im - t2re * x(iv,2,i,2,j) - t2im * x(iv,1,i,2,j)
         
         x0re = 2 * x0re          -       x1re ; x0im = 2 * x0im          -       x1im
         x2re =     x(iv,1,i,0,j) + C31 * x0re ; x2im =     x(iv,2,i,0,j) + C31 * x0im
         
         x(iv,1,i,0,j) =     x0re +       x(iv,1,i,0,j) ; x(iv,2,i,0,j) =     x0im +       x(iv,2,i,0,j)
         x(iv,1,i,2,j) =     x2re + C32 * x1im          ; x(iv,2,i,2,j) =     x2im - C32 * x1re
         x(iv,1,i,1,j) = 2 * x2re -       x(iv,1,i,2,j) ; x(iv,2,i,1,j) = 2 * x2im -       x(iv,2,i,2,j)
       end do
       
       ij = ij + 2
    end do
    
  end subroutine fxzm3a
  
  module pure subroutine fxzm3b(m, l, x)
    integer,        intent(in)    :: m, l
    real(kind=dbl), intent(inout) :: x(m,2,l/3,0:2)
    integer                       :: i, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im
    
    do concurrent ( i = 1:l/3 , iv = 1:m )
      x1re = x(iv,1,i,1) -       x(iv,1,i,2) ; x1im = x(iv,2,i,1) -       x(iv,2,i,2)
      x0re = x(iv,1,i,1) +       x(iv,1,i,2) ; x0im = x(iv,2,i,1) +       x(iv,2,i,2)
      x2re = x(iv,1,i,0) + C31 * x0re        ; x2im = x(iv,2,i,0) + C31 * x0im
      
      x(iv,1,i,0) =     x0re +       x(iv,1,i,0) ; x(iv,2,i,0) =     x0im +       x(iv,2,i,0)
      x(iv,1,i,2) =     x2re + C32 * x1im        ; x(iv,2,i,2) =     x2im - C32 * x1re
      x(iv,1,i,1) = 2 * x2re -       x(iv,1,i,2) ; x(iv,2,i,1) = 2 * x2im -       x(iv,2,i,2)
    end do
    
  end subroutine fxzm3b
  
  module pure subroutine fxzm4a(m, k, l, x, t)
    integer,        intent(in)    :: m, k, l
    REAL(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,l/4,0:3,0:k-1)
    integer                       :: i, j, ij, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im
    real(kind=dbl)                :: t1re, t1im, t2re, t2im
    
    ij = 0
    
    do j = 0, k-1
      t1re = t(1,ij  ) ; t1im = t(2,ij  )            
      t2re = t(1,ij+1) ; t2im = t(2,ij+1)
      
      do concurrent ( i = 1:l/4 , iv = 1:m )
        x2re = ( x(iv,1,i,0,j) - t2re * x(iv,1,i,2,j) ) + t2im * x(iv,2,i,2,j)
        x2im = ( x(iv,2,i,0,j) - t2im * x(iv,1,i,2,j) ) - t2re * x(iv,2,i,2,j)
        x3re = ( x(iv,1,i,1,j) - t2re * x(iv,1,i,3,j) ) + t2im * x(iv,2,i,3,j)
        x3im = ( x(iv,2,i,1,j) - t2im * x(iv,1,i,3,j) ) - t2re * x(iv,2,i,3,j)
        
        x0re = 2 * x(iv,1,i,0,j) - x2re ; x0im = 2 * x(iv,2,i,0,j) - x2im
        x1re = 2 * x(iv,1,i,1,j) - x3re ; x1im = 2 * x(iv,2,i,1,j) - x3im
        
        x(iv,1,i,2,j) = ( x0re - t1re * x1re ) + t1im * x1im ; x(iv,2,i,2,j) = ( x0im - t1im * x1re ) - t1re * x1im
        x(iv,1,i,0,j) = 2 * x0re - x(iv,1,i,2,j)             ; x(iv,2,i,0,j) = 2 * x0im - x(Iv,2,i,2,j)
        x(iv,1,i,1,j) = ( x2re - t1re * x3im ) - t1im * x3re ; x(iv,2,i,1,j) = ( x2im + t1re * x3re ) - t1im * x3im
        x(iv,1,i,3,j) = 2 * x2re - x(iv,1,i,1,j)             ; x(iv,2,i,3,j) = 2 * x2im - x(iv,2,i,1,j)
      end do
      
      ij = ij + 3
    end do
    
  end subroutine fxzm4a
  
  module pure subroutine fxzm4b(m, l, x)
    integer,        intent(in)    :: m, l
    real(kind=dbl), intent(inout) :: x(m,2,l/4,0:3)
    integer                       :: i, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im
  
    do concurrent ( i = 1:l/4 , iv = 1:m )
      x2re = x(iv,1,i,0) - x(iv,1,i,2) ; x2im = x(iv,2,i,0) - x(iv,2,i,2)
      x0re = x(iv,1,i,0) + x(iv,1,i,2) ; x0im = x(iv,2,i,0) + x(iv,2,i,2)
      x3re = x(iv,1,i,1) - x(iv,1,i,3) ; x3im = x(iv,2,i,1) - x(iv,2,i,3)
      x1re = x(iv,1,i,1) + x(iv,1,i,3) ; x1im = x(iv,2,i,1) + x(iv,2,i,3)
      
      x(iv,1,i,2) =     x0re - x1re        ; x(iv,2,i,2) =     x0im - x1im
      x(iv,1,i,0) = 2 * x0re - x(iv,1,i,2) ; x(iv,2,i,0) = 2 * x0im - x(iv,2,i,2)       
      x(iv,1,i,1) =     x2re - x3im        ; x(iv,2,i,1) =     x2im + x3re
      x(iv,1,i,3) = 2 * x2re - x(iv,1,i,1) ; x(iv,2,i,3) = 2 * x2im - x(iv,2,i,1)       
    end do
  
  end subroutine fxzm4b
  
  module pure subroutine fxzm5a(m, k, l, x, t)
    integer,        intent(in)    :: m, k, l
    real(kind=dbl), intent(in)    :: t(2,0:*)
    real(kind=dbl), intent(inout) :: x(m,2,l/5,0:4,0:k-1)
    integer                       :: i, j, ij, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, x4re, x4im
    real(kind=dbl)                :: t1re, t1im, t2re, t2im, t3re, t3im, t4re, t4im
    
    ij = 0
    
    do j = 0, k-1
       t1re = t(1,ij  ) ; t1im = t(2,ij  )
       t2re = t(1,ij+1) ; t2im = t(2,ij+1)
       t3re = t(1,ij+2) ; t3im = t(2,ij+2)
       t4re = t(1,ij+3) ; t4im = t(2,ij+3)
  
       do concurrent ( i = 1:l/5 , iv = 1:m )
         x1re =        t1re * x(iv,1,i,1,j) - t1im * x(iv,2,i,1,j) ; x1im =        t1re * x(iv,2,i,1,j) + t1im * x(iv,1,i,1,j)
         x2re =        t2re * x(iv,1,i,2,j) - t2im * x(iv,2,i,2,j) ; x2im =        t2re * x(iv,2,i,2,j) + t2im * x(iv,1,i,2,j)
         x3re = x2re - t3re * x(iv,1,i,3,j) + t3im * x(iv,2,i,3,j) ; x3im = x2im - t3re * x(iv,2,i,3,j) - t3im * x(iv,1,i,3,j)
         x0re = x1re - t4re * x(iv,1,i,4,j) + t4im * x(iv,2,i,4,j) ; x0im = x1im - t4re * x(iv,2,i,4,j) - t4im * x(iv,1,i,4,j)
         
         x1re =  2  * x1re -       x0re ; x1im =  2  * x1im -       x0im
         x4re =  2  * x2re -       x3re ; x4im =  2  * x2im -       x3im
         x2re =       x0re + C53 * x3re ; x2im =       x0im + C53 * x3im
         x3re = C53 * x0re -       x3re ; x3im = C53 * x0im -       x3im
         x0re =       x1re +       x4re ; x0im =       x1im +       x4im
         
         x1re = x1re          -       x4re ; x1im = x1im          -       x4im
         x4re = x(iv,1,i,0,j) - C51 * x0re ; x4im = x(iv,2,i,0,j) - C51 * x0im
         
         x1re =     x4re - C52 * x1re ; x1im =     x4im - C52 * x1im
         x4re = 2 * x4re -       x1re ; x4im = 2 * x4im -       x1im
         
         x(iv,1,i,0,j) =     x(iv,1,i,0,j) +       x0re          ; x(iv,2,i,0,j) =     x(iv,2,i,0,j) +       x0im
         x(iv,1,i,3,j) =     x1re          - C54 * x3im          ; x(iv,2,i,3,j) =     x1im          + C54 * x3re
         x(iv,1,i,2,j) = 2 * x1re          -       x(iv,1,i,3,j) ; x(iv,2,i,2,j) = 2 * x1im          -       x(iv,2,i,3,j)
         x(iv,1,i,4,j) =     x4re          - C54 * x2im          ; x(iv,2,i,4,j) =     x4im          + C54 * x2re
         x(iv,1,i,1,j) = 2 * x4re          -       x(iv,1,i,4,j) ; x(iv,2,i,1,j) = 2 * x4im          -       x(iv,2,i,4,j)
       end do
       
       ij = ij + 4
    end do
  
  end subroutine fxzm5a
  
  module pure subroutine fxzm5b(m, l, x)
    integer,        intent(in)    :: m, l
    real(kind=dbl), intent(inout) :: x(m,2,l/5,0:4)
    integer                       :: i, iv
    real(kind=dbl)                :: x0re, x0im, x1re, x1im, x2re, x2im, x3re, x3im, x4re, x4im
    
    do concurrent (i = 1:l/5 , iv = 1:m)
      x0re = x(iv,1,i,1) - x(iv,1,i,4) ; x0im = x(iv,2,i,1) - x(iv,2,i,4)
      x1re = x(iv,1,i,1) + x(iv,1,i,4) ; x1im = x(iv,2,i,1) + x(iv,2,i,4)
      x3re = x(iv,1,i,2) - x(iv,1,i,3) ; x3im = x(iv,2,i,2) - x(iv,2,i,3)
      x4re = x(iv,1,i,2) + x(iv,1,i,3) ; x4im = x(iv,2,i,2) + x(iv,2,i,3)
      
      x2re =       x0re + C53 * x3re ; x2im =       x0im + C53 * x3im
      x3re = C53 * x0re -       x3re ; x3im = C53 * x0im -       x3im
      x0re =       x1re +       x4re ; x0im =       x1im +       x4im
      x1re =       x1re -       x4re ; x1im =       x1im -       x4im
      
      x4re =     x(iv,1,i,0) - C51 * x0re ; x4im =     x(iv,2,i,0) - C51 * x0im
      x1re =     x4re        - C52 * x1re ; x1im =     x4im        - C52 * x1im
      x4re = 2 * x4re        -       x1re ; x4im = 2 * x4im        -       x1im
      
      x(iv,1,i,0) =     x(iv,1,i,0) +       x0re        ; x(iv,2,i,0) =     x(iv,2,i,0) +       x0im
      x(iv,1,i,3) =     x1re        - C54 * x3im        ; x(iv,2,i,3) =     x1im        + C54 * x3re
      x(iv,1,i,2) = 2 * x1re        -       x(iv,1,i,3) ; x(iv,2,i,2) = 2 * x1im        -       x(iv,2,i,3)
      x(iv,1,i,4) =     x4re        - C54 * x2im        ; x(iv,2,i,4) =     x4im        + C54 * x2re
      x(iv,1,i,1) = 2 * x4re        -       x(iv,1,i,4) ; x(iv,2,i,1) = 2 * x4im        -       x(iv,2,i,4)
    end do
    
  end subroutine fxzm5b
  
end submodule fft_mods