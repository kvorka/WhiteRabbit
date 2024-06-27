submodule (Fourier_transform) fxinit
  implicit none; contains
  
  module pure subroutine fxzini(n, it, t)
    integer,        intent(in)  :: n
    integer,        intent(out) :: it(n)
    real(kind=dbl), intent(out) :: t(2,0:n-1)
    integer                     :: j, k, l, i, is, ic, i4, it1, isj, ipc, ip, ipd
    integer,        allocatable :: ipn(:), itc(:), itw(:), iw(:)
    
    allocate( ipn(4), itc(2:5), itw(0:n-1) , iw(n-2) ) ; ipn = [5,4,3,2] ; t = zero ; iw = 0 ; itw = 0
    
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
  
end submodule fxinit