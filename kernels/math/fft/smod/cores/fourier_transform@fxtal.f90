submodule (fourier_transform) fxtal
  implicit none; contains
  
  module procedure fxztal
    integer :: itsum, it1, l, ip, isd, ld, k1, icdd, it1d
    
    l     = this%n/2
    it1   = this%it(l-1)
    itsum = this%it(l)
    ip    = mod(it1,4)+2
    
    select case (ip)
      case (4)
        call fxzm4b(m, l, x)
      case (2)
        call fxzm2b(m, l, x)
      case (3)
        call fxzm3b(m, l, x)
      case (5)
        call fxzm5b(m, l, x)
    end select
    
    k1   = ip
    ld   = l / ip
    isd  = ip-1
    it1d = it1 / 4
    
    do icdd = 2, itsum
      ip = mod(it1d,4)+2
      
      select case (ip)
        case (4)
          call fxzm4a(m, k1, ld, x, this%t(1+2*isd))
        case (2)
          call fxzm2a(m, k1, ld, x, this%t(1+2*isd))
        case (3)
          call fxzm3a(m, k1, ld, x, this%t(1+2*isd))
        case (5)
          call fxzm5a(m, k1, ld, x, this%t(1+2*isd))
      end select
      
      isd  = isd + k1 * (ip-1)
      k1   = k1 * ip
      ld   = ld / ip
      it1d = it1d / 4
    end do
    
  end procedure fxztal
  
end submodule fxtal