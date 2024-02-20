submodule (SphericalHarmonics) lege_setup
  implicit none; contains
  
  module pure subroutine lege_setup_2_sub(this, it, cosx, sinx, weight)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: it
    real(kind=dbl),       intent(out) :: cosx(*), sinx(*), weight(*)
    integer                           :: i
    
    do concurrent ( i = 1:2 )
      cosx(i)   = this%roots(it+i)
      sinx(i)   = sqrt(1-cosx(i)**2)
      weight(i) = this%fftLege(it+i)
    end do
    
  end subroutine lege_setup_2_sub
  
  module pure subroutine lege_setup_4_sub(this, it, cosx, sinx, weight)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: it
    real(kind=dbl),       intent(out) :: cosx(*), sinx(*), weight(*)
    integer                           :: i
    
    do concurrent ( i = 1:4 )
      cosx(i)   = this%roots(it+i)
      sinx(i)   = sqrt(1-cosx(i)**2)
      weight(i) = this%fftLege(it+i)
    end do
    
  end subroutine lege_setup_4_sub
  
  module pure subroutine lege_setup_8_sub(this, it, cosx, sinx, weight)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: it
    real(kind=dbl),       intent(out) :: cosx(*), sinx(*), weight(*)
    integer                           :: i
    
    do concurrent ( i = 1:8 )
      cosx(i)   = this%roots(it+i)
      sinx(i)   = sqrt(1-cosx(i)**2)
      weight(i) = this%fftLege(it+i)
    end do
    
  end subroutine lege_setup_8_sub
  
  module pure subroutine lege_setup_16_sub(this, it, cosx, sinx, weight)
    class(T_lateralGrid), intent(in)  :: this
    integer,              intent(in)  :: it
    real(kind=dbl),       intent(out) :: cosx(*), sinx(*), weight(*)
    integer                           :: i
    
    do concurrent ( i = 1:16 )
      cosx(i)   = this%roots(it+i)
      sinx(i)   = sqrt(1-cosx(i)**2)
      weight(i) = this%fftLege(it+i)
    end do
    
  end subroutine lege_setup_16_sub
  
end submodule lege_setup