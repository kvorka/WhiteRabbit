submodule (SphericalHarmonics) Init_SphericalHarmonics
  implicit none ; contains
  
  module subroutine init_harmonics_sub(this, jmax)
    class(T_lateralGrid), intent(inout) :: this
    integer,              intent(in)    :: jmax
    integer                             :: i, k, j, m, n, ncnt
    real(kind=dbl)                      :: xincr, x, y, fx, fy
    
    if ( .not. ( any( addmissible_jmax == jmax ) ) ) then
      write(*,*) 'Due to FFT, this value of jmax is prohibited. Please, see table of admissible values in SphericalHarmonics.f90'
      stop
    end if
    
    this%jmax = jmax
    this%jms  =     ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jms1 =     ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    this%jms2 =     ((jmax+2)*(jmax+3)/2 + (jmax+2)) + 1
    this%jmv  = 3 * ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jmv1 = 3 * ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    
    this%maxj      = jmax+2
    this%nFourier  = 3*(this%maxj+1)

    this%nLegendre = ( 3*(this%maxj+1)/2+1 ) / 2 + 1
    if ( mod(this%nLegendre,2) /= 0 ) this%nLegendre = this%nLegendre+1
    
    this%scale = 1 / ( 8 * this%nLegendre**2 * sqrt(pi) )

    allocate( this%roots(this%nLegendre), this%fftLege(this%nLegendre), this%amjrr(this%jms2), this%bmjrr(this%jms2) )
    
    n = this%nLegendre
      do
        n = 2*n; xincr = 1._dbl/n
        if (xincr < 1.0d-15) exit

        ncnt = 0
        x = 0._dbl; fx = lege_fn(2*this%nLegendre, x)
        do i = 1, n
          y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
          if (fx*fy < 0._dbl) ncnt = ncnt+1
          x = y; fx = fy
        end do

        if (ncnt == this%nLegendre) exit
      end do

    i = 0
    x = 0._dbl   ; fx = lege_fn(2*this%nLegendre, x)
    y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
      do
        if (fx*fy < 0._dbl) then
          i = i+1
            this%roots(i) = xnode_fn(this%nLegendre, x, y, fx, fy)
            if (i == this%nLegendre) exit
        end if

        x = y;         fx = fy
        y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
      end do
    
    do m = 0, this%maxj
      do j = m+1, this%maxj
        this%amjrr(m*(this%maxj+1)-m*(m+1)/2+j+1) = sqrt((2*j-1._dbl)*(2*j+1._dbl)                          /(        (j-m)*(j+m)))
        this%bmjrr(m*(this%maxj+1)-m*(m+1)/2+j+1) = sqrt(             (2*j+1._dbl)*(j-m-1._dbl)*(j+m-1._dbl)/((2*j-3)*(j-m)*(j+m)))
      end do
    end do
      
    do i = 1, this%nLegendre
      this%fftLege(i) = (1-this%roots(i)**2) / lege_fn(2*this%nLegendre-1, this%roots(i))**2
    end do
    
    call this%fourtrans%init_sub( this%nFourier )
    
    this%tolm = 0.1_dbl; call this%vctol_sub();
    
  end subroutine init_harmonics_sub
  
  module pure subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call this%fourtrans%deallocate_sub()
    
    if ( allocated(this%roots)   ) deallocate( this%roots   )
    if ( allocated(this%fftLege) ) deallocate( this%fftLege )
    if ( allocated(this%amjrr)   ) deallocate( this%amjrr   )
    if ( allocated(this%bmjrr)   ) deallocate( this%bmjrr   )
    
  end subroutine deallocate_harmonics_sub
  
end submodule Init_SphericalHarmonics
