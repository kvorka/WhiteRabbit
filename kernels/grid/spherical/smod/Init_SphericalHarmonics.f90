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
    
    this%jmax  = jmax
    this%jmax1 = jmax+1
    this%jmax2 = jmax+2
    this%jmax3 = jmax+3

    this%jms  =     ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jms1 =     ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    this%jms2 =     ((jmax+2)*(jmax+3)/2 + (jmax+2)) + 1

    this%jmv  = 3 * ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jmv1 = 3 * ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    
    this%nFourier  = 3*(this%jmax2+1)
    this%nLegendre = ( 3*(this%jmax2+1)/2+1 ) / 2 + 1
    if ( mod(this%nLegendre,4) /= 0 ) this%nLegendre = this%nLegendre+4-mod(this%nLegendre,4)

    allocate( this%cosx(this%nLegendre) )
    
    n = this%nLegendre
      do
        n = 2*n; xincr = one/n
        if (xincr < 1.0d-15) exit

        ncnt = 0
        x = zero; fx = lege_fn(2*this%nLegendre, x)
        do i = 1, n
          y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
          if (fx*fy < zero) ncnt = ncnt+1
          x = y; fx = fy
        end do

        if (ncnt == this%nLegendre) exit
      end do

    i = 0
    x = zero     ; fx = lege_fn(2*this%nLegendre, x)
    y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
      do
        if (fx*fy < zero) then
          i = i+1
            this%cosx(i) = xnode_fn(this%nLegendre, x, y, fx, fy)
            if (i == this%nLegendre) exit
        end if

        x = y;         fx = fy
        y = x + xincr; fy = lege_fn(2*this%nLegendre, y)
      end do
    
    allocate( this%amj(this%jms2), this%bmj(this%jms2), this%cmm(this%jms2) )
    
    do m = 0, this%jmax2
      if ( m == 0 ) then
        this%cmm(m) = one / s4pi
      else
        this%cmm(m) = -sqrt( (2*m+one) / (2*m) )
      end if
      
      do j = m+1, this%jmax2
        this%amj(m*this%jmax3-m*(m+1)/2+j+1) = sqrt((2*j-1)*(2*j+one)                /(        (j-m)*(j+m)))
        this%bmj(m*this%jmax3-m*(m+1)/2+j+1) = sqrt(        (2*j+one)*(j-m-1)*(j+m-1)/((2*j-3)*(j-m)*(j+m)))
      end do
    end do
    
    allocate( this%weight(this%nLegendre) )
    
    do i = 1, this%nLegendre
      this%weight(i) = pi * (1-this%cosx(i)**2) / ( this%nLegendre * lege_fn(2*this%nLegendre-1, this%cosx(i)) )**2
    end do
    
    call this%fourtrans%init_sub( this%nFourier )
    
  end subroutine init_harmonics_sub
  
  module pure subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call this%fourtrans%deallocate_sub()
    
    if ( allocated(this%cosx)   ) deallocate( this%cosx   )
    if ( allocated(this%weight) ) deallocate( this%weight )
    if ( allocated(this%amj)    ) deallocate( this%amj   )
    if ( allocated(this%bmj)    ) deallocate( this%bmj   )
    if ( allocated(this%cmm)    ) deallocate( this%cmm     )
    
  end subroutine deallocate_harmonics_sub
  
end submodule Init_SphericalHarmonics
