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
    
    this%jmt  = 5 * ((jmax  )*(jmax+1)/2 + (jmax  )) + 1
    
    this%nFourier  = 3*(this%jmax2+1)
    this%nLegendre = ( 3*(this%jmax2+1)/2+1 ) / 2 + 1
    if ( mod(this%nLegendre,2) /= 0 ) this%nLegendre = this%nLegendre+1
    
    this%scale = 1 / ( 8 * this%nLegendre**2 * sqrt(pi) )

    allocate( this%roots(this%nLegendre) )
    
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
    
    allocate( this%amjrr(this%jms2), this%bmjrr(this%jms2) )
    
    do m = 0, this%jmax2
      do j = m+1, this%jmax2
        this%amjrr(m*this%jmax3-m*(m+1)/2+j+1) = sqrt((2*j-1._dbl)*(2*j+1._dbl)                          /(        (j-m)*(j+m)))
        this%bmjrr(m*this%jmax3-m*(m+1)/2+j+1) = sqrt(             (2*j+1._dbl)*(j-m-1._dbl)*(j+m-1._dbl)/((2*j-3)*(j-m)*(j+m)))
      end do
    end do
    
    allocate( this%pmm(this%nLegendre,0:this%jmax2) )
    
    do i = 1, this%nLegendre
      this%pmm(i,0) = 1._dbl
      
      do m = 1, this%jmax2
        this%pmm(i,m) = -sqrt( ( 2*m+1 ) / ( 2._dbl * m ) * (1-this%roots(i)**2) ) * this%pmm(i,m-1)
      end do
    end do
    
    allocate( this%fftLege(this%nLegendre) )
    
    do i = 1, this%nLegendre
      this%fftLege(i) = (1-this%roots(i)**2) / lege_fn(2*this%nLegendre-1, this%roots(i))**2
    end do
    
    call this%fourtrans%init_sub( this%nFourier )
    
    this%tolm = 0.1_dbl; call this%vctol_sub()
    
    allocate( this%maxm(this%nLegendre) ) ; this%maxm = this%jmax2
    
    do i = 1, (this%nLegendre/16)*16, 16
      do m = 0, this%jmax2
        if ( maxval(abs(this%pmm(i:i+15,m))) < this%tolm ) then
          this%maxm(i:i+15) = m-1
          exit
        end if
      end do
    end do
    
    do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
      do m = 0, this%jmax2
        if ( maxval(abs(this%pmm(i:i+7,m))) < this%tolm ) then
          this%maxm(i:i+7) = m-1
          exit
        end if
      end do
    end do
    
    do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
      do m = 0, this%jmax2
        if ( maxval(abs(this%pmm(i:i+3,m))) < this%tolm ) then
          this%maxm(i:i+3) = m-1
          exit
        end if
      end do
    end do
    
    do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
      do m = 0, this%jmax2
        if ( maxval(abs(this%pmm(i:i+1,m))) < this%tolm ) then
          this%maxm(i:i+1) = m-1
          exit
        end if
      end do
    end do
    
  end subroutine init_harmonics_sub
  
  module pure subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call this%fourtrans%deallocate_sub()
    
    if ( allocated(this%roots)   ) deallocate( this%roots   )
    if ( allocated(this%fftLege) ) deallocate( this%fftLege )
    if ( allocated(this%amjrr)   ) deallocate( this%amjrr   )
    if ( allocated(this%bmjrr)   ) deallocate( this%bmjrr   )
    if ( allocated(this%pmm)     ) deallocate( this%pmm     )
    
  end subroutine deallocate_harmonics_sub
  
end submodule Init_SphericalHarmonics
