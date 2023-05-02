submodule (SphericalHarmonics) Init_SphericalHarmonics
  implicit none
  
  contains
  
  subroutine init_harmonics_sub(this, jmax)
    class(T_lateralGrid), intent(inout) :: this
    integer,              intent(in)    :: jmax
    integer                             :: i, k, j, m, n, ncnt
    real(kind=dbl)                      :: xincr, x, y, fx, fy

    this%jmax = jmax
    this%jms  =     ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jms1 =     ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    this%jms2 =     ((jmax+2)*(jmax+3)/2 + (jmax+2)) + 1
    this%jmv  = 3 * ( jmax   *(jmax+1)/2 +  jmax   ) + 1
    this%jmv1 = 3 * ((jmax+1)*(jmax+2)/2 + (jmax+1)) + 1
    
    this%maxj      = jmax+2
    this%nFourier  = 3*(this%maxj+1)
    this%nLegendre = (((3*(this%maxj+1)/2+1)/2+1+step)/step)*step

    allocate( this%ish(this%jms2), this%roots(this%nLegendre), this%fftLege(this%nLegendre) )

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
        this%ish(m*(this%maxj+1)-m*(m+1)/2+j+1) = sqrt((j**2-m**2) / (4*j**2-1._dbl))
      end do
    end do
      
    do i = 1, this%nLegendre
      this%fftLege(i) = (1-this%roots(i)**2) / lege_fn(2*this%nLegendre-1, this%roots(i))**2
    end do
      
  end subroutine init_harmonics_sub
  
  subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call destroy_plan_sub(this%fftw_01_r2c)
    call destroy_plan_sub(this%fftw_03_r2c)
    call destroy_plan_sub(this%fftw_04_r2c)
    
    call destroy_plan_sub(this%fftw_06_c2r)
    call destroy_plan_sub(this%fftw_16_c2r)
    call destroy_plan_sub(this%fftw_19_c2r)
    
    deallocate(this%roots, this%fftLege, this%ish)
    
  end subroutine deallocate_harmonics_sub
  
  subroutine destroy_plan_sub(plan)
    type(C_ptr), intent(inout) :: plan
    
    if ( c_associated(plan) ) call fftw_destroy_plan(plan)
    
  end subroutine destroy_plan_sub
  
end submodule Init_SphericalHarmonics