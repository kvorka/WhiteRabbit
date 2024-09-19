submodule (IceTidesMod) Potentials_iceTides
  implicit none; contains
  
  module pure complex(kind=dbl) function Vdelta_iceTides_fn(this, ir, ijm)
    class(T_iceTides), intent(in) :: this
    integer,           intent(in) :: ir, ijm
    real(kind=dbl)                :: ri
    integer                       :: j, m
    
    j  = this%j_indx(ijm)
    m  = ijm - ( j*(j+1)/2 + 1 )
    ri = this%rad_grid%r(ir)
    
    Vdelta_iceTides_fn = this%gravity%V_tide_fn(j, m, ri, 2*pi*this%t/this%period)                             + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW -this%rhoI , this%sol%u_dn(ijm)) + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI            , this%sol%u_up(ijm)) + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW , this%sol%u_I2(ijm)) + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC -this%rhoI2, this%sol%u_C(ijm) )
    
    Vdelta_iceTides_fn = Vdelta_iceTides_fn / this%gravity%g_fn( ri )
    
  end function Vdelta_iceTides_fn
  
  module subroutine set_layers_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    integer                          :: i, j, m
    real(kind=dbl)                   :: a11, a12, a21, a22, det
    complex(kind=dbl)                :: rhs1, rhs2
    
    j = 2
    
    a11 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g - this%gravity%g_fn(this%rI2)
    a12 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g * (this%rC/this%rI2)**(j+2._dbl)
    a21 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g * (this%rC/this%rI2)**(j-1._dbl)
    a22 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g - this%gravity%g_fn(this%rC)
    
    det = a11 * a22 - a12 * a21
      a11 = a11 / det; a12 = a12 / det
      a21 = a21 / det; a22 = a22 / det
      
    do m = 0, j, 2
      rhs1 = -( this%gravity%V_bnd_fn( j, m, this%rI2, this%rd, this%rhoW-this%rhoI, this%sol%u_dn(jm(j,m))) + &
              & this%gravity%V_bnd_fn( j, m, this%rI2, this%ru, this%rhoI          , this%sol%u_up(jm(j,m))) + &
              & this%gravity%V_tide_fn(j, m, this%rI2, 2*pi*this%t/this%period) ) 
              
      rhs2 = -( this%gravity%V_bnd_fn( j, m, this%rC, this%rd, this%rhoW-this%rhoI, this%sol%u_dn(jm(j,m))) + &
              & this%gravity%V_bnd_fn( j, m, this%rC, this%ru, this%rhoI          , this%sol%u_up(jm(j,m))) + &
              & this%gravity%V_tide_fn(j, m, this%rC, 2*pi*this%t/this%period) )
    
      this%sol%u_I2(jm(j,m)) = a22 * rhs1 - a12 * rhs2
      this%sol%u_C(jm(j,m))  = a11 * rhs2 - a21 * rhs1
    end do
      
  end subroutine set_layers_iceTides_sub
  
end submodule Potentials_iceTides