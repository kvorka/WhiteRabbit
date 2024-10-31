submodule (icetides) potentials
  implicit none; contains
  
  module procedure Vdelta_iceTides_fn
    real(kind=dbl) :: ri
    integer        :: j, m
    
    j  = this%j_indx(ijm)
    m  = ijm - ( j*(j+1)/2 + 1 )
    ri = this%rad_grid%r(ir)
    
    Vdelta_iceTides_fn = this%gravity%V_tide_fn(j, m, ri, 2*pi*this%t/this%period)                             + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW -this%rhoI , this%bnd%u_dn(ijm)) + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI            , this%bnd%u_up(ijm)) + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW , this%bnd%u_I2(ijm)) + &
                       & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC -this%rhoI2, this%bnd%u_C(ijm) )
    
    Vdelta_iceTides_fn = Vdelta_iceTides_fn / this%gravity%g_fn( ri )
    
  end procedure Vdelta_iceTides_fn
  
  module procedure set_layers_iceTides_sub
    integer           :: i, j, m, ijm
    real(kind=dbl)    :: a11, a12, a21, a22, det
    complex(kind=dbl) :: rhs1, rhs2
    
    j = 2
    
    a11 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g - this%gravity%g_fn(this%rI2)
    a12 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g * (this%rC/this%rI2)**(j+2._dbl)
    a21 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g * (this%rC/this%rI2)**(j-1._dbl)
    a22 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g - this%gravity%g_fn(this%rC)
    
    det = a11 * a22 - a12 * a21
      a11 = a11 / det; a12 = a12 / det
      a21 = a21 / det; a22 = a22 / det
      
    do m = 0, j, 2
      ijm = j*(j+1)/2+m+1
      
      rhs1 = -( this%gravity%V_bnd_fn( j, m, this%rI2, this%rd, this%rhoW-this%rhoI, this%bnd%u_dn(ijm)) + &
              & this%gravity%V_bnd_fn( j, m, this%rI2, this%ru, this%rhoI          , this%bnd%u_up(ijm)) + &
              & this%gravity%V_tide_fn(j, m, this%rI2, 2*pi*this%t/this%period) ) 
              
      rhs2 = -( this%gravity%V_bnd_fn( j, m, this%rC, this%rd, this%rhoW-this%rhoI, this%bnd%u_dn(ijm)) + &
              & this%gravity%V_bnd_fn( j, m, this%rC, this%ru, this%rhoI          , this%bnd%u_up(ijm)) + &
              & this%gravity%V_tide_fn(j, m, this%rC, 2*pi*this%t/this%period) )
    
      this%bnd%u_I2(ijm) = a22 * rhs1 - a12 * rhs2
      this%bnd%u_C(ijm)  = a11 * rhs2 - a21 * rhs1
    end do
      
  end procedure set_layers_iceTides_sub
  
end submodule potentials