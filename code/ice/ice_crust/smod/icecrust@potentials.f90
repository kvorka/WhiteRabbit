submodule (icecrust) potentials
  implicit none; contains
  
  module pure complex(kind=dbl) function Vdelta_iceCrust_fn(this, ir, ijm)
    class(T_iceCrust),  intent(in) :: this
    integer,            intent(in) :: ir, ijm
    integer                        :: k, j, m
    real(kind=dbl)                 :: ri, fac
    complex(kind=dbl), allocatable :: field(:)
    
    j   = this%j_indx(ijm)
    m   = ijm - (j*(j+1)/2+1)
    ri  = this%rad_grid%r(ir)
    fac = -this%rhoI * this%alphaU * (this%Td-this%Tu)
    
    allocate( field(this%nd+1) ); field = czero
    
    if (ir == 1) then
      do concurrent ( k = 1:this%nd+1 )
        field(k) = fac * this%alpha_rr_fn(k) * this%temp_rr_fn(k,ijm) * ( this%rd / this%rad_grid%rr(k) ) ** (j-1)
      end do
    else if (ir == this%nd) then
      do concurrent ( k = 1:this%nd+1 )
        field(k) = fac * this%alpha_rr_fn(k) * this%temp_rr_fn(k,ijm) * ( this%rad_grid%rr(k) / this%ru ) ** (j+2)
      end do
    end if
    
    Vdelta_iceCrust_fn = ( this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI           , this%bnd%u_up(ijm) ) + &
                         & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW-this%rhoI , this%bnd%u_dn(ijm) ) + &
                         & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW, this%bnd%u_I2(ijm) ) + &
                         & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC-this%rhoI2, this%bnd%u_C(ijm)  ) + &
                         & this%gravity%V_rho_fn( j, m, ri, field, this%rad_grid)                                + &
                         & this%gravity%V_rt_fn(  j, m, ri ) ) / this%gravity%g_fn( ri )
    
    deallocate( field )
    
    if ( m == 0 ) Vdelta_iceCrust_fn%im = zero
    
  end function Vdelta_iceCrust_fn
  
  module subroutine set_layers_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
    integer                          :: ir, j, m, ijm
    real(kind=dbl)                   :: a11, a12, a21, a22, det, fac
    complex(kind=dbl)                :: rhs1, rhs2
    complex(kind=dbl), allocatable   :: field(:)
    
    fac = -this%rhoI * this%alphaU * (this%Td-this%Tu) ; rhs1 = czero ; rhs2 = czero
    
    allocate( field(this%nd+1) ) ; field = czero
    
    do j = 1, this%jmax
      a11 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g - this%gravity%g_fn(this%rI2)
      a12 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g * (this%rC/this%rI2)**(j+2)
      a21 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g * (this%rC/this%rI2)**(j-1)
      a22 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g - this%gravity%g_fn(this%rC)
      
      det = a11 * a22 - a12 * a21
        a11 = a11 / det; a12 = a12 / det
        a21 = a21 / det; a22 = a22 / det
      
      do m = 0, j
        ijm = j*(j+1)/2+m+1
        
        do concurrent ( ir = 1:this%nd+1 )
          field(ir) = fac * this%alpha_rr_fn(ir) * this%temp_rr_fn(ir,ijm) * (this%rI2/this%rad_grid%rr(ir))**(j-1)
        end do
        
        rhs1 = -( this%gravity%V_bnd_fn(j,m,this%rI2,this%ru,this%rhoI          ,this%bnd%u_up(ijm)+this%bnd%v_up(ijm)*this%dt ) + &
                & this%gravity%V_bnd_fn(j,m,this%rI2,this%rd,this%rhoW-this%rhoI,this%bnd%u_dn(ijm)+this%bnd%v_dn(ijm)*this%dt ) + &
                & this%gravity%V_rho_fn(j,m,this%rI2, field, this%rad_grid) + &
                & this%gravity%V_rt_fn( j,m,this%rI2) )
        
        do concurrent ( ir = 1:this%nd+1 )
          field(ir) = fac * this%alpha_rr_fn(ir) * this%temp_rr_fn(ir,ijm) * (this%rC/this%rad_grid%rr(ir))**(j-1)
        end do
        
        rhs2 = -( this%gravity%V_bnd_fn(j,m,this%rC,this%ru,this%rhoI          ,this%bnd%u_up(ijm)+this%bnd%v_up(ijm)*this%dt ) + &
                & this%gravity%V_bnd_fn(j,m,this%rC,this%rd,this%rhoW-this%rhoI,this%bnd%u_dn(ijm)+this%bnd%v_dn(ijm)*this%dt ) + &
                & this%gravity%V_rho_fn(j,m,this%rC,field,this%rad_grid) + &
                & this%gravity%V_rt_fn( j,m,this%rC) )
        
        this%bnd%u_I2(ijm) = a22 * rhs1 - a12 * rhs2
        this%bnd%u_C(ijm)  = a11 * rhs2 - a21 * rhs1  
      end do
    end do
    
    deallocate(field)
    
  end subroutine set_layers_iceCrust_sub
  
end submodule potentials