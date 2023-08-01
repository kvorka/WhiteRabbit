submodule(IceMod) GravDeformIce
  implicit none
  
  contains

  pure complex(kind=dbl) function Vdelta_ice_fn(this, ir, ijm)
    class(T_ice),       intent(in) :: this
    integer,            intent(in) :: ir, ijm
    integer                        :: k, j, m
    real(kind=dbl)                 :: ri
    complex(kind=dbl), allocatable :: field(:)
    
    j  = this%j_indx(ijm)
    m  = ijm - (j*(j+1)/2+1)
    ri = this%rad_grid%r(ir)
    
    allocate( field(this%nd+1) ); field = czero
    
    do k = 1, this%nd+1
      field(k) = -this%rhoI * this%alphaU * this%alpha_fn(k) * (this%Td-this%Tu) * this%sol%temp_fn(k,ijm)
    end do
    
    if (ir == 1) then
      field = field * ( this%rad_grid%r(1) / this%rad_grid%rr(:) ) ** (j-1)
    else if (ir == this%nd) then
      field = field * ( this%rad_grid%rr(:) / this%rad_grid%r(this%nd) ) ** (j+2)
    end if
    
    Vdelta_ice_fn = ( this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI           , this%sol%u_up(ijm) ) + &
                    & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW-this%rhoI , this%sol%u_dn(ijm) ) + &
                    & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW, this%sol%u_I2(ijm) ) + &
                    & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC-this%rhoI2, this%sol%u_C(ijm)  ) + &
                    & this%gravity%V_rho_fn( j, m, ri, field, this%rad_grid)                                + &
                    & this%gravity%V_rt_fn(  j, m, ri ) ) / this%gravity%g_fn( ri )
    
    deallocate( field )
    
    if ( m == 0 ) Vdelta_ice_fn%im = 0._dbl
    
  end function Vdelta_ice_fn
  
  subroutine set_layers_ice_sub(this)
    class(T_ice),      intent(inout) :: this
    integer                          :: i, j, m, jm_int
    real(kind=dbl)                   :: a11, a12, a21, a22, det
    complex(kind=dbl)                :: rhs1, rhs2
    complex(kind=dbl), allocatable   :: field(:)
    
    rhs1 = czero
    rhs2 = czero
    
    associate( rd => this%rad_grid%r(1), ru => this%rad_grid%r(this%nd), grv => this%gravity  )
    allocate( field(this%nd+1) )

    do j = 1, this%jmax
      a11 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g - grv%g_fn(this%rI2)
      a12 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g * (this%rC/this%rI2)**(j+2)
      a21 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g * (this%rC/this%rI2)**(j-1)
      a22 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g - grv%g_fn(this%rC)
      
      det = a11 * a22 - a12 * a21
        a11 = a11 / det; a12 = a12 / det
        a21 = a21 / det; a22 = a22 / det
      
      do m = 0, j
        jm_int = jm(j,m)
        
        do i = 1, this%nd+1
          field(i) = -this%rhoI * this%alphaU * this%alpha_fn(i) * &
                          & (this%Td-this%Tu) * this%sol%temp_fn(i,jm_int) * (this%rI2/this%rad_grid%rr(i))**(j-1)
        end do
        
        rhs1 = -( grv%V_bnd_fn(j,m,this%rI2,ru,this%rhoI          ,this%sol%u_up(jm_int)+this%sol%v_up(jm_int)*this%dt ) + &
                & grv%V_bnd_fn(j,m,this%rI2,rd,this%rhoW-this%rhoI,this%sol%u_dn(jm_int)+this%sol%v_dn(jm_int)*this%dt ) + &
                & grv%V_rho_fn(j,m,this%rI2, field, this%rad_grid) + &
                & grv%V_rt_fn( j,m,this%rI2) )
        
        do i = 1, this%nd+1
          field(i) = -this%rhoI * this%alphaU * this%alpha_fn(i) * &
                            & (this%Td-this%Tu) * this%sol%temp_fn(i,jm_int) * (this%rC/this%rad_grid%rr(i))**(j-1)
        end do
        
        rhs2 = -( grv%V_bnd_fn(j,m,this%rC,ru,this%rhoI          ,this%sol%u_up(jm_int)+this%sol%v_up(jm_int)*this%dt ) + &
                & grv%V_bnd_fn(j,m,this%rC,rd,this%rhoW-this%rhoI,this%sol%u_dn(jm_int)+this%sol%v_dn(jm_int)*this%dt ) + &
                & grv%V_rho_fn(j,m,this%rC,field,this%rad_grid) + &
                & grv%V_rt_fn( j,m,this%rC) )
    
        this%sol%u_I2(jm_int) = a22 * rhs1 - a12 * rhs2
        this%sol%u_C(jm_int)  = a11 * rhs2 - a21 * rhs1  
      end do
    end do

    deallocate(field)
    end associate
    
  end subroutine set_layers_ice_sub

end submodule GravDeformIce