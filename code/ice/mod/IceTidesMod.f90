module IceTidesMod
  use IceMod
  implicit none
    
  type, extends(T_ice), public :: T_iceTides
    
    contains
    
    procedure :: init_sub       => init_iceTides_sub
    procedure :: compute_sub    => compute_iceTides_sub
    procedure :: visc_r_fn      => visc_iceTides_fn
    procedure :: Vdelta_fn      => Vdelta_iceTides_fn
    procedure :: set_layers_sub => set_layers_iceTides_sub
    
    procedure, private :: EE_mech_sub => EE_mech_iceTides_sub
    
  end type T_iceTides
  
  private :: init_iceTides_sub, compute_iceTides_sub, EE_mech_iceTides_sub, Vdelta_iceTides_fn, set_layers_iceTides_sub, &
           & visc_iceTides_fn
  
  contains
  
  subroutine init_iceTides_sub(this, latvisc)
    class(T_iceTides), intent(inout) :: this
    logical,           intent(in)    :: latvisc
    
    call this%init_ice_sub(jmax_in=2, rheol_in='viscel', n_iter=n_iter_tides, noharm=.true.)
      this%cf = one
    
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    
    if ( latvisc ) then
      call this%mparams%init_visc_sub()
    else
      call this%mparams%init_visc_radial_sub()
    end if
    
  end subroutine init_iceTides_sub
    
  subroutine compute_iceTides_sub(this, visc_prof)
    class(T_iceTides), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: visc_prof(this%nd,*)
    integer                          :: ijm, ir, n
    real(kind=dbl)                   :: P, Pglobal
    
    call this%sol%nulify_sub() ; Pglobal = zero
    
    if ( this%mparams%initvisc ) then
      do concurrent ( ijm = 1:this%jms, ir = 1:this%nd )
        this%mparams%visc(ir,ijm) = visc_prof(ir,ijm)
      end do
    else
      do concurrent ( ir = 1:this%nd )
        this%mparams%visc_radial(ir) = visc_prof(ir,1)
      end do
    end if
    
    this%t = zero
    this%dt = this%period / this%n_iter ; this%htide = czero
    
    do
      do n = 1, this%n_iter
        this%t = this%t + this%dt
        
        call this%EE_mech_sub()
        call this%tidal_heating_sub()
      end do
          
      P = this%rad_grid%intV_fn( real(this%htide(:,1), kind=dbl) )
        if ( abs(P-Pglobal) / P < 1.0d-3 ) then
          exit
        else
          Pglobal    = P
          this%htide = czero
        end if
    end do
    
  end subroutine compute_iceTides_sub
    
    subroutine EE_mech_iceTides_sub(this)
      class(T_iceTides), intent(inout) :: this
      integer                          :: ijm
      
      do concurrent ( ijm = 4:6:2 )
        this%rsph1(          1,ijm) = -this%sol%u_dn(ijm) + this%Vdelta_fn(1,ijm)
        this%rsph1(2:this%nd+1,ijm) = czero
        
        this%rsph2(1:this%nd,ijm) = czero
        this%rsph2(this%nd+1,ijm) = -this%sol%u_up(ijm) + this%Vdelta_fn(this%nd,ijm)
      end do
      
      call this%solve_mech_sub( ijmstart=4, ijmend=6, ijmstep=2, rematrix=.true., matxsol=.false. )
      
      do concurrent ( ijm = 4:6:2 )
        this%sol%u_dn(ijm) = this%sol%u_dn(ijm) + this%vr_r_fn(1      ,ijm) * this%dt
        this%sol%u_up(ijm) = this%sol%u_up(ijm) + this%vr_r_fn(this%nd,ijm) * this%dt
      end do
      
    end subroutine EE_mech_iceTides_sub
    
    pure complex(kind=dbl) function Vdelta_iceTides_fn(this, ir, ijm)
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
    
    subroutine set_layers_iceTides_sub(this)
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
    
    subroutine andrade_visc_iceTides_sub(this, visc_prof_jm)
      class(T_iceTides), intent(in) :: this
      complex(kind=dbl), intent(in) :: visc_prof_jm(:,:)
      
      !***********************************!
      ! Treba napisat prechod na grid pre !
      ! slapovu viskozitu.                !
      !***********************************!
      
    end subroutine andrade_visc_iceTides_sub
    
    pure real(kind=dbl) function visc_iceTides_fn(this, ir)
      class(T_iceTides), intent(in) :: this
      integer,           intent(in) :: ir
      
      if ( this%mparams%initvisc ) then
        visc_iceTides_fn = c2r_fn( this%mparams%visc(1,ir) ) / s4pi
      else
        visc_iceTides_fn = c2r_fn( this%visc_radial(ir) )
      end if
      
    end function visc_iceTides_fn
  
end module IceTidesMod