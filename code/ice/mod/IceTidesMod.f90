module IceTidesMod
  use IceMod
  implicit none
    
  type, extends(T_ice), public :: T_iceTides
    real(kind=dbl), allocatable, private :: stress_prof(:)
    
    contains
    
    procedure, public,  pass :: init_sub       => init_iceTides_sub
    procedure, public,  pass :: compute_sub    => compute_iceTides_sub
    procedure, private, pass :: EE_temp_sub    => EE_temp_iceTides_sub
    procedure, private, pass :: EE_mech_sub    => EE_mech_iceTides_sub
    procedure, public,  pass :: visc_fn        => visc_iceTides_fn
    procedure, public,  pass :: Vdelta_fn      => Vdelta_iceTides_fn
    procedure, public,  pass :: set_layers_sub => set_layers_iceTides_sub
    
  end type T_iceTides
  
  private :: init_iceTides_sub, compute_iceTides_sub
  private :: EE_temp_iceTides_sub, EE_mech_iceTides_sub
  private :: Vdelta_iceTides_fn , set_layers_iceTides_sub , visc_iceTides_fn
  private :: vypis_iceTides_sub
  
  contains
  
  subroutine init_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    
    call this%init_ice_sub(jmax_in=2, rheol_in='viscel', n_iter=n_iter_tides, noharm=.true.)
      this%cf = 1._dbl ; this%andrade = .true.
    
    call this%init_eq_temp_sub( rhs=.true. , nl=.false. )
    call this%init_eq_mech_sub( rhs=.true. , nl=.false. )
    

  end subroutine init_iceTides_sub
    
  subroutine compute_iceTides_sub(this, stress_prof_i)
    class(T_iceTides),           intent(inout) :: this
    real(kind=dbl),    optional, intent(in)    :: stress_prof_i(:)
    integer                                    :: n
    real(kind=dbl)                             :: P, Pglobal
    complex(kind=dbl), allocatable             :: Temp(:)
    
    call this%sol%nulify_sub()
    Pglobal = 0._dbl
    
    if ( present(stress_prof_i) ) this%stress_prof = stress_prof_i(:)
    
    allocate( Temp(this%nd+1) ) ; Temp = czero ; this%sol%temp(1:3*this%nd+1:3,1) = cone
    
    do
      this%t = 0._dbl
      
      this%dt = huge(0._dbl)
        Temp = this%sol%temp_i_fn(1)
        call this%EE_temp_sub() ; if ( maxval( abs( this%sol%temp_i_fn(1)-Temp ) / abs(Temp) ) < 1e-6 ) exit
      
      this%dt = this%period / this%n_iter ; this%htide = czero
        do
          do n = 1, this%n_iter
            this%t = this%t + this%dt
              call this%EE_mech_sub()
              call this%tidal_heating_sub()
          end do
          
          P = this%rad_grid%intV_fn( real(this%htide(:,1), kind=dbl) )
            if ( abs(P-Pglobal) / P < 1.0d-6 ) then
              write(*,*) P ; exit
            else
              Pglobal    = P
              this%htide = czero
            end if
        end do
    end do
      
    deallocate( Temp )
    
    !call vypis_iceTides_sub(this)
    
  end subroutine compute_iceTides_sub
  
    subroutine EE_temp_iceTides_sub(this)
      class(T_iceTides), intent(inout) :: this
      integer                          :: ijm, ir
      complex(kind=dbl), allocatable   :: Temp1(:)
      
      allocate( Temp1(this%nd+1) )
      
      ijm = 1
        do
          Temp1(:) = this%sol%temp_i_fn(ijm)
          
          ir = 1
            this%rtemp(ir,ijm) = cs4pi
           
          do concurrent ( ir = 2:this%nd )
            this%rtemp(ir,ijm) = this%htide_fn(ir,1)
          end do
          
          ir = this%nd+1
            this%rtemp(ir,ijm) = czero
          
          call this%solve_temp_sub( ijmstart=ijm, ijmend=ijm, ijmstep=1, rematrix=.true., matxsol=.false. )
          
          if ( maxval( abs( this%sol%temp_i_fn(ijm) - Temp1 ) / abs( Temp1 ) ) < 1e-8 ) exit
        end do
        
      deallocate( Temp1 )
      
    end subroutine EE_temp_iceTides_sub
    
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
        this%sol%u_dn(ijm) = this%sol%u_dn(ijm) + this%vr_fn(1      ,ijm) * this%dt
        this%sol%u_up(ijm) = this%sol%u_up(ijm) + this%vr_fn(this%nd,ijm) * this%dt
      end do
      
    end subroutine EE_mech_iceTides_sub
    
    subroutine vypis_iceTides_sub(this)
      class(T_iceTides), intent(in) :: this
      integer                       :: i
    
      open(unit=7, file='data/data_ice_tides/Temp_tides.dat', status='new', action='write')
        do i = 1, this%nd+1
          write(7,*) this%rad_grid%rr(i), this%sol%temp_fn(i,1)
        end do
      close(7)
    
      open(unit=7, file='data/data_ice_tides/tidal_heating.dat', status='new', action='write')
        do i = 1, this%nd
          write(7,*) this%rad_grid%r(i), this%htide(i,:)
        end do
      close(7)
    
      write(*,*) 'Tides initialized'
    
    end subroutine vypis_iceTides_sub
    
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
  
    pure real(kind=dbl) function visc_iceTides_fn(this, i)
      class(T_iceTides), intent(in) :: this
      integer,           intent(in) :: i
      real(kind=dbl)                :: visc, temp
      
      temp = this%Tu + (this%Td-this%Tu) * real( this%rad_grid%c(i,-1) * this%sol%temp_fn(i  ,1) +          &
                                               & this%rad_grid%c(i,+1) * this%sol%temp_fn(i+1,1) , kind=dbl ) / sqrt(4*pi)
      
      if ( allocated(this%stress_prof) ) then
        visc = min( goldsby_visc_fn(this%diam, temp, this%stress_prof(i)), this%cutoff )
      else
        visc = min( goldsby_diffvisc_fn(this%diam, temp), this%cutoff )
      end if
      
      if ( .not. this%andrade ) then
        visc_iceTides_fn = visc
      else
        visc_iceTides_fn = andrade_visc_fn(this%mu, this%omega, visc)
      end if
      
      visc_iceTides_fn = visc_iceTides_fn / this%viscU
      
    end function visc_iceTides_fn
  
end module IceTidesMod