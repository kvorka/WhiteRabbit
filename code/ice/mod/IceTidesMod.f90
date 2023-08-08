module IceTidesMod
  use IceMod
  implicit none
    
  type, extends(T_ice), public :: T_iceTides
    contains
    
    procedure, public, pass :: init_sub       => init_iceTides_sub
    procedure, public, pass :: Vdelta_fn      => Vdelta_iceTides_fn
    procedure, public, pass :: set_layers_sub => set_layers_iceTides_sub
    
  end type T_iceTides
  
  private :: init_iceTides_sub
  private :: vypis_iceTides_sub , Vdelta_iceTides_fn , set_layers_iceTides_sub
  
  contains
  
  subroutine init_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    integer                          :: n, ir, is, ijm
    real(kind=dbl)                   :: P, Pglobal
    complex(kind=dbl), allocatable   :: Temp(:)
    
    call this%init_ice_sub(jmax_in=2, rheol_in='viscel', n_iter=n_iter_tides, noharm=.true.)
      this%cf = 1._dbl ; this%andrade = .true.
    
    call this%init_eq_temp_sub( rhs=.false. , nl=.false. )
    call this%init_eq_mech_sub( rhs=.true.  , nl=.false. )
    
    allocate( Temp(this%nd+1) )
    
    Temp    = czero
    Pglobal = 0._dbl
    this%sol%temp(1:3*this%nd+1:3,1) = cone
    
    do
      this%dt = huge(0._dbl) ; Temp = this%sol%temp_i_fn(1)
      
      call this%mat%temp(0)%fill_sub( this%matica_temp_fn(j_in=0, a_in=1._dbl), &
                                    & this%matica_temp_fn(j_in=0, a_in=0._dbl)  )
      
      ir = 1
        this%sol%temp(1,1) = cs4pi
        this%sol%temp(2,1) = czero
        this%sol%temp(3,1) = czero
      
      do ir = 2, this%nd
        is = 3*(ir-1)+1
        
        this%sol%temp(is  ,1) = this%htide_fn(ir,1)
        this%sol%temp(is+1,1) = czero
        this%sol%temp(is+2,1) = czero
      end do
      
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,1) = czero
        
      call this%mat%temp(0)%luSolve_sub( this%sol%temp(:,1) ) ; this%dt = this%period / this%n_iter
      
      do
        do n = 1, this%n_iter
          do concurrent ( ijm = 4:6:2 )
            this%rsph1(          1,ijm) = -this%sol%u_dn(ijm) + this%Vdelta_fn(1,ijm)
            this%rsph1(2:this%nd+1,ijm) = czero
            
            this%rsph2(1:this%nd,ijm) = czero
            this%rsph2(this%nd+1,ijm) = -this%sol%u_up(ijm) + this%Vdelta_fn(this%nd,ijm)
          end do
          
          call this%solve_mech_sub(ijmstart=4, ijmend=6, ijmstep=2, rematrix=.true.) ; call this%tidal_heating_sub()
          
          do concurrent ( ijm = 4:6:2 )
            this%sol%u_dn(ijm) = this%sol%u_dn(ijm) + this%vr_fn(1      ,ijm) * this%dt
            this%sol%u_up(ijm) = this%sol%u_up(ijm) + this%vr_fn(this%nd,ijm) * this%dt
          end do
          
          this%t = this%t + this%dt
        end do
        
        P = this%rad_grid%intV_fn( real(this%htide(:,1), kind=dbl) ) ; write(*,*) P
        
        if ( abs(P-Pglobal) / P < 1.0d-6 ) then
          write(*,*) P ; exit
        else
          Pglobal = P ; this%htide = czero
        end if
      end do
      
      if ( maxval( abs( this%sol%temp_i_fn(1)-Temp ) / abs(Temp) ) < 1e-8 ) exit
    end do
    
    deallocate( Temp )
    
    call vypis_iceTides_sub(this)

  end subroutine init_iceTides_sub
  
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
    
end module IceTidesMod