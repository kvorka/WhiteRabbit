module IceTidesMod
  use Math
  use IceMod
  use MatrixDefinitions
  implicit none
    
  type, extends(T_ice), public :: T_iceTides
    real(kind=dbl)                                  , private :: Pglobal
    complex(kind=dbl), dimension(:,:,:), allocatable, private :: Rtide
    
    contains
    
    procedure, public  :: init_sub       => init_iceTides_sub
    procedure, public  :: deallocate_sub => deallocate_iceTides_sub
    
  end type T_iceTides
  
  private :: init_iceTides_sub
  private :: EE_mech_iceTides_sub
  private :: EE_temp_iceTides_sub
  private :: Vdelta_iceTides_fn
  private :: set_layers_iceTides_sub
  private :: vypis_iceTides_sub
  private :: deallocate_iceTides_sub
  
  contains
  
  subroutine init_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    logical                          :: exitCode
    integer                          :: n, i
    real(kind=dbl)                   :: P
    
    !Vseobecna mobilizacia a nastavenie pociatocneho stavu
    call this%init_ice_sub(jmax_in = 2, rheol_in = 'viscel', n_iter = n_iter_tides, noharm = .true.)

    allocate( this%Rtide(3,this%nd,this%jms) )

    do i = 1, this%nd+1
      this%sol%temp(3*(i-1)+1,1) = cmplx( 0.5_dbl, 0._dbl, kind=dbl)
    end do
    
    !Vypocet
    do
      !Ustalenie teploty
      call EE_temp_iceTides_sub(this, exitCode); if ( exitCode ) exit
      
      !Mechanicka cast riesenia
      this%Rtide    = cmplx(0._dbl, 0._dbl, kind=dbl) ; this%sol%mech = cmplx(0._dbl, 0._dbl, kind=dbl)
      this%sol%u_C  = cmplx(0._dbl, 0._dbl, kind=dbl) ; this%sol%u_I2 = cmplx(0._dbl, 0._dbl, kind=dbl)
      this%sol%u_dn = cmplx(0._dbl, 0._dbl, kind=dbl) ; this%sol%u_up = cmplx(0._dbl, 0._dbl, kind=dbl)
      this%htide    = cmplx(0._dbl, 0._dbl, kind=dbl) ; this%t = 0._dbl
      
      do
        do n = 1, this%n_iter
          call EE_mech_iceTides_sub(this)
        end do
        
        P = this%rad_grid%intV_fn( real(this%htide(:,1), kind=dbl) )
        
        if ( abs(P-this%Pglobal) / P < 1.0d-6 ) then
          exit
        else
          this%Pglobal = P
          this%htide   = cmplx(0._dbl, 0._dbl, kind=dbl)
        end if
      end do
    end do
    
    call vypis_iceTides_sub(this)
    
  end subroutine init_iceTides_sub
  
    subroutine vypis_iceTides_sub(this)
      class(T_iceTides), intent(in) :: this
      integer                       :: i
    
      open(unit=7, file='data/data_ice_tides/Temp_tides.dat', status='new', action='write')
        do i = 1, this%nd+1
          write(7,*) this%rad_grid%rr(i), this%sol%temp_fn(i,0,0)
        end do
      close(7)
    
      open(unit=7, file='data/data_ice_tides/tidal_heating.dat', status='new', action='write')
        do i = 1, this%nd
          write(7,*) this%rad_grid%r(i), this%htide(i,:)
        end do
      close(7)
    
      write(*,*) 'Tides initialized'
    
    end subroutine vypis_iceTides_sub

    subroutine EE_mech_iceTides_sub(this)
      class(T_iceTides),             intent(inout) :: this
      integer                                      :: i, j, m, jm1
      complex(kind=dbl), dimension(:), allocatable :: Dstrss, H
    
      !Matice a casovy krok
      this%dt = this%period / this%n_iter
      call this%mat%mech(2)%fill_sub( matica_mech_fn(this, j_in=2, a_in=1._dbl), &
                                    & matica_mech_fn(this, j_in=2, a_in=0._dbl)  )

      !Selfgravitacia  
      !call set_layers_iceTides_sub(this)
    
      !Priprava pravych stran
      do jm1 = 4, 6, 2
        i = 1
          this%Rtide(1,i,jm1) = this%Rtide(1,i,jm1) - this%Ramu * this%dt / this%visc_fn(i) * this%sol%mech(6*(i-1)+3,jm1)
          this%Rtide(2,i,jm1) = this%Rtide(2,i,jm1) - this%Ramu * this%dt / this%visc_fn(i) * this%sol%mech(6*(i-1)+5,jm1)
          this%Rtide(3,i,jm1) = this%Rtide(3,i,jm1) - this%Ramu * this%dt / this%visc_fn(i) * this%sol%mech(6*(i-1)+6,jm1)
        
          this%sol%mech( 6*(i-1)+1             , jm1 ) = Vdelta_iceTides_fn(this, 2, jm1-4, this%rd)
          this%sol%mech( 6*(i-1)+2 : 6*(i-1)+3 , jm1 ) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%mech( 6*(i-1)+4 : 6*(i-1)+6 , jm1 ) = this%Rtide(1:3,i,jm1)
      
        do i = 2, this%nd
          this%Rtide(1,i,jm1) = this%Rtide(1,i,jm1) - this%Ramu * this%dt / this%visc_fn(i) * this%sol%mech(6*(i-1)+3,jm1)
          this%Rtide(2,i,jm1) = this%Rtide(2,i,jm1) - this%Ramu * this%dt / this%visc_fn(i) * this%sol%mech(6*(i-1)+5,jm1)
          this%Rtide(3,i,jm1) = this%Rtide(3,i,jm1) - this%Ramu * this%dt / this%visc_fn(i) * this%sol%mech(6*(i-1)+6,jm1)
        
          this%sol%mech( 6*(i-1)+1 : 6*(i-1)+3 , jm1 ) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%mech( 6*(i-1)+4 : 6*(i-1)+6 , jm1 ) = this%Rtide(1:3,i,jm1)
        end do
      
        i = this%nd+1
          this%sol%mech(6*(i-1)+1,jm1) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%mech(6*(i-1)+2,jm1) = Vdelta_iceTides_fn(this, 2, jm1-4, this%ru)
        
        call this%mat%mech(2)%luSolve_sub(this%sol%mech(:,jm1))
      
        this%sol%u_dn(jm1) = this%vr_fn(       1 , 2 , jm1-4 )   
        this%sol%u_up(jm1) = this%vr_fn( this%nd , 2 , jm1-4 )
      end do
    
      !Slapove zahrievanie
      allocate( Dstrss(jml2(this%jmax,this%jmax,this%jmax+2)), H(jms4) )
        do i = 1, this%nd
          H = cmplx(0._dbl, 0._dbl, kind=dbl)
          Dstrss = this%sol%deviatoric_stress_jml2_fn(i)
        
          H(1) =     Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0,-2))  / sqrt(4*pi) + &
               & 2 * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2,-2))) / sqrt(4*pi) + &
               &     Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  / sqrt(4*pi) + &
               & 2 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) / sqrt(4*pi) + &
               &     Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  / sqrt(4*pi) + &
               & 2 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2))) / sqrt(4*pi)
        
          H(1) = cmplx(real( H(1), kind=dbl ), 0._dbl, kind=dbl)
        
          H(4) = - 2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0, 0))  + &
               &   2 / sqrt(14*pi)        * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2, 0))) + &
               &   2 / sqrt(14*pi)        * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,-2))) - &
               &   3 * sqrt(5/pi/49) / 14 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  + &
               &   3 * sqrt(5/pi/49) /  7 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) - &
               &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0,+2))  + &
               &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,+2))) + &
               &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2, 0))) + &
               &   5 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  - &
               &  10 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2)))
        
          H(4) = cmplx(real( H(4), kind=dbl ), 0._dbl, kind=dbl)
                      
          H(6) =  2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) * Dstrss(jml2(2,2, 0)) + &
               &  2 / sqrt(14*pi)        * Dstrss(jml2(2,2,-2)) * Dstrss(jml2(2,0, 0)) + &
               &  3 * sqrt(5/pi/49) /  7 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2, 0)) + &
               & 12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2,+2)) + &
               & 12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,0,+2)) - &
               & 10 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,0,+2)) * Dstrss(jml2(2,2,+2))
        
          H(11) =  2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0,+2))  + &
                &  1 / sqrt(14*pi)   /  3 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,-2))) + &
                &  1 / sqrt(14*pi)   /  3 * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2,+2))) + &
                &  6 / sqrt(   pi)   / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  + &
                &  2 / sqrt(   pi)   / 49 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) - &
                & 50 / sqrt( 5*pi)   / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0,+2))  - &
                & 25 / sqrt( 5*pi)   /147 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,+2))) - &
                & 25 / sqrt( 5*pi)   /147 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2, 0))) + &
                &  9 / sqrt(   pi)   / 98 * Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  + &
                &  6 / sqrt(   pi)   /196 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2)))
        
          H(11) = cmplx(real( H(11), kind=dbl ), 0._dbl, kind=dbl)
        
          H(13) =      sqrt( 5/pi/42)      * Dstrss(jml2(2,0,-2)) * Dstrss(jml2(2,2,+2)) + &
                &      sqrt( 5/pi/42)      * Dstrss(jml2(2,2,-2)) * Dstrss(jml2(2,0,+2)) + &
                &  2 * sqrt(15/pi)    / 49 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2, 0)) - &
                & 25 * sqrt( 3/pi)    /147 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2,+2)) - &
                & 25 * sqrt( 3/pi)    /147 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,0,+2)) + &
                &  3 * sqrt(15/pi)    / 98 * Dstrss(jml2(2,0,+2)) * Dstrss(jml2(2,2,+2))
        
          H(15) =      sqrt( 5/pi)    /  3 * Dstrss(jml2(2,2,-2)) * Dstrss(jml2(2,2,+2)) + &
                &      sqrt(10/pi/7)  /  7 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,2, 0)) - &
                & 50 / sqrt(14*pi)    / 21 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,2,+2)) + &
                &  3 * sqrt( 5/pi/14) / 14 * Dstrss(jml2(2,2,+2)) * Dstrss(jml2(2,2,+2))
        
          this%htide(i,:) = this%htide(i,:) + H(:) / this%visc_fn(i) / 2 / this%n_iter
        end do    
      deallocate( Dstrss, H )

      !Cas
      this%t = this%t + this%dt
      
    end subroutine EE_mech_iceTides_sub
  
    subroutine EE_temp_iceTides_sub(this, exitCode)
      class(T_iceTides), intent(inout) :: this
      logical,           intent(out)   :: exitCode
      integer                          :: i
      complex(kind=dbl), allocatable   :: Temp(:)
    
      this%dt = huge(0._dbl)
      call this%mat%temp(0)%fill_sub( matica_temp_fn(this, j_in=0, a_in=1._dbl), &
                                      matica_temp_fn(this, j_in=0, a_in=0._dbl)  )
      
      allocate( Temp(this%nd+1) ); Temp = this%sol%temp_i_fn(0,0)
        i = 1
          this%sol%temp(3*(i-1)+1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
          this%sol%temp(3*(i-1)+2,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%temp(3*(i-1)+3,1) = cmplx(0._dbl, 0._dbl, kind=dbl)

        do i = 2, this%nd
          this%sol%temp(3*(i-1)+1,1) = this%Ds/this%Ra * this%htide_fn(i,1) / this%cp_fn(i)
          this%sol%temp(3*(i-1)+2,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%temp(3*(i-1)+3,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
        end do

        i = this%nd+1
          this%sol%temp(3*(i-1)+1,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
    
        call this%mat%temp(0)%luSolve_sub(this%sol%temp(:,1))

        if ( maxval(abs(this%sol%temp_i_fn(0,0)-Temp)/abs(Temp)) < 1e-5 ) then
          exitCode = .true.
        else
          exitCode = .false.
        end if

      deallocate( Temp )

    end subroutine EE_temp_iceTides_sub

      pure complex(kind=dbl) function Vdelta_iceTides_fn(this, j, m, ri)
        class(T_iceTides), intent(in) :: this
        integer,           intent(in) :: j, m
        real(kind=dbl),    intent(in) :: ri

        Vdelta_iceTides_fn = this%gravity%V_tide_fn(j, m, ri, 2*pi*this%t/this%period)                                 + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW -this%rhoI , this%sol%u_dn(jm(j,m))) + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI            , this%sol%u_up(jm(j,m))) + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW , this%sol%u_I2(jm(j,m))) + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC -this%rhoI2, this%sol%u_C(jm(j,m)) )
        
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
  
  subroutine deallocate_iceTides_sub(this)
    class(T_iceTides), intent(inout) :: this
    
    deallocate( this%Rtide )
    call this%deallocate_ice_sub()
    
  end subroutine deallocate_iceTides_sub
  
end module IceTidesMod