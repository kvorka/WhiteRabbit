module IceCrustMod
  use Math
  use IceMod
  use IceTidesMod
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    
    contains
    
    procedure, public,  pass :: init_sub        => init_iceCrust_sub
    procedure, public,  pass :: time_scheme_sub => iter_iceCrust_sub
    procedure, public,  pass :: htide_fn        => htide_iceCrust_fn
    procedure, private, pass :: Vdelta_fn       => Vdelta_iceCrust_fn
    procedure, private, pass :: setLayers_sub   => set_layers_iceCrust_sub
    procedure, public,  pass :: deallocate_sub  => deallocate_iceCrust_sub
    
  end type T_iceCrust
  
  private :: init_iceCrust_sub
  private :: find_hydrostatic_iceCrust_sub
  private :: EE_iceCrust_sub
  private :: Vdelta_iceCrust_fn
  private :: set_layers_iceCrust_sub
  private :: htide_iceCrust_fn
  private :: vypis_iceCrust_sub
  private :: deallocate_iceCrust_sub
  
  contains
  
  subroutine init_iceCrust_sub(this, notides)
    class(T_iceCrust), intent(inout) :: this
    logical, optional, intent(in)    :: notides
    type(T_iceTides)                 :: tides
    integer                          :: i, j, m
    real(kind=dbl)                   :: radius
   
    call this%init_ice_sub(jmax_in = jmax_ice, rheol_in = 'viscos', n_iter = n_iter_ice)
    call this%lat_grid%init_vcvv_sub()

    allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp = czero
    
    call tides%init_sub()
    call tides%deallocate_sub()

    if (.not. notides) then
      open(unit=1, file='data/data_ice_tides/Temp_tides.dat', status='old', action='read')
        do i = 1, this%nd+1; read(1,*) radius, this%sol%temp(3*(i-1)+1,1); end do
      close(1)
      
    else
      open(unit=1, file='data/data_ice_tides/Temp_tides.dat', status='old', action='read')
        do i = 1, this%nd+1; read(1,*) radius, this%sol%temp(3*(i-1)+1,1); end do
      close(1)
    
      open(unit=1, file='data/data_ice_tides/tidal_heating.dat', status='old', action='read')
        do i = 1, this%nd; read(1,*) radius, this%htide(i,:); end do
      close(1)
    end if

    call find_hydrostatic_iceCrust_sub(this)
    call vypis_iceCrust_sub(this)

  end subroutine init_iceCrust_sub

    subroutine find_hydrostatic_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
      complex(kind=dbl), allocatable   :: flux_bnd(:)
      
      allocate( flux_bnd(this%jms) ); flux_bnd = czero
      
      do
        call EE_iceCrust_sub(this, flux_bnd)
        write(*,*) realpart(this%sol%u_dn(4)) * this%D_ud , realpart(this%sol%u_up(4)) * this%D_ud
        if ( abs(this%sol%v_up(4) * this%dt / this%sol%u_up(4)) < 1e-5 ) exit
      end do
      
      deallocate( flux_bnd )
      
      write(*,*) 'Icy crust quasi hydrostatic'
      
    end subroutine find_hydrostatic_iceCrust_sub
  
  subroutine iter_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux_bnd(:)
    integer                          :: n
    
    do n = 1, this%n_iter
      call EE_iceCrust_sub(this, flux_bnd)
    end do

    call vypis_iceCrust_sub(this)
    
  end subroutine iter_iceCrust_sub
  
  subroutine EE_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux_bnd(:)
    integer                          :: ir, is, ij, ijm
    real(kind=dbl)                   :: qConv
    complex(kind=dbl), allocatable   :: Temp(:), Temp1(:)
    
    this%t = this%t + this%dt
    
    !$omp parallel do
    do ir = 2, this%nd
      this%ntemp(:,ir) = -this%vgradT_fn(ir)
    end do
    !$omp end parallel do
    
    allocate( Temp(this%nd+1), Temp1(this%nd+1) ); Temp = this%sol%temp_i_fn(0,0)
      do
        Temp1 = this%sol%temp_i_fn(0,0)
        call this%mat%temp(0)%fill_sub( this%matica_temp_fn(j_in=0, a_in=1._dbl), this%matica_temp_fn(j_in=0, a_in=0._dbl)  )
      
        ir = 1
          this%sol%temp(1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
          this%sol%temp(2,1) = czero
          this%sol%temp(3,1) = czero
    
        do ir = 2, this%nd
          is = 3*(ir-1)+1
          
          this%sol%temp(is  ,1) = Temp(ir) / this%dt + this%ntemp(1,ir) + this%htide_fn(ir,1)
          this%sol%temp(is+1,1) = czero
          this%sol%temp(is+2,1) = czero
        end do
    
        ir = this%nd+1
          this%sol%temp(3*this%nd+1,1) = czero
    
        call this%mat%temp(0)%luSolve_sub(this%sol%temp(:,1))
    
        if ( maxval(abs(this%sol%temp_i_fn(0,0) - Temp1)/abs(Temp1)) < 1e-5) exit       
      end do
    deallocate( Temp, Temp1 )
    
    qConv = real(-this%sol%flux_fn(1,0,0,1), kind=dbl) / sqrt(4*pi)
    this%sol%v_dn = this%vr_jm_fn(1) + this%Raf * qConv * flux_bnd(1:this%jms) ; this%sol%v_dn(1) = czero
    
    !$omp parallel do
    do ij = 1, this%jmax
      call this%mat%temp(ij)%fill_sub( this%matica_temp_fn(j_in=ij, a_in=1._dbl), this%matica_temp_fn(j_in=ij, a_in=0._dbl) )
    end do
    !$omp end parallel do
  
    !$omp parallel do private (is,ir)
    do ijm = 2, this%jms
      ir = 1
        this%sol%temp(1,ijm) = -( this%sol%u_dn(ijm) + this%sol%v_dn(ijm) * this%dt )
        this%sol%temp(2,ijm) = czero
        this%sol%temp(3,ijm) = czero
      
      do ir = 2, this%nd
        is = 3*(ir-1)+1
        
        this%sol%temp(is  ,ijm) = this%sol%temp(is,ijm) / this%dt + this%ntemp(ijm,ir) + this%htide_fn(ir,ijm)
        this%sol%temp(is+1,ijm) = czero
        this%sol%temp(is+2,ijm) = czero
      end do
  
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = -(this%sol%u_up(ijm) + this%sol%v_up(ijm) * this%dt)
      
      call this%mat%temp( this%j_indx(ijm) )%luSolve_sub(this%sol%temp(:,ijm))
    end do
    !$omp end parallel do
    
    !$omp parallel do
    do ij = 1, this%jmax
      call this%mat%mech(ij)%fill_sub( this%matica_mech_fn(j_in=ij, a_in=1._dbl), this%matica_mech_fn(j_in=ij, a_in=0._dbl) )
    end do
    !$omp end parallel do

    this%sol%v_dn = - this%Raf * ( this%qr_jm_fn(1) - qConv * flux_bnd(1:this%jms) ) ; this%sol%v_dn(1) = czero

    !$omp parallel do private (ir,is,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)
    
      ir = 1
        this%sol%mech(1,ijm) = -( this%sol%u_dn(ijm) + this%sol%v_dn(ijm) * this%dt - this%Vdelta_fn(1,ijm) )
        this%sol%mech(2,ijm) = czero
        this%sol%mech(3,ijm) = czero
        this%sol%mech(4,ijm) = czero
        this%sol%mech(5,ijm) = czero
        this%sol%mech(6,ijm) = czero
  
      do ir = 2, this%nd
        is = 6*(ir-1)+1
        
        this%sol%mech(is  ,ijm) = -sqrt((ij  )/(2*ij+1._dbl)) * this%buoy_rr_fn(ir,ijm)
        this%sol%mech(is+1,ijm) = +sqrt((ij+1)/(2*ij+1._dbl)) * this%buoy_rr_fn(ir,ijm)
        this%sol%mech(is+2,ijm) = czero
        this%sol%mech(is+3,ijm) = czero
        this%sol%mech(is+4,ijm) = czero
        this%sol%mech(is+5,ijm) = czero
      end do
  
      ir = this%nd+1
        this%sol%mech(6*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+2,ijm) = -( this%sol%u_up(ijm) - this%Vdelta_fn(this%nd,ijm) )
      
      call this%mat%mech(ij)%luSolve_sub(this%sol%mech(:,ijm))
    end do
    !$omp end parallel do
    
    this%sol%v_dn = this%sol%v_dn + this%vr_jm_fn(1)      ; this%sol%v_dn(1) = czero
    this%sol%v_up =                 this%vr_jm_fn(this%nd); this%sol%v_up(1) = czero
    
    this%sol%u_dn = this%sol%u_dn + this%sol%v_dn * this%dt
    this%sol%u_up = this%sol%u_up + this%sol%v_up * this%dt
    
    do ijm = 2, this%jms
      this%sol%t_dn(ijm) = this%sol%u_dn(ijm) - this%Vdelta_fn(1      ,ijm)
      this%sol%t_up(ijm) = this%sol%u_up(ijm) - this%Vdelta_fn(this%nd,ijm)
    end do
    
  end subroutine EE_iceCrust_sub

    complex(kind=dbl) function Vdelta_iceCrust_fn(this, ir, ijm)
      class(T_iceCrust),  intent(in) :: this
      integer,            intent(in) :: ir, ijm
      integer                        :: k, j, m
      real(kind=dbl)                 :: ri
      complex(kind=dbl), allocatable :: field(:)
      
      j  = this%j_indx(ijm)
      m  = ijm - (j*(j+1)/2+1)
      ri = this%rad_grid%r(ir)
      
      allocate( field(this%nd+1) ); field = czero
      
      do k = 1, this%nd+1
        field(k) = -this%rhoI * this%alphaU * this%alpha_fn(k) * (this%Td-this%Tu) * this%sol%temp_fn(k,j,m)
      end do
      
      if (ir == 1) then
        field = field * ( this%rad_grid%r(1) / this%rad_grid%rr(:) ) ** (j-1)
      else if (ir == this%nd) then
        field = field * ( this%rad_grid%rr(:) / this%rad_grid%r(this%nd) ) ** (j+2)
      end if
      
      Vdelta_iceCrust_fn = ( this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI           , this%sol%u_up(ijm) ) + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW-this%rhoI , this%sol%u_dn(ijm) ) + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW, this%sol%u_I2(ijm) ) + &
                           & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC-this%rhoI2, this%sol%u_C(ijm)  ) + &
                           & this%gravity%V_rho_fn( j, m, ri, field, this%rad_grid)                                   + &
                           & this%gravity%V_rt_fn(  j, m, ri ) ) / this%gravity%g_fn( ri )
      
      deallocate( field )
      
    end function Vdelta_iceCrust_fn
    
    subroutine set_layers_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
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
                            & (this%Td-this%Tu) * this%sol%temp_fn(i,j,m) * (this%rI2/this%rad_grid%rr(i))**(j-1)
          end do
          
          rhs1 = -( grv%V_bnd_fn(j,m,this%rI2,ru,this%rhoI          ,this%sol%u_up(jm_int)+this%sol%v_up(jm_int)*this%dt ) + &
                  & grv%V_bnd_fn(j,m,this%rI2,rd,this%rhoW-this%rhoI,this%sol%u_dn(jm_int)+this%sol%v_dn(jm_int)*this%dt ) + &
                  & grv%V_rho_fn(j,m,this%rI2, field, this%rad_grid) + &
                  & grv%V_rt_fn( j,m,this%rI2) )
          
          do i = 1, this%nd+1
            field(i) = -this%rhoI * this%alphaU * this%alpha_fn(i) * &
                              & (this%Td-this%Tu) * this%sol%temp_fn(i,j,m) * (this%rC/this%rad_grid%rr(i))**(j-1)
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
      
    end subroutine set_layers_iceCrust_sub
    
    pure function htide_iceCrust_fn(this, ir, ijm) result(HI)
      class(T_iceCrust), intent(in) :: this
      integer,           intent(in) :: ir, ijm
      complex(kind=dbl)             :: HI
      
      if ( ijm <= jms4 ) then
        HI = this%rad_grid%cc(ir,-1) * this%htide(ir-1,ijm) + this%rad_grid%cc(ir,+1) * this%htide(ir,ijm)
      else
        HI = czero
      end if
      
      HI = this%Ds/this%Ra * HI / this%cp_fn(ir)
      
    end function htide_iceCrust_fn
    
  subroutine vypis_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
      
    call this%vypis_sub(8, 'data/data_ice_topo' , 'topo' )
    call this%vypis_sub(8, 'data/data_ice_shape', 'shape')
    
    this%poc = this%poc + 1
    
  end subroutine vypis_iceCrust_sub
  
  subroutine deallocate_iceCrust_sub(this)
    class(T_iceCrust), intent(inout) :: this
    
    call this%lat_grid%deallocate_fftw_vcvv_sub()
    call this%deallocate_ice_sub()
    
  end subroutine deallocate_iceCrust_sub
  
end module IceCrustMod