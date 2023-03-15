module IceCrustMod
  use Math
  use IceMod
  use IceTidesMod
  use NonLinearTerms
  use MatrixDefinitions
  implicit none
  
  type, extends(T_ice), public :: T_iceCrust
    complex(kind=dbl), allocatable :: ntemp(:,:)
    
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
    
    call tides%init_sub()
    call tides%deallocate_sub()

    allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp = cmplx(0._dbl, 0._dbl, kind=dbl)

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

    call find_hydrostatic_iceCrust_sub(this) ; call vypis_iceCrust_sub(this)

  end subroutine init_iceCrust_sub

    subroutine find_hydrostatic_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
      complex(kind=dbl), allocatable   :: flux_bnd(:)
      
      allocate( flux_bnd(this%jms) ); flux_bnd = cmplx(0._dbl, 0._dbl, kind=dbl)
      
      do
        call EE_iceCrust_sub(this, flux_bnd)
        if ( maxval(abs(this%sol%v_up * this%dt / this%sol%u_up)) < 1e-3 ) exit
      end do
      
      deallocate( flux_bnd )
      
      write(*,*) 'Icy crust quasi hydrostatic'
      
    end subroutine find_hydrostatic_iceCrust_sub
  
  subroutine iter_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust),               intent(inout) :: this
    complex(kind=dbl), dimension(:), intent(in)    :: flux_bnd
    integer                                        :: n
    
    do n = 1, this%n_iter
      call EE_iceCrust_sub(this, flux_bnd)
    end do

    call vypis_iceCrust_sub(this)
    
  end subroutine iter_iceCrust_sub
  
  subroutine EE_iceCrust_sub(this, flux_bnd)
    class(T_iceCrust), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: flux_bnd(:)
    integer                          :: i, j, m, jm_int, jm1
    real(kind=dbl)                   :: qConv, Ra_g_alpha
    complex(kind=dbl), allocatable   :: Temp(:)
    
    this%t = this%t + this%dt

    !$omp parallel do private (j, jm1)
    do i = 2, this%nd
      j = 0 ; this%ntemp(1,i) = this%mat%temp(0)%multipl_fn(3*(i-1)+1,this%sol%temp(:,1))
      
      do j = 1, this%jmax
        jm1 = j*(j+1)/2+1 ; this%ntemp(jm1:jm1+j,i) = this%mat%temp(j)%multipl2_fn(3*(i-1)+1,this%sol%temp(:,jm1:jm1+j))
      end do

      this%ntemp(:,i) = this%ntemp(:,i) - vgradT_fn(this,i)
    end do
    !$omp end parallel do
      
    allocate( Temp(this%nd+1) )
      do
        Temp = this%sol%temp_i_fn(0,0)
        call this%mat%temp(0)%fill_sub( matica_temp_fn(this, j_in=0, a_in=1._dbl), matica_temp_fn(this, j_in=0, a_in=0._dbl) )
        
        i = 1
          this%sol%temp(3*(i-1)+1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
          this%sol%temp(3*(i-1)+2,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%temp(3*(i-1)+3,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
      
        do i = 2, this%nd
          this%sol%temp(3*(i-1)+1,1) = this%ntemp(1,i) + this%Ds/this%Ra * this%htide_fn(i,1) / this%cp_fn(i)
          this%sol%temp(3*(i-1)+2,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
          this%sol%temp(3*(i-1)+3,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
        end do
      
        i = this%nd+1
          this%sol%temp(3*(i-1)+1,1) = cmplx(0._dbl, 0._dbl, kind=dbl)
      
        call this%mat%temp(0)%luSolve_sub(this%sol%temp(:,1))
        if ( maxval(abs(this%sol%temp_i_fn(0,0) - Temp)/abs(Temp)) < 1e-5 ) exit       
      end do
    deallocate( Temp )
    
    qConv = real(-this%sol%flux_fn(1,0,0,1), kind=dbl) / sqrt(4*pi)

    !$omp parallel do
    do j = 1, this%jmax
      call this%mat%temp(j)%fill_sub( matica_temp_fn(this, j_in=j, a_in=1._dbl), matica_temp_fn(this, j_in=j, a_in=0._dbl) )
      call this%mat%mech(j)%fill_sub( matica_mech_fn(this, j_in=j, a_in=1._dbl), matica_mech_fn(this, j_in=j, a_in=0._dbl) )
    end do
    !$omp end parallel do

    this%sol%v_dn = this%vr_jm_fn(1) + this%Raf * qConv * flux_bnd(1:this%jms)
    
    !$omp parallel do private (i)
    do jm_int = 2, this%jms
      i = 1
        this%sol%temp( 3*(i-1)+1           , jm_int ) = -( this%sol%u_dn(jm_int) + this%sol%v_dn(jm_int) * this%dt )
        this%sol%temp( 3*(i-1)+2:3*(i-1)+3 , jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
    
      do i = 2, this%nd
        this%sol%temp( 3*(i-1)+1          , jm_int ) = this%ntemp(jm_int,i) + this%Ds/this%Ra*this%htide_fn(i,jm_int)/this%cp_fn(i)
        this%sol%temp( 3*(i-1)+2:3*(i-1)+3, jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
      end do
    
      i = this%nd+1
        this%sol%temp( 3*(i-1)+1, jm_int ) = -(this%sol%u_up(jm_int) + this%sol%v_up(jm_int) * this%dt)
        
      call this%mat%temp(this%j_indx(jm_int))%luSolve_sub(this%sol%temp(:,jm_int))
    end do
    !$omp end parallel do

    this%sol%v_dn = - this%Raf * ( this%qr_jm_fn(1) - qConv * flux_bnd(1:this%jms) )

    !$omp parallel do private (i,j,m,Ra_g_alpha)
    do jm_int = 2, this%jms
      j  = this%j_indx(jm_int)
      m  = jm_int - (j*(j+1)/2+1)
    
      i = 1
        this%sol%mech( 6*(i-1)+1 ,jm_int ) = -( this%sol%u_dn(jm_int) + this%sol%v_dn(jm_int)*this%dt - this%Vdelta_fn(j,m,1) )
        this%sol%mech( 6*(i-1)+2 :   &
                     & 6*(i-1)+6 ,jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
  
      do i = 2, this%nd
        Ra_g_alpha = this%Ra * this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i))
      
        this%sol%mech( 6*(i-1)+1 ,jm_int ) = -Ra_g_alpha * sqrt((j  )/(2*j+1._dbl)) * this%sol%temp_fn(i,j,m)
        this%sol%mech( 6*(i-1)+2 ,jm_int ) = +Ra_g_alpha * sqrt((j+1)/(2*j+1._dbl)) * this%sol%temp_fn(i,j,m)
        this%sol%mech( 6*(i-1)+3 :   &
                     & 6*(i-1)+6 ,jm_int ) = cmplx(0._dbl, 0._dbl, kind=dbl)
      end do
  
      i = this%nd+1
        this%sol%mech(6*(i-1)+1,jm_int) = cmplx(0._dbl, 0._dbl, kind=dbl)
        this%sol%mech(6*(i-1)+2,jm_int) = -( this%sol%u_up(jm_int) - this%Vdelta_fn(j,m,this%nd) )
      
      call this%mat%mech(j)%luSolve_sub(this%sol%mech(:,jm_int))
    end do
    !$omp end parallel do
    
    this%sol%v_dn = this%sol%v_dn + this%vr_jm_fn(1);       this%sol%v_dn(1) = cmplx(0._dbl, 0._dbl, kind=dbl)
    this%sol%v_up =                 this%vr_jm_fn(this%nd); this%sol%v_up(1) = cmplx(0._dbl, 0._dbl, kind=dbl)
    
      this%sol%u_dn = this%sol%u_dn + this%sol%v_dn * this%dt
      this%sol%u_up = this%sol%u_up + this%sol%v_up * this%dt
    
    do j = 1, this%jmax
      do m = 0, j
        jm_int = j*(j+1)/2+m+1
        
        this%sol%t_dn(jm_int) = this%sol%u_dn(jm_int) - this%Vdelta_fn(j,m,      1)
        this%sol%t_up(jm_int) = this%sol%u_up(jm_int) - this%Vdelta_fn(j,m,this%nd)
      end do
    end do

  end subroutine EE_iceCrust_sub
  
    complex(kind=dbl) function Vdelta_iceCrust_fn(this, j, m, i)
      class(T_iceCrust),  intent(in) :: this
      integer,            intent(in) :: j, m, i
      integer                        :: k
      complex(kind=dbl), allocatable :: field(:)
      
      allocate( field(this%nd+1) ); field = cmplx(0._dbl, 0._dbl, kind=dbl)

        do k = 1, this%nd+1
          field(k) = -this%rhoI * this%alphaU * this%alpha_fn(k) * (this%Td-this%Tu) * this%sol%temp_fn(k,j,m)
        end do
        
        if (i == 1) then
          field = field * ( this%rd / this%rad_grid%rr(:) ) ** (j-1)
        else if (i == this%nd) then
          field = field * ( this%rad_grid%rr(:) / this%ru ) ** (j+2)
        end if

        associate( ri => this%rad_grid%r(i) )
          Vdelta_iceCrust_fn = this%gravity%V_bnd_fn( j, m, ri, this%ru , this%rhoI           , this%sol%u_up(jm(j,m)) ) + &
                             & this%gravity%V_bnd_fn( j, m, ri, this%rd , this%rhoW-this%rhoI , this%sol%u_dn(jm(j,m)) ) + &
                             & this%gravity%V_bnd_fn( j, m, ri, this%rI2, this%rhoI2-this%rhoW, this%sol%u_I2(jm(j,m)) ) + &
                             & this%gravity%V_bnd_fn( j, m, ri, this%rC , this%rhoC-this%rhoI2, this%sol%u_C(jm(j,m))  ) + &
                             & this%gravity%V_rho_fn( j, m, ri, field, this%rad_grid)                                    + &
                             & this%gravity%V_rt_fn(  j, m, ri )
        end associate

      deallocate( field )

      Vdelta_iceCrust_fn = Vdelta_iceCrust_fn / this%gravity%g_fn(this%rad_grid%r(i))
      
    end function Vdelta_iceCrust_fn
    
    subroutine set_layers_iceCrust_sub(this)
      class(T_iceCrust), intent(inout) :: this
      integer                          :: i, j, m, jm_int
      real(kind=dbl)                   :: a11, a12, a21, a22, det
      complex(kind=dbl)                :: rhs1, rhs2
      complex(kind=dbl), allocatable   :: field(:)
      
      rhs1 = cmplx(0._dbl, 0._dbl, kind=dbl)
      rhs2 = cmplx(0._dbl, 0._dbl, kind=dbl)
      
      associate( rd => this%rad_grid%r(1), ru => this%rad_grid%r(this%nd), grv => this%gravity  )
      allocate( field(this%nd+1) )

      do j = 1, this%jmax
        a11 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g - this%gravity%g_fn(this%rI2)
        a12 = 4 * pi * kappa * this%rI2 * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g * (this%rC/this%rI2)**(j+2)
        a21 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoI2-this%rhoW)  / (2*j+1) / this%g * (this%rC/this%rI2)**(j-1)
        a22 = 4 * pi * kappa * this%rC  * this%D_ud * (this%rhoC -this%rhoI2) / (2*j+1) / this%g - this%gravity%g_fn(this%rC)
        
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
    
    pure function htide_iceCrust_fn(this, i, jm_int) result(HI)
      class(T_iceCrust), intent(in) :: this
      integer,           intent(in) :: i, jm_int
      complex(kind=dbl)             :: HI
      
      HI = cmplx(0._dbl, 0._dbl, kind=dbl)
      if (jm_int <= jms4) HI = this%rad_grid%cc(i,-1) * this%htide(i-1,jm_int) + &
                            &  this%rad_grid%cc(i,+1) * this%htide(i  ,jm_int)
    
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