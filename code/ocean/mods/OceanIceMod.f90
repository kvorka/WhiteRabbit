module OceanIceMod
  use OceanMod
  implicit none

  type, extends(T_ocean), public :: T_oceanice
    real(kind=dbl), private :: ClRoc2

    contains

    procedure, public, pass :: init_sub        => init_oceanice_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanice_sub
  end type T_oceanice

  contains

  subroutine init_oceanice_sub(this)
    class(T_oceanice), intent(inout) :: this
    
    call this%init_ocean_sub() ; call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()

    this%ClRoc2 = Cl_ocean / ( this%Ra * this%Ek**2 / this%Pr )
    
    call this%init_bnd_deformation_sub()
    call this%init_eq_temp_sub(rhs=.true., nl=.true.)
    call this%init_eq_torr_sub(rhs=.true., nl=.true.)
    call this%init_eq_mech_sub(rhs=.true., nl=.true.)
    
    call this%init_state_sub()

  end subroutine init_oceanice_sub

  subroutine time_scheme_oceanice_sub(this, cf)
    class(T_oceanice), intent(inout) :: this
    real(kind=dbl),    intent(in)    :: cf
    integer                          :: ir, ir1, ir2, ij, ijm
    real(kind=dbl)                   :: q

    ij = 0
      do ir = 2, this%nd
        this%rtemp(ir,1) = (1-cf) * this%ntemp(1,ir) + this%mat%temp(0)%multipl_fn(3*(ir-1)+1,this%sol%temp(:,1))
      end do
    
    !$omp parallel do private (ir1,ir2,ij) collapse(2)
    do ijm = 2, this%jms
      do ir = 2, this%nd
        ij = this%j_indx(ijm)
        
        ir1 = 3*(ir-1)+1
        ir2 = 6*(ir-1)+1

        this%rtemp(ir,ijm) = (1-cf) * this%ntemp(ijm,ir) + this%mat%temp(ij)%multipl_fn(ir1  ,this%sol%temp(:,ijm))
        this%rtorr(ir,ijm) = (1-cf) * this%ntorr(ijm,ir) + this%mat%torr(ij)%multipl_fn(ir1  ,this%sol%torr(:,ijm))
        this%rsph1(ir,ijm) = (1-cf) * this%nsph1(ijm,ir) + this%mat%mech(ij)%multipl_fn(ir2  ,this%sol%mech(:,ijm))
        this%rsph2(ir,ijm) = (1-cf) * this%nsph2(ijm,ir) + this%mat%mech(ij)%multipl_fn(ir2+1,this%sol%mech(:,ijm))
      end do
    end do
    !$omp end parallel do

    !$omp parallel do
    do ir = 2, this%nd
      call fullnl_sub(this, ir, this%ntemp(:,ir), this%nsph1(:,ir), this%ntorr(:,ir), this%nsph2(:,ir))
    end do
    !$omp end parallel do
    
    ij = 0
      ir = 1
        this%sol%temp(1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
        this%sol%temp(2,1) = czero
        this%sol%temp(3,1) = czero
      
      do ir = 2, this%nd
        ir1 = 3*(ir-1)+1
          
        this%sol%temp(ir1  ,1) = this%rtemp(ir,1) + cf * this%ntemp(1,ir)
        this%sol%temp(ir1+1,1) = czero
        this%sol%temp(ir1+2,1) = czero
      end do
      
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,1) = czero
      
      call this%mat%temp(0)%luSolve_sub( this%sol%temp(:,1) ) ; q = real(-this%sol%flux_fn(this%nd,0,0,+1), kind=dbl) / sqrt(4*pi)
      
    !$omp parallel do private (ir,ir1,ir2,ij)
    do ijm = 2, this%jms
      ir = 1
        this%sol%temp(1,ijm) = czero
        this%sol%torr(1,ijm) = czero
        this%sol%mech(1,ijm) = czero
        this%sol%mech(2,ijm) = czero
        
        this%sol%temp(2:3,ijm) = czero
        this%sol%torr(2:3,ijm) = czero
        this%sol%mech(3:6,ijm) = czero
        
      do ir = 2, this%nd
        ir1 = 3*(ir-1)+1
        ir2 = 6*(ir-1)+1
        
        this%sol%temp(ir1  ,ijm) = this%rtemp(ir,ijm) + cf * this%ntemp(ijm,ir)
        this%sol%torr(ir1  ,ijm) = this%rtorr(ir,ijm) + cf * this%ntorr(ijm,ir)
        this%sol%mech(ir2  ,ijm) = this%rsph1(ir,ijm) + cf * this%nsph1(ijm,ir)
        this%sol%mech(ir2+1,ijm) = this%rsph2(ir,ijm) + cf * this%nsph2(ijm,ir)
        
        this%sol%temp(ir1+1:ir1+2,ijm) = czero
        this%sol%torr(ir1+1:ir1+2,ijm) = czero
        this%sol%mech(ir2+2:ir2+5,ijm) = czero
      end do
      
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = this%ClRoc2 * this%sol%t_up(ijm) + q * this%sol%u_up(ijm)
        this%sol%torr(3*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+2,ijm) = czero
        
      ij = this%j_indx(ijm)
        call this%mat%temp(ij)%luSolve_sub( this%sol%temp(:,ijm) )
        call this%mat%torr(ij)%luSolve_sub( this%sol%torr(:,ijm) )
        call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
    end do
    !$omp end parallel do
      
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()

  end subroutine time_scheme_oceanice_sub

end module OceanIceMod