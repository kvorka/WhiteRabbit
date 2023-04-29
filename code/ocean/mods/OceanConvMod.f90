module OceanConvMod
  use OceanMod
  implicit none
  
  type, extends(T_ocean), public :: T_oceanConv
    contains
    
    procedure, public, pass :: init_sub        => init_oceanConv_sub
    procedure, public, pass :: time_scheme_sub => time_scheme_oceanConv_sub
  end type T_oceanConv
  
  contains
  
  subroutine init_oceanConv_sub(this)
    class(T_oceanConv), intent(inout) :: this
    
    call this%init_ocean_sub() ; call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()
    
    call this%init_eq_temp_sub(rhs=.true., nl=.true.)
    call this%init_eq_torr_sub(rhs=.true., nl=.true.)
    call this%init_eq_mech_sub(rhs=.true., nl=.true.)
    
    call this%init_state_sub()
    
  end subroutine init_oceanConv_sub
  
  subroutine time_scheme_oceanConv_sub(this, cf)
    class(T_oceanConv), intent(inout) :: this
    real(kind=dbl),     intent(in)    :: cf
    integer                           :: ir, ir1, ij, ijm

    ij = 0
      do ir = 2, this%nd
        this%rtemp(ir,1) = (1-cf) * this%ntemp(1,ir) + this%mat%temp(0)%multipl_fn(3*(ir-1)+1,this%sol%temp(:,1))
      end do
    
    !$omp parallel do private (ir1,ij) collapse (2)
    do ijm = 2, this%jms
      do ir = 2, this%nd
        ij = this%j_indx(ijm)
        
        ir1 = 3*(ir-1)+1
          this%rtemp(ir,ijm) = (1-cf) * this%ntemp(ijm,ir) + this%mat%temp(ij)%multipl_fn(ir1  ,this%sol%temp(:,ijm))
          this%rtorr(ir,ijm) = (1-cf) * this%ntorr(ijm,ir) + this%mat%torr(ij)%multipl_fn(ir1  ,this%sol%torr(:,ijm))
        
        ir1 = 6*(ir-1)+1
          this%rsph1(ir,ijm) = (1-cf) * this%nsph1(ijm,ir) + this%mat%mech(ij)%multipl_fn(ir1  ,this%sol%mech(:,ijm))
          this%rsph2(ir,ijm) = (1-cf) * this%nsph2(ijm,ir) + this%mat%mech(ij)%multipl_fn(ir1+1,this%sol%mech(:,ijm))
      end do
    end do
    !$omp end parallel do

    !$omp parallel do
    do ir = 2, this%nd
      call fullnl_sub(this, ir, this%ntemp(:,ir), this%nsph1(:,ir), this%ntorr(:,ir), this%nsph2(:,ir))
    end do
    !$omp end parallel do

    !$omp parallel do collapse (2)
    do ijm = 2, this%jms
      do ir = 2, this%nd
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + cf * this%ntemp(ijm,ir)
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + cf * this%ntorr(ijm,ir)
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + cf * this%nsph1(ijm,ir)
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + cf * this%nsph2(ijm,ir)
      end do
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
      
      call this%mat%temp(0)%luSolve_sub( this%sol%temp(:,1) )
      
    !$omp parallel do private (ir,ir1,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)
      
      ir = 1
        this%sol%temp(1,ijm) = czero
        this%sol%temp(2,ijm) = czero
        this%sol%temp(3,ijm) = czero
      
      do concurrent ( ir=2:this%nd )
        ir1 = 3 * (ir-1) + 1
        
        this%sol%temp(ir1  ,ijm) = this%rtemp(ir,ijm)
        this%sol%temp(ir1+1,ijm) = czero
        this%sol%temp(ir1+2,ijm) = czero
      end do
      
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = czero
      
      call this%mat%temp(ij)%luSolve_sub( this%sol%temp(:,ijm) )

      ir = 1
        this%sol%torr(1,ijm) = czero
        this%sol%torr(2,ijm) = czero
        this%sol%torr(3,ijm) = czero
      
      do concurrent ( ir=2:this%nd )
        ir1 = 3 * (ir-1) + 1

        this%sol%torr(ir1  ,ijm) = this%rtorr(ir,ijm)
        this%sol%torr(ir1+1,ijm) = czero
        this%sol%torr(ir1+2,ijm) = czero
      end do

      ir = this%nd+1
        this%sol%torr(3*this%nd+1,ijm) = czero
      
      call this%mat%torr(ij)%luSolve_sub( this%sol%torr(:,ijm) )

      ir = 1
        this%sol%mech(1,ijm) = czero
        this%sol%mech(2,ijm) = czero
        this%sol%mech(3,ijm) = czero
        this%sol%mech(4,ijm) = czero
        this%sol%mech(5,ijm) = czero
        this%sol%mech(6,ijm) = czero
        
      do concurrent ( ir=2:this%nd )
        ir1 = 6*(ir-1)+1

        this%sol%mech(ir1  ,ijm) = this%rsph1(ir,ijm)
        this%sol%mech(ir1+1,ijm) = this%rsph2(ir,ijm)
        this%sol%mech(ir1+2,ijm) = czero
        this%sol%mech(ir1+3,ijm) = czero
        this%sol%mech(ir1+4,ijm) = czero
        this%sol%mech(ir1+5,ijm) = czero
      end do
      
      ir = this%nd+1
        this%sol%mech(6*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+2,ijm) = czero
      
      call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
    end do
    !$omp end parallel do
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()
    
  end subroutine time_scheme_oceanConv_sub
  
end module OceanConvMod