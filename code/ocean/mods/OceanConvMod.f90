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
    
    call this%init_ocean_sub()
      call this%lat_grid%init_vcsv_vcvv_vcvgv_sub()
    
    call this%init_eq_temp_sub()
      allocate( this%ntemp(2:this%nd,this%jms) ) ; this%ntemp   = czero
      allocate( this%rtemp(2:this%nd,this%jms) ) ; this%rtemp   = czero
      allocate( this%flux_up(this%jms)         ) ; this%flux_up = czero

    call this%init_eq_torr_sub()
    call this%init_eq_mech_sub()
      allocate( this%nmech(2:this%nd,this%jmv) ) ; this%nmech   = czero
      allocate( this%rmech(2:this%nd,this%jmv) ) ; this%rmech   = czero
    
    call this%init_state_sub()
    
  end subroutine init_oceanConv_sub
  
  subroutine time_scheme_oceanConv_sub(this, cf)
    class(T_oceanConv), intent(inout) :: this
    real(kind=dbl),     intent(in)    :: cf
    integer                           :: ir, ir1, ir2, ij, ijm, ijml

    ij = 0
      do ir = 2, this%nd
        this%rtemp(ir,1) = (1-cf) * this%ntemp(ir,1) + this%mat%temp(0)%multipl_fn(3*(ir-1)+1,this%sol%temp(:,1))
      end do
    
    !$omp parallel do private (ir1,ir2,ij,ijml) collapse(2)
    do ijm = 2, this%jms
      do ir = 2, this%nd
        ij   = this%j_indx(ijm)
        ijml = 3*(ijm-1)
        
        ir1 = 3*(ir-1)+1
        ir2 = 6*(ir-1)+1

        this%rtemp(ir,ijm   ) = (1-cf) * this%ntemp(ir,ijm   ) + this%mat%temp(ij)%multipl_fn(ir1  ,this%sol%temp(:,ijm))
        this%rmech(ir,ijml-1) = (1-cf) * this%nmech(ir,ijml-1) + this%mat%mech(ij)%multipl_fn(ir2  ,this%sol%mech(:,ijm))
        this%rmech(ir,ijml  ) = (1-cf) * this%nmech(ir,ijml  ) + this%mat%torr(ij)%multipl_fn(ir1  ,this%sol%torr(:,ijm))
        this%rmech(ir,ijml+1) = (1-cf) * this%nmech(ir,ijml+1) + this%mat%mech(ij)%multipl_fn(ir2+1,this%sol%mech(:,ijm))
      end do
    end do
    !$omp end parallel do

    !$omp parallel do
    do ir = 2, this%nd
      call fullnl_sub(this, ir, this%ntemp(ir,:), this%nmech(ir,:))
    end do
    !$omp end parallel do
    
    ij = 0
      ir = 1
        this%sol%temp(1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
        this%sol%temp(2,1) = czero
        this%sol%temp(3,1) = czero
      
      do ir = 2, this%nd
        ir1 = 3*(ir-1)+1
          
        this%sol%temp(ir1  ,1) = this%rtemp(ir,1) + cf * this%ntemp(ir,1)
        this%sol%temp(ir1+1,1) = czero
        this%sol%temp(ir1+2,1) = czero
      end do
      
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,1) = czero
      
      call this%mat%temp(0)%luSolve_sub( this%sol%temp(:,1) )
      
    !$omp parallel do private (ir,ir1,ir2,ij,ijml)
    do ijm = 2, this%jms
      ij   = this%j_indx(ijm)
      ijml = 3*(ijm-1)
      
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

        this%sol%temp(ir1  ,ijm) = this%rtemp(ir,ijm   ) + cf * this%ntemp(ir,ijm   )
        this%sol%torr(ir1  ,ijm) = this%rmech(ir,ijml  ) + cf * this%nmech(ir,ijml  )
        this%sol%mech(ir2  ,ijm) = this%rmech(ir,ijml-1) + cf * this%nmech(ir,ijml-1)
        this%sol%mech(ir2+1,ijm) = this%rmech(ir,ijml+1) + cf * this%nmech(ir,ijml+1)
        
        this%sol%temp(ir1+1:ir1+2,ijm) = czero
        this%sol%torr(ir1+1:ir1+2,ijm) = czero
        this%sol%mech(ir2+2:ir2+5,ijm) = czero
      end do
      
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = czero
        this%sol%torr(3*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+1,ijm) = czero
        this%sol%mech(6*this%nd+2,ijm) = czero
        
      call this%mat%temp(ij)%luSolve_sub( this%sol%temp(:,ijm) )
      call this%mat%torr(ij)%luSolve_sub( this%sol%torr(:,ijm) )
      call this%mat%mech(ij)%luSolve_sub( this%sol%mech(:,ijm) )
    end do
    !$omp end parallel do
    
    if (this%mechanic_bnd == 'frees') call this%global_rotation_sub()
    
  end subroutine time_scheme_oceanConv_sub
  
end module OceanConvMod