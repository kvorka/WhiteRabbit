submodule(OceanMod) OceanEquations
  implicit none

  contains

  subroutine init_eq_temp_sub(this, rhs, nl)
    class(T_ocean), intent(inout) :: this
    logical,        intent(in)    :: rhs, nl
    integer                       :: j
    
    call this%sol%init_stemp_sub()
    call this%mat%init_mtemp_sub()

    allocate( this%flux_up(this%jms) ) ; this%flux_up = czero

    if (rhs) then
      allocate( this%rtemp(this%nd+1,this%jms) ) ; this%rtemp = czero

      this%rtemp(1,1) = cmplx(sqrt(4*pi), 0._dbl, kind=dbl)
    end if

    if (nl) then
      allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp = czero
    end if

    do j=0, this%jmax
      call this%mat%temp(j)%fill_sub( matica_temp_fn(this,j,+0.6_dbl), matica_temp_fn(this,j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_temp_sub
  
  subroutine solve_temp_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: ij, ir, ir1, ijm
    
    !$omp parallel do private (ir,ir1,ij)
    do ijm = 1, this%jms
      ij = this%j_indx(ijm)
      
      do concurrent ( ir=2:this%nd )
        this%rtemp(ir,ijm) = this%rtemp(ir,ijm) + this%mat%temp(ij)%multipl_fn(3*(ir-1)+1, this%sol%temp(:,ijm))
      end do
      
      do concurrent ( ir=1:this%nd )
        ir1 = 3*(ir-1)+1
        
        this%sol%temp(ir1  ,ijm) = this%rtemp(ir,ijm)
        this%sol%temp(ir1+1,ijm) = czero
        this%sol%temp(ir1+2,ijm) = czero
      end do
        
      ir = this%nd+1
        this%sol%temp(3*this%nd+1,ijm) = this%rtemp(ir,ijm)
        
      call this%mat%temp(ij)%luSolve_sub( this%sol%temp(:,ijm) )
    end do
    !$omp end parallel do
    
  end subroutine solve_temp_sub
  
  subroutine init_eq_torr_sub(this, rhs, nl)
    class(T_ocean), intent(inout) :: this
    logical,        intent(in)    :: rhs, nl
    integer                       :: j
    
    call this%sol%init_storr_sub()
    call this%mat%init_mtorr_sub()

    if (rhs) then
      allocate( this%rtorr(this%nd+1,this%jms) ) ; this%rtorr = czero
    end if

    if (nl) then
      allocate( this%ntorr(this%jms,2:this%nd) ) ; this%ntorr = czero
    end if

    do j=1, this%jmax
      call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+0.6_dbl), matica_torr_fn(this,j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_torr_sub

  subroutine solve_torr_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: ij, ijm, ir, ir1

    !$omp parallel do private (ir,ir1,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      do concurrent ( ir=2:this%nd )
        this%rtorr(ir,ijm) = this%rtorr(ir,ijm) + this%mat%torr(ij)%multipl_fn(3*(ir-1)+1,this%sol%torr(:,ijm))
      end do

      do concurrent ( ir=1:this%nd )
        ir1 = 3*(ir-1) + 1
        
        this%sol%torr(ir1  ,ijm) = this%rtorr(ir,ijm)
        this%sol%torr(ir1+1,ijm) = czero
        this%sol%torr(ir1+2,ijm) = czero
      end do
  
      ir = this%nd+1
        this%sol%torr(3*this%nd+1,ijm) = this%rtorr(ir,ijm)
        
      call this%mat%torr(ij)%luSolve_sub( this%sol%torr(:,ijm) )
    end do
    !$omp end parallel do

  end subroutine solve_torr_sub

  subroutine init_eq_mech_sub(this, rhs, nl)
    class(T_ocean), intent(inout) :: this
    logical,        intent(in)    :: rhs, nl
    integer                       :: j
    
    call this%sol%init_smech_sub()
    call this%mat%init_mmech_sub()

    if (rhs) then
      allocate( this%rsph1(this%nd+1,this%jms) ) ; this%rsph1 = czero
      allocate( this%rsph2(this%nd+1,this%jms) ) ; this%rsph2 = czero
    end if

    if (nl) then
      allocate( this%nsph1(this%jms,2:this%nd) ) ; this%nsph1 = czero
      allocate( this%nsph2(this%jms,2:this%nd) ) ; this%nsph2 = czero
    end if

    do j=1, this%jmax
      call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+0.6_dbl), matica_mech_fn(this,j,-0.4_dbl) )
    end do

  end subroutine init_eq_mech_sub

  subroutine solve_mech_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: ij, ijm, ir, ir1

    !$omp parallel do private (ir,ir1,ij)
    do ijm = 2, this%jms
      ij = this%j_indx(ijm)

      do concurrent ( ir=2:this%nd )
        this%rsph1(ir,ijm) = this%rsph1(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+1,this%sol%mech(:,ijm))
        this%rsph2(ir,ijm) = this%rsph2(ir,ijm) + this%mat%mech(ij)%multipl_fn(6*(ir-1)+2,this%sol%mech(:,ijm))
      end do
      
      do concurrent ( ir=1:this%nd )
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

  end subroutine solve_mech_sub

  subroutine init_bnd_deformation_sub(this)
    class(T_ocean), intent(inout) :: this
    
    call this%sol%init_layer_u_sub()
    
  end subroutine init_bnd_deformation_sub

end submodule OceanEquations