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
      allocate( this%rtemp(2:this%nd,this%jms) ) ; this%rtemp = czero
    end if

    if (nl) then
      allocate( this%ntemp(this%jms,2:this%nd) ) ; this%ntemp = czero
    end if

    do j=0, this%jmax
      call this%mat%temp(j)%fill_sub( matica_temp_fn(this,j,+0.6_dbl), matica_temp_fn(this,j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_temp_sub

  subroutine init_eq_torr_sub(this, rhs, nl)
    class(T_ocean), intent(inout) :: this
    logical,        intent(in)    :: rhs, nl
    integer                       :: j
    
    call this%sol%init_storr_sub()
    call this%mat%init_mtorr_sub()

    if (rhs) then
      allocate( this%rtorr(2:this%nd,this%jms) ) ; this%rtorr = czero
    end if

    if (nl) then
      allocate( this%ntorr(this%jms,2:this%nd) ) ; this%ntorr = czero
    end if

    do j=1, this%jmax
      call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+0.6_dbl), matica_torr_fn(this,j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_torr_sub

  subroutine init_eq_mech_sub(this, rhs, nl)
    class(T_ocean), intent(inout) :: this
    logical,        intent(in)    :: rhs, nl
    integer                       :: j
    
    call this%sol%init_smech_sub()
    call this%mat%init_mmech_sub()

    if (rhs) then
      allocate( this%rsph1(2:this%nd,this%jms) ) ; this%rsph1 = czero
      allocate( this%rsph2(2:this%nd,this%jms) ) ; this%rsph2 = czero
    end if

    if (nl) then
      allocate( this%nsph1(this%jms,2:this%nd) ) ; this%nsph1 = czero
      allocate( this%nsph2(this%jms,2:this%nd) ) ; this%nsph2 = czero
    end if

    do j=1, this%jmax
      call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+0.6_dbl), matica_mech_fn(this,j,-0.4_dbl) )
    end do

  end subroutine init_eq_mech_sub

  subroutine init_bnd_deformation_sub(this)
    class(T_ocean), intent(inout) :: this
    
    call this%sol%init_layer_u_sub()
    
  end subroutine init_bnd_deformation_sub

end submodule OceanEquations