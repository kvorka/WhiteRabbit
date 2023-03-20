submodule(OceanMod) OceanEquations
  implicit none

  contains

  subroutine init_eq_temp_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: j
    
    call this%sol%init_stemp_sub()
    call this%mat%init_mtemp_sub()

    do j=0, this%jmax
      call this%mat%temp(j)%fill_sub( matica_temp_fn(this,j,+0.6_dbl), matica_temp_fn(this,j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_temp_sub

  subroutine init_eq_torr_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: j
    
    call this%sol%init_storr_sub()
    call this%mat%init_mtorr_sub()

    do j=1, this%jmax
      call this%mat%torr(j)%fill_sub( matica_torr_fn(this,j,+0.6_dbl), matica_torr_fn(this,j,-0.4_dbl) )
    end do
    
  end subroutine init_eq_torr_sub

  subroutine init_eq_mech_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: j
    
    call this%sol%init_smech_sub()
    call this%mat%init_mmech_sub()

    do j=1, this%jmax
      call this%mat%mech(j)%fill_sub( matica_mech_fn(this,j,+0.6_dbl), matica_mech_fn(this,j,-0.4_dbl) )
    end do

  end subroutine init_eq_mech_sub

  subroutine init_bnd_deformation_sub(this)
    class(T_ocean), intent(inout) :: this
    
    call this%sol%init_layer_u_sub()
    
  end subroutine init_bnd_deformation_sub

end submodule OceanEquations