submodule(PhysicalObject) Forces
  implicit none
  
  contains
  
  pure function buoy_rr_jml_fn(this, i) result(gdrho)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: gdrho(this%jmv)
    
    gdrho = ersv_fn(this%jmax, this%Ra * this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i)) * this%sol%temp_jm_fn(i))
    
  end function buoy_rr_jml_fn

  pure function coriolis_rr_jml_fn(this, i, v) result(coriolis)
    class(T_physicalObject),     intent(in) :: this
    integer,           optional, intent(in) :: i
    complex(kind=dbl), optional, intent(in) :: v(:)
    complex(kind=dbl)                       :: coriolis(this%jmv)
  
    if ( present(i) ) coriolis = ezvv_fn( this%jmax, this%sol%velocity_jml_fn(i) )  * 2 / this%Ek
    if ( present(v) ) coriolis = ezvv_fn( this%jmax, v )  * 2 / this%Ek

  end function coriolis_rr_jml_fn
  
end submodule Forces