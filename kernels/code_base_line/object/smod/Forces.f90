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
    real(kind=dbl)                          :: fac
    complex(kind=dbl)                       :: coriolis(this%jmv)
  
    if ( present(i) ) coriolis = ezvv_fn( this%jmax, this%sol%velocity_jml_fn(i) )
    if ( present(v) ) coriolis = ezvv_fn( this%jmax, v )

    fac = 2 / this%Ek ; coriolis = coriolis * fac

  end function coriolis_rr_jml_fn

  subroutine coriolis_rr_jml_sub(this, v, coriolis) 
    class(T_physicalObject), intent(in)    :: this
    complex(kind=dbl),       intent(in)    :: v(:)
    complex(kind=dbl),       intent(inout) :: coriolis(:,:)
    real(kind=dbl)                         :: fac
    
    call ezvv_sub(this%jmax, 2/this%Ek, v, coriolis)

  end subroutine coriolis_rr_jml_sub

  subroutine buoy_rr_jml_sub(this, i, force)
    class(T_physicalObject), intent(in)    :: this
    integer,                 intent(in)    :: i
    complex(kind=dbl),       intent(inout) :: force(:,:)
    integer                                :: ijm, j
    real(kind=dbl)                         :: fac, fac1, fac2
    complex(kind=dbl),       allocatable   :: gdrho(:)
    
    allocate( gdrho(this%jms) ) ; call this%sol%temp_jm_sub(i, gdrho)
      
      fac = this%Ra * this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i))
      
      do j = 1, this%jmax
        fac1 = -sqrt( (j  ) / (2*j+1._dbl) ) * fac
        fac2 = +sqrt( (j+1) / (2*j+1._dbl) ) * fac

        do ijm = j*(j+1)/2+1, j*(j+1)/2+j+1
          force(1,ijm) = force(1,ijm) + fac1 * gdrho(ijm)
          force(3,ijm) = force(3,ijm) + fac2 * gdrho(ijm)
        end do
      end do
    
    deallocate( gdrho )
    
  end subroutine buoy_rr_jml_sub
  
end submodule Forces
