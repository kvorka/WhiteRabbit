submodule(PhysicalObject) Variables
  implicit none

  contains

  pure complex(kind=dbl) function htide_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm

    htide_fn = this%rad_grid%cc(ir,-1) * this%htide(ir-1,ijm) + this%rad_grid%cc(ir,+1) * this%htide(ir,ijm)

  end function htide_fn

  pure complex(kind=dbl) function vr_fn(this, i, j, m)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i, j, m
    real(kind=dbl)                      :: jj
    
    jj = real(j, kind=dbl)

    vr_fn = sqrt( (jj  ) / (2*jj+1) ) * this%sol%velocity_fn(i  ,j,m,-1) * this%rad_grid%c(i,-1) - &
          & sqrt( (jj+1) / (2*jj+1) ) * this%sol%velocity_fn(i  ,j,m,+1) * this%rad_grid%c(i,-1) + &
          & sqrt( (jj  ) / (2*jj+1) ) * this%sol%velocity_fn(i+1,j,m,-1) * this%rad_grid%c(i,+1) - &
          & sqrt( (jj+1) / (2*jj+1) ) * this%sol%velocity_fn(i+1,j,m,+1) * this%rad_grid%c(i,+1)

  end function vr_fn

  pure function qr_jm_fn(this, i) result(qr)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: qr(this%jms)
    
    qr = ervs_fn( this%jmax, this%sol%flux_jml_fn(i) )
    
  end function qr_jm_fn
  
  pure function vr_jm_fn(this, i) result(vr)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: vr(this%jms)
    
    vr = ervs_fn( this%jmax, this%rad_grid%c(i,-1) * this%sol%velocity_jml_fn(i  ) + &
                &            this%rad_grid%c(i,+1) * this%sol%velocity_jml_fn(i+1)   )
    
  end function vr_jm_fn

  pure function dv_dr_rrjml_fn(this, i, v) result(dv)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl),       intent(in) :: v(:)
    complex(kind=dbl)                   :: dv(this%jmv)
    
    dv = this%rad_grid%drr(i,-1) * this%sol%velocity_jml_fn(i-1) + &
       & this%rad_grid%drr(i, 0) * v                             + &
       & this%rad_grid%drr(i,+1) * this%sol%velocity_jml_fn(i+1)
    
  end function dv_dr_rrjml_fn

  pure function mgradT_rrjml_fn(this, i) result(gradT)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: gradT(this%jmv)

    gradT = this%rad_grid%cc(i,-1) * this%sol%flux_jml_fn(i-1) / this%lambda_fn(i-1) + &
          & this%rad_grid%cc(i,+1) * this%sol%flux_jml_fn(i  ) / this%lambda_fn(i)

  end function mgradT_rrjml_fn

  subroutine dv_dr_rr_jml_sub(this, i, v, dv)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: i
    complex(kind=dbl),       intent(out) :: dv(:), v(:)
    integer                              :: ijml
    real(kind=dbl)                       :: fac1, fac2, fac3
    complex(kind=dbl), allocatable       :: v1(:), v2(:)
    
    allocate( v1(this%jmv), v2(this%jmv) ) ; call this%sol%velocity_jml_many_sub(i-1, v1, v, v2)

    fac1 = this%rad_grid%drr(i,-1)
    fac2 = this%rad_grid%drr(i, 0)
    fac3 = this%rad_grid%drr(i,+1)

    do concurrent ( ijml = 1:this%jmv )
      dv(ijml) = fac1 * v1(ijml) + fac2 * v(ijml) + fac3 * v2(ijml)
    end do

    deallocate( v1, v2 )
    
  end subroutine dv_dr_rr_jml_sub

  subroutine mgradT_rr_jml_sub(this, i, T, gradT)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: i
    complex(kind=dbl),       intent(out) :: T(:), gradT(:)
    integer                              :: ijml
    real(kind=dbl)                       :: fac1, fac2
    complex(kind=dbl), allocatable       :: flux1(:), flux2(:)
    
    allocate( flux1(this%jmv), flux2(this%jmv) ) ; call this%sol%flux_jml_many_sub(i-1, T, flux1, flux2)
    
    fac1 = this%rad_grid%cc(i,-1) / this%lambda_fn(i-1)
    fac2 = this%rad_grid%cc(i,+1) / this%lambda_fn(i)
    
    do concurrent ( ijml = 1:this%jmv )
      gradT(ijml) = fac1 * flux1(ijml) + fac2 * flux2(ijml)
    end do
    
    deallocate( flux1, flux2 )
    
  end subroutine mgradT_rr_jml_sub
  
end submodule Variables