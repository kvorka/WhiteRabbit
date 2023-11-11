submodule(PhysicalObject) Variables
  implicit none

  contains

  pure complex(kind=dbl) function htide_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    htide_fn = this%Ds/this%Ra * ( this%rad_grid%cc(ir,-1) * this%htide(ir-1,ijm) + &
                                 & this%rad_grid%cc(ir,+1) * this%htide(ir  ,ijm)   ) / this%cp_fn(ir)
    
  end function htide_fn

  pure complex(kind=dbl) function vr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    real(kind=dbl)                      :: j, cj1, cj2, cr1, cr2
    
    j = real(this%j_indx(ijm), kind=dbl)
    
    cj1 = sqrt( (j  ) / (2*j+1) ) ; cr1 = this%rad_grid%c(ir,-1)
    cj2 = sqrt( (j+1) / (2*j+1) ) ; cr2 = this%rad_grid%c(ir,+1)
    
    vr_fn = cr1 * ( cj1 * this%sol%velocity_fn(ir  ,-1,ijm) - cj2 * this%sol%velocity_fn(ir  ,+1,ijm) ) + &
          & cr2 * ( cj1 * this%sol%velocity_fn(ir+1,-1,ijm) - cj2 * this%sol%velocity_fn(ir+1,+1,ijm) )
    
  end function vr_fn

  pure complex(kind=dbl) function qr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    real(kind=dbl)                      :: j
    
    j = real(this%j_indx(ijm), kind=dbl)
    
    qr_fn = sqrt( j / (2*j+1) ) * this%sol%flux_fn(ir,-1,ijm) - sqrt( (j+1) / (2*j+1) ) * this%sol%flux_fn(ir,+1,ijm)
    
  end function qr_fn

  pure function qr_jm_fn(this, ir) result(qr)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    complex(kind=dbl),      allocatable :: qr(:)
    
    allocate( qr(this%jms) ) ; qr = ervs_fn( this%jmax, this%sol%flux_jml_fn(ir) )
    
  end function qr_jm_fn
  
  pure function vr_jm_fn(this, ir) result(vr)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    complex(kind=dbl),      allocatable :: vr(:)
    
    allocate( vr(this%jms) ) ; vr = ervs_fn( this%jmax, this%rad_grid%c(ir,-1) * this%sol%velocity_jml_fn(ir  ) + &
                                  &                     this%rad_grid%c(ir,+1) * this%sol%velocity_jml_fn(ir+1)   )
    
  end function vr_jm_fn

  pure function dv_dr_rrjml_fn(this, ir, v) result(dv)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    complex(kind=dbl),       intent(in) :: v(:)
    complex(kind=dbl),      allocatable :: dv(:)
    
    allocate( dv(this%jmv) ) ; dv = this%rad_grid%drr(ir,-1) * this%sol%velocity_jml_fn(ir-1) + &
                                  & this%rad_grid%drr(ir, 0) * v                              + &
                                  & this%rad_grid%drr(ir,+1) * this%sol%velocity_jml_fn(ir+1)
    
  end function dv_dr_rrjml_fn

  pure function mgradT_rrjml_fn(this, ir) result(gradT)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir
    complex(kind=dbl),      allocatable :: gradT(:)
    
    allocate( gradT(this%jmv) ) ; gradT = this%rad_grid%cc(ir,-1) * this%sol%flux_jml_fn(ir-1) / this%lambda_fn(ir-1) + &
                                        & this%rad_grid%cc(ir,+1) * this%sol%flux_jml_fn(ir  ) / this%lambda_fn(ir)
    
  end function mgradT_rrjml_fn

  pure subroutine dv_dr_rr_jml_sub(this, ir, v, dv)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: dv(:), v(:)
    integer                              :: ijml
    real(kind=dbl)                       :: fac1, fac2, fac3
    complex(kind=dbl), allocatable       :: v1(:), v2(:)
    
    fac1 = this%rad_grid%drr(ir,-1)
    fac2 = this%rad_grid%drr(ir, 0)
    fac3 = this%rad_grid%drr(ir,+1)
    
    allocate( v1(this%jmv), v2(this%jmv) )
      
      call this%sol%velocity_jml_many_sub(ir-1, v1, v, v2)
      
      do concurrent ( ijml = 1:this%jmv )
        dv(ijml) = fac1 * v1(ijml) + fac2 * v(ijml) + fac3 * v2(ijml)
      end do
      
    deallocate( v1, v2 )
    
  end subroutine dv_dr_rr_jml_sub

  pure subroutine mgradT_rr_jml_sub(this, ir, T, gradT)
    class(T_physicalObject),     intent(in)  :: this
    integer,                     intent(in)  :: ir
    complex(kind=dbl), optional, intent(out) :: T(:)
    complex(kind=dbl),           intent(out) :: gradT(:)
    integer                                  :: ijml
    real(kind=dbl)                           :: fac1, fac2
    complex(kind=dbl), allocatable           :: flux1(:), flux2(:)
    
    fac1 = this%rad_grid%cc(ir,-1) / this%lambda_fn(ir-1)
    fac2 = this%rad_grid%cc(ir,+1) / this%lambda_fn(ir  )
    
    allocate( flux1(this%jmv), flux2(this%jmv) )
    
    if ( present(T) ) then
      call this%sol%flux_jml_many_sub(ir=ir-1, temp2=T, flux1=flux1, flux2=flux2)
    else
      call this%sol%flux_jml_many_sub(ir=ir-1, flux1=flux1, flux2=flux2)
    end if
    
    do concurrent ( ijml = 1:this%jmv )
      gradT(ijml) = fac1 * flux1(ijml) + fac2 * flux2(ijml)
    end do
    
    deallocate( flux1, flux2 )
    
  end subroutine mgradT_rr_jml_sub
  
end submodule Variables