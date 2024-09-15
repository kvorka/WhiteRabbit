submodule(PhysicalObject) Variables_temperature
  implicit none ; contains
  
  module pure complex(kind=dbl) function temp_r_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    temp_r_fn = this%rad_grid%c(ir,-1) * this%sol%temp_fn(ir  ,ijm) + &
             & this%rad_grid%c(ir,+1) * this%sol%temp_fn(ir+1,ijm)
    
  end function temp_r_fn
  
  module pure complex(kind=dbl) function temp_rr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    temp_rr_fn = this%sol%temp_fn(ir,ijm)
    
  end function temp_rr_fn
  
  module pure complex(kind=dbl) function qr_r_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    integer                             :: j
    
    j = this%j_indx(ijm)
    
    qr_r_fn = sqrt( (j  ) / (2*j+one) ) * this%sol%flux_fn(ir,-1,ijm) - &
            & sqrt( (j+1) / (2*j+one) ) * this%sol%flux_fn(ir,+1,ijm)
    
  end function qr_r_fn
  
  module pure subroutine mgradT_rr_jml_sub(this, ir, T, gradT)
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
  
end submodule Variables_temperature