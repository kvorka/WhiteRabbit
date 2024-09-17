submodule (PhysicalObject) Variables_heatflux
  implicit none; contains
  
  module pure complex(kind=dbl) function qr_r_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    real(kind=dbl)                      :: j
    
    j = i2r_fn( this%j_indx(ijm) )
    
    qr_r_fn = sqrt((j  )/(2*j+1)) * this%sol%flux_fn(ir,-1,ijm) - &
            & sqrt((j+1)/(2*j+1)) * this%sol%flux_fn(ir,+1,ijm)
    
  end function qr_r_fn
  
  module pure subroutine qr_r_ijm_sub(this, ir, qr_r_ijm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: qr_r_ijm(*)
    integer                              :: ijm
    real(kind=dbl)                       :: j, cj1, cj2
    
    do concurrent ( ijm = 1:this%jms )
      j = i2r_fn( this%j_indx(ijm) )
      
      cj1 = +sqrt((j  )/(2*j+1))
      cj2 = -sqrt((j+1)/(2*j+1))
      
      qr_r_ijm(ijm) = cj1 * this%sol%flux_fn(ir,-1,ijm) + &
                    & cj2 * this%sol%flux_fn(ir,+1,ijm)
    end do
    
  end subroutine qr_r_ijm_sub
  
  module pure subroutine qr_ir_jm_sub(this, ijm, qr_ir_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: qr_ir_jm(*)
    integer                              :: ir
    real(kind=dbl)                       :: j, cj1, cj2
    
    j = i2r_fn( this%j_indx(ijm) )
    
    cj1 = +sqrt((j  )/(2*j+1))
    cj2 = -sqrt((j+1)/(2*j+1))
    
    do concurrent ( ir = 1:this%nd )
      qr_ir_jm(ijm) = cj1 * this%sol%flux_fn(ir,-1,ijm) + &
                    & cj2 * this%sol%flux_fn(ir,+1,ijm)
    end do
    
  end subroutine qr_ir_jm_sub
  
  module pure complex(kind=dbl) function qr_rr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    real(kind=dbl)                      :: j, cj1, cj2, cr1, cr2
    
    j = i2r_fn( this%j_indx(ijm) )
    
    cj1 = +sqrt((j  )/(2*j+1)) ; cr1 = this%rad_grid%cc(ir,-1)
    cj2 = -sqrt((j+1)/(2*j+1)) ; cr2 = this%rad_grid%cc(ir,+1)
    
    qr_rr_fn = cr1 * ( cj1 * this%sol%flux_fn(ir-1,-1,ijm) + cj2 * this%sol%flux_fn(ir-1,+1,ijm) ) + &
             & cr2 * ( cj1 * this%sol%flux_fn(ir  ,-1,ijm) + cj2 * this%sol%flux_fn(ir  ,+1,ijm) )
    
  end function qr_rr_fn
  
  module pure subroutine qr_rr_ijm_sub(this, ir, qr_rr_ijm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: qr_rr_ijm(*)
    integer                              :: ijm
    real(kind=dbl)                       :: j, cr1, cr2, cj1, cj2
    
    cr1 = this%rad_grid%cc(ir,-1)
    cr2 = this%rad_grid%cc(ir,+1)
    
    do concurrent ( ijm = 1:this%jms )
      j = i2r_fn( this%j_indx(ijm) )
      
      cj1 = +sqrt((j  )/(2*j+1))
      cj2 = -sqrt((j+1)/(2*j+1))
      
      qr_rr_ijm(ijm) = cr1 * ( cj1 * this%sol%flux_fn(ir-1,-1,ijm) + cj2 * this%sol%flux_fn(ir-1,+1,ijm) ) + &
                     & cr2 * ( cj1 * this%sol%flux_fn(ir  ,-1,ijm) + cj2 * this%sol%flux_fn(ir  ,+1,ijm) )
    end do
    
  end subroutine qr_rr_ijm_sub
  
  module pure subroutine qr_irr_jm_sub(this, ijm, qr_irr_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: qr_irr_jm(2:this%nd)
    integer                              :: ir
    real(kind=dbl)                       :: j, cj1, cj2, cr1, cr2
    
    j = i2r_fn( this%j_indx(ijm) )
    
    cj1 = +sqrt((j  )/(2*j+1))
    cj2 = -sqrt((j+1)/(2*j+1))
    
    do concurrent ( ir = 2:this%nd )
      cr1 = this%rad_grid%cc(ir,-1)
      cr2 = this%rad_grid%cc(ir,+1)
      
      qr_irr_jm(ijm) = cr1 * ( cj1 * this%sol%flux_fn(ir-1,-1,ijm) + cj2 * this%sol%flux_fn(ir-1,+1,ijm) ) + &
                     & cr2 * ( cj1 * this%sol%flux_fn(ir  ,-1,ijm) + cj2 * this%sol%flux_fn(ir  ,+1,ijm) )
    end do
    
  end subroutine qr_irr_jm_sub
  
end submodule Variables_heatflux