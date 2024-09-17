submodule(PhysicalObject) Variables_temperature
  implicit none ; contains
  
  module pure complex(kind=dbl) function temp_r_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    temp_r_fn = this%rad_grid%c(ir,-1) * this%sol%temp_fn(ir  ,ijm) + &
              & this%rad_grid%c(ir,+1) * this%sol%temp_fn(ir+1,ijm)
    
  end function temp_r_fn
  
  module pure subroutine temp_r_ijm_sub(this, ir, temp_r_ijm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: temp_r_ijm(*)
    integer                              :: ijm
    real(kind=dbl)                       :: cr1, cr2
    
    cr1 = this%rad_grid%c(ir,-1)
    cr2 = this%rad_grid%c(ir,+1)
    
    do concurrent ( ijm = 1:this%jms )
      temp_r_ijm(ijm) = cr1 * this%sol%temp_fn(ir  ,ijm) + &
                      & cr2 * this%sol%temp_fn(ir+1,ijm)
    end do
    
  end subroutine temp_r_ijm_sub
  
  module pure subroutine temp_ir_jm_sub(this, ijm, temp_ir_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: temp_ir_jm(*)
    integer                              :: ir
    
    do concurrent ( ir = 1:this%nd )
      temp_ir_jm(ijm) = this%rad_grid%c(ir,-1) * this%sol%temp_fn(ir  ,ijm) + &
                      & this%rad_grid%c(ir,+1) * this%sol%temp_fn(ir+1,ijm)
    end do
    
  end subroutine temp_ir_jm_sub
  
  module pure complex(kind=dbl) function dT_dr_r_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    real(kind=dbl)                      :: fac1, fac2, fac3, fac4
    
    select case (this%grid_type)
      case ('homog')
        fac1 = zero
        fac2 = hd(ir,-1)
        fac3 = hd(ir,+1)
        fac4 = zero
      
      case ('chebv')
        fac1 = d(ir,-2)
        fac2 = d(ir,-1)
        fac3 = d(ir,+1)
        fac4 = d(ir,+2)
        
    end select
    
    if ( (ir > 1) .and. (ir < this%nd) ) then
      dT_dr_r_fn = fac1 * this%sol%temp_fn(ir-1,ijm) + &
                 & fac2 * this%sol%temp_fn(ir  ,ijm) + &
                 & fac3 * this%sol%temp_fn(ir+1,ijm) + &
                 & fac4 * this%sol%temp_fn(ir+2,ijm)
    else if ( ir == 1) then
      dT_dr_r_fn = fac2 * this%sol%temp_fn(ir  ,ijm) + &
                 & fac3 * this%sol%temp_fn(ir+1,ijm) + &
                 & fac4 * this%sol%temp_fn(ir+2,ijm)
    
    else
      dT_dr_r_fn = fac1 * this%sol%temp_fn(ir-1,ijm) + &
                 & fac2 * this%sol%temp_fn(ir  ,ijm) + &
                 & fac3 * this%sol%temp_fn(ir+1,ijm)
    
    end if
    
  end function dT_dr_r_fn
  
  module pure complex(kind=dbl) function temp_rr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    temp_rr_fn = this%sol%temp_fn(ir,ijm)
    
  end function temp_rr_fn
  
  module pure subroutine temp_rr_ijm_sub(this, ir, temp_rr_ijm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: temp_rr_ijm(*)
    integer                              :: ijm
    
    do concurrent ( ijm = 1:this%jms )
      temp_rr_ijm(ijm) = this%sol%temp_fn(ir,ijm)
    end do
    
  end subroutine temp_rr_ijm_sub
  
  module pure subroutine temp_irr_jm_sub(this, ijm, temp_irr_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: temp_irr_jm(*)
    integer                              :: ir
    
    do concurrent ( ir = 1:this%nd+1 )
      temp_irr_jm(ir) = this%sol%temp_fn(ir,ijm)
    end do
    
  end subroutine temp_irr_jm_sub
  
  module pure complex(kind=dbl) function dT_dr_rr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    real(kind=dbl)                      :: fac1, fac2, fac3
    
    select case ( this%grid_type )
      case ('homog')
        fac1 = this%rad_grid%hdrr(ir,-1)
        fac2 = this%rad_grid%hdrr(ir, 0)
        fac3 = this%rad_grid%hdrr(ir,+1)
      
      case('chebv')
        fac1 = this%rad_grid%drr(ir,-1)
        fac2 = this%rad_grid%drr(ir, 0)
        fac3 = this%rad_grid%drr(ir,+1)
    end select
    
    dT_dr_rr_fn = fac1 * this%sol%temp_fn(ir-1,ijm) + &
                & fac2 * this%sol%temp_fn(ir  ,ijm) + &
                & fac3 * this%sol%temp_fn(ir+1,ijm)
    
  end function dT_dr_rr_fn
  
  module pure subroutine dT_dr_rr_jm_sub(this, ir, T, dT)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: T(*), dT(*)
    integer                              :: ijm
    real(kind=dbl)                       :: fac1, fac2, fac3
    complex(kind=dbl), allocatable       :: temp1(:), temp3(:)
    
    select case ( this%grid_type )
      case ('homog')
        fac1 = this%rad_grid%hdrr(ir,-1)
        fac2 = this%rad_grid%hdrr(ir, 0)
        fac3 = this%rad_grid%hdrr(ir,+1)
      
      case('chebv')
        fac1 = this%rad_grid%drr(ir,-1)
        fac2 = this%rad_grid%drr(ir, 0)
        fac3 = this%rad_grid%drr(ir,+1)
    end select
    
    allocate( temp1(this%jms), temp3(this%jms) )
      
      call this%sol%temp_jm_many_sub( ir-1, temp1(1), T(1), temp3(1) )
      
      do concurrent ( ijm = 1:this%jms )
        dT(ijm) = fac1 * temp1(ijm) + &
                & fac2 * T(ijm)     + &
                & fac3 * temp3(ijm)
      end do
    
    deallocate( temp1, temp3 )
    
  end subroutine dT_dr_rr_jm_sub
  
  module pure subroutine mgradT_rr_jml_sub(this, ir, T, gradT)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: T(*), gradT(*)
    integer                              :: ij, im, ijm
    real(kind=dbl)                       :: cj1, cj2
    complex(kind=dbl),       allocatable :: dT_dr(:)
    
    allocate( dT_dr(this%jms) )
      
      call this%dT_dr_rr_jm_sub( ir, T(1), dT_dr(1) )
      
      ij = 0
        im = 0
          gradT(1) = -dT_dr(1)
      
      do ij = 1, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) )
        cj2 = -sqrt( (ij+1) / (2*ij+one) )
        
        do im = 0, ij
          ijm = ij*(ij+1)/2+im+1
          
          gradT(3*(ijm-1)-1) = -cj1 * ( dT_dr(ijm) + (ij+1) / this%rad_grid%rr(ir) * T(ijm) )
          gradT(3*(ijm-1)  ) = -czero
          gradT(3*(ijm+1)+1) = -cj2 * ( dT_dr(ijm) - (ij  ) / this%rad_grid%rr(ir) * T(ijm) )
        end do
      end do
      
    deallocate( dT_dr )
    
  end subroutine mgradT_rr_jml_sub
  
end submodule Variables_temperature