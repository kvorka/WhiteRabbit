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
    complex(kind=dbl),       intent(out) :: temp_r_ijm(:)
    integer                              :: ijm
    real(kind=dbl)                       :: cr1, cr2
    complex(kind=dbl),       allocatable :: temp1(:), temp2(:)
    
    cr1 = this%rad_grid%c(ir,-1)
    cr2 = this%rad_grid%c(ir,+1)
    
    allocate( temp1(this%jms), temp2(this%jms) )
    
      call this%sol%temp_jm_many2_sub(ir, temp1, temp2)
      
      do concurrent ( ijm = 1:this%jms )
        temp_r_ijm(ijm) = cr1 * temp1(ijm) + &
                        & cr2 * temp2(ijm)
      end do
    
    deallocate( temp1, temp2 )
    
  end subroutine temp_r_ijm_sub
  
  module pure subroutine temp_ir_jm_sub(this, ijm, temp_ir_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: temp_ir_jm(:)
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
        fac2 = this%rad_grid%hd(ir,-1)
        fac3 = this%rad_grid%hd(ir,+1)
        fac4 = zero
      
      case ('chebv')
        fac1 = this%rad_grid%d(ir,-2)
        fac2 = this%rad_grid%d(ir,-1)
        fac3 = this%rad_grid%d(ir,+1)
        fac4 = this%rad_grid%d(ir,+2)
        
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
  
  module pure subroutine dT_dr_r_ijm_sub(this, ir, dT_dr_r)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: dT_dr_r(:)
    integer                              :: ijm
    real(kind=dbl)                       :: fac1, fac2, fac3, fac4
    complex(kind=dbl),       allocatable :: temp1(:), temp2(:), temp3(:), temp4(:)
    
    select case (this%grid_type)
      case ('homog')
        fac1 = zero
        fac2 = this%rad_grid%hd(ir,-1)
        fac3 = this%rad_grid%hd(ir,+1)
        fac4 = zero
      
      case ('chebv')
        fac1 = this%rad_grid%d(ir,-2)
        fac2 = this%rad_grid%d(ir,-1)
        fac3 = this%rad_grid%d(ir,+1)
        fac4 = this%rad_grid%d(ir,+2)
        
    end select
    
    if ( (ir > 1) .and. (ir < this%nd) ) then
      allocate( temp1(this%jms), temp2(this%jms), temp3(this%jms), temp4(this%jms) )
      
        call this%sol%temp_jm_many4_sub( ir-1, temp1, temp2, temp3, temp4 )
        
        do concurrent ( ijm = 1:this%jms )
          dT_dr_r(ijm) = fac1 * temp1(ijm) + &
                       & fac2 * temp2(ijm) + &
                       & fac3 * temp3(ijm) + &
                       & fac4 * temp4(ijm)
        end do
      
      deallocate( temp1, temp2, temp3, temp4 )
      
    else if ( ir == 1) then
      allocate( temp2(this%jms), temp3(this%jms), temp4(this%jms) )
      
        call this%sol%temp_jm_many3_sub( ir, temp2, temp3, temp4 )
        
        do concurrent ( ijm = 1:this%jms )
          dT_dr_r(ijm) = fac2 * temp2(ijm) + &
                       & fac3 * temp3(ijm) + &
                       & fac4 * temp4(ijm)
        end do
      
      deallocate( temp2, temp3, temp4 )
    
    else
      allocate( temp1(this%jms), temp2(this%jms), temp3(this%jms) )
      
        call this%sol%temp_jm_many3_sub( ir-1, temp1, temp2, temp3 )
        
        do concurrent ( ijm = 1:this%jms )
          dT_dr_r(ijm) = fac1 * temp1(ijm) + &
                       & fac2 * temp2(ijm) + &
                       & fac3 * temp3(ijm)
        end do
      
      deallocate( temp1, temp2, temp3 )
    
    end if
    
  end subroutine dT_dr_r_ijm_sub
  
  module pure subroutine gradT_r_ijml_sub(this, ir, gradT, sgn)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir, sgn
    complex(kind=dbl),       intent(out) :: gradT(:)
    integer                              :: ij, im, ijm, ijml
    real(kind=dbl)                       :: cj1, cj2, cjr1, cjr2
    complex(kind=dbl),       allocatable :: dT_dr(:), T(:)
    
    allocate( dT_dr(this%jms) ) ; call this%dT_dr_r_ijm_sub( ir, dT_dr )
    allocate( T(this%jms) )     ; call this%temp_r_ijm_sub( ir, T )
      
      ij = 0
        im = 0
          gradT(1) = -dT_dr(1)
      
      do ij = 1, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) ) * sgn
        cj2 = -sqrt( (ij+1) / (2*ij+one) ) * sgn
        
        cjr1 = +(ij+1) / this%rad_grid%r(ir)
        cjr2 = -(ij  ) / this%rad_grid%r(ir)
        
        do im = 0, ij
          ijm  = ij*(ij+1)/2+im+1
          ijml = 3*(ijm-1)-1
          
          gradT(ijml  ) = cj1 * ( dT_dr(ijm) + cjr1 * T(ijm) )
          gradT(ijml+1) = czero
          gradT(ijml+2) = cj2 * ( dT_dr(ijm) + cjr2 * T(ijm) )
        end do
      end do
      
    deallocate( dT_dr, T )
    
  end subroutine gradT_r_ijml_sub
  
  module pure complex(kind=dbl) function temp_rr_fn(this, ir, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, ijm
    
    temp_rr_fn = this%sol%temp_fn(ir,ijm)
    
  end function temp_rr_fn
  
  module pure subroutine temp_rr_ijm_sub(this, ir, temp_rr_ijm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: temp_rr_ijm(:)
    
    call this%sol%temp_jm_many1_sub( ir, temp_rr_ijm )
    
  end subroutine temp_rr_ijm_sub
  
  module pure subroutine temp_irr_jm_sub(this, ijm, temp_irr_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: temp_irr_jm(:)
    
    call this%sol%temp_rr_many1_sub( ijm, temp_irr_jm )
    
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
  
  module pure subroutine dT_dr_rr_ijm_sub(this, ir, T, dT)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: T(:), dT(:)
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
      
      call this%sol%temp_jm_many3_sub( ir-1, temp1, T, temp3 )
      
      do concurrent ( ijm = 1:this%jms )
        dT(ijm) = fac1 * temp1(ijm) + &
                & fac2 * T(ijm)     + &
                & fac3 * temp3(ijm)
      end do
    
    deallocate( temp1, temp3 )
    
  end subroutine dT_dr_rr_ijm_sub
  
  module pure subroutine gradT_rr_ijml_sub(this, ir, T, gradT, sgn)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir, sgn
    complex(kind=dbl),       intent(out) :: T(:), gradT(:)
    integer                              :: ij, im, ijm
    real(kind=dbl)                       :: cj1, cj2, cjr1, cjr2
    complex(kind=dbl),       allocatable :: dT_dr(:)
    
    allocate( dT_dr(this%jms) ) ; call this%dT_dr_rr_ijm_sub( ir, T, dT_dr )
      
      ij = 0
        im = 0
          gradT(1) = -dT_dr(1)
      
      do ij = 1, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) ) * sgn
        cj2 = -sqrt( (ij+1) / (2*ij+one) ) * sgn
        
        cjr1 = +(ij+1) / this%rad_grid%rr(ir)
        cjr2 = -(ij  ) / this%rad_grid%rr(ir)
        
        do im = 0, ij
          ijm = ij*(ij+1)/2+im+1
          
          gradT(3*(ijm-1)-1) = cj1 * ( dT_dr(ijm) + cjr1 * T(ijm) )
          gradT(3*(ijm-1)  ) = czero
          gradT(3*(ijm-1)+1) = cj2 * ( dT_dr(ijm) + cjr2 * T(ijm) )
        end do
      end do
      
    deallocate( dT_dr )
    
  end subroutine gradT_rr_ijml_sub
  
end submodule Variables_temperature