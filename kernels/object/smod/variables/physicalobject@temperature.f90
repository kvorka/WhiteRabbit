submodule (physicalobject) temperature
  implicit none ; contains
  
  module procedure temp_r_fn
    
    temp_r_fn = this%rad_grid%c(ir,-1) * this%sol%temp_fn(ir  ,ijm) + &
              & this%rad_grid%c(ir,+1) * this%sol%temp_fn(ir+1,ijm)
    
  end procedure temp_r_fn
  
  module procedure temp_r_ijm_sub
    integer                        :: ijm
    real(kind=dbl)                 :: cr1, cr2
    complex(kind=dbl), allocatable :: temp1(:), temp2(:)
    
    cr1 = this%rad_grid%c(ir,-1)
    cr2 = this%rad_grid%c(ir,+1)
    
    allocate( temp1(this%jms), temp2(this%jms) )
    
      call this%sol%temp_jm_many2_sub(ir, temp1, temp2)
      
      do concurrent ( ijm = 1:this%jms )
        temp_r_ijm(ijm) = cr1 * temp1(ijm) + &
                        & cr2 * temp2(ijm)
      end do
    
    deallocate( temp1, temp2 )
    
  end procedure temp_r_ijm_sub
  
  module procedure temp_ir_jm_sub
    integer                        :: ir
    complex(kind=dbl), allocatable :: temp1(:)
    
    allocate( temp1(this%nd+1) )
      
      call this%sol%temp_rr_many1_sub( ijm, temp1 )
      
      do concurrent ( ir = 1:this%nd )
        temp_ir_jm(ijm) = this%rad_grid%c(ir,-1) * temp1(ir  ) + &
                        & this%rad_grid%c(ir,+1) * temp1(ir+1)
      end do
    
    deallocate( temp1 )
    
  end procedure temp_ir_jm_sub
  
  module procedure dT_dr_r_fn
    real(kind=dbl) :: fac1, fac2, fac3, fac4
    
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
    
  end procedure dT_dr_r_fn
  
  module procedure dT_dr_r_ijm_sub
    integer                        :: ijm
    real(kind=dbl)                 :: fac1, fac2, fac3, fac4
    complex(kind=dbl), allocatable :: temp1(:), temp2(:), temp3(:), temp4(:)
    
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
    
  end procedure dT_dr_r_ijm_sub
  
  module procedure gradT_r_ijml_sub
    integer                        :: ij, im, ijm, ijml
    real(kind=dbl)                 :: cj1, cj2, cjr1, cjr2
    complex(kind=dbl), allocatable :: dT_dr(:), T(:)
    
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
    
  end procedure gradT_r_ijml_sub
  
  module procedure temp_rr_fn
    
    temp_rr_fn = this%sol%temp_fn(ir,ijm)
    
  end procedure temp_rr_fn
  
  module procedure temp_rr_ijm_sub
    
    call this%sol%temp_jm_many1_sub( ir, temp_rr_ijm )
    
  end procedure temp_rr_ijm_sub
  
  module procedure temp_irr_jm_sub
    
    call this%sol%temp_rr_many1_sub( ijm, temp_irr_jm )
    
  end procedure temp_irr_jm_sub
  
  module procedure dT_dr_rr_fn
    real(kind=dbl) :: fac1, fac2, fac3
    
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
    
  end procedure dT_dr_rr_fn
  
  module procedure dT_dr_rr_ijm_sub
    integer                        :: ijm
    real(kind=dbl)                 :: fac1, fac2, fac3
    complex(kind=dbl), allocatable :: temp1(:), temp3(:)
    
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
    
  end procedure dT_dr_rr_ijm_sub
  
  module procedure gradT_rr_ijml_sub
    integer                        :: ij, im, ijm, ijml
    real(kind=dbl)                 :: cj1, cj2, cjr1, cjr2
    complex(kind=dbl), allocatable :: dT_dr(:)
    
    allocate( dT_dr(this%jms) )
      
      call this%dT_dr_rr_ijm_sub( ir, T, dT_dr )
      
      ij = 0
        im = 0
          gradT(1) = -sgn * dT_dr(1)
      
      do ij = 1, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) ) * sgn
        cj2 = -sqrt( (ij+1) / (2*ij+one) ) * sgn
        
        cjr1 = +(ij+1) / this%rad_grid%rr(ir)
        cjr2 = -(ij  ) / this%rad_grid%rr(ir)
        
        do im = 0, ij
          ijm  = ij*(ij+1)/2+im+1
          ijml = 3*(ijm-1)-1
          
          gradT(ijml  ) = cj1 * ( dT_dr(ijm) + cjr1 * T(ijm) )
          gradT(ijml+1) = czero
          gradT(ijml+2) = cj2 * ( dT_dr(ijm) + cjr2 * T(ijm) )
        end do
      end do
      
    deallocate( dT_dr )
    
  end procedure gradT_rr_ijml_sub
  
end submodule temperature