submodule (PhysicalObject) Variables_heatflux
  implicit none; contains
  
  module pure complex(kind=dbl) function q_r_fn(this, ir, il, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, il, ijm
    
    q_r_fn = this%sol%flux_fn(ir, il, ijm)
    
  end function q_r_fn
  
  module pure subroutine q_r_ijml_sub(this, ir, q_r_ijml)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: q_r_ijml(:)
    
    call this%sol%flux_jml_many1_sub( ir, q_r_ijml )
    
  end subroutine q_r_ijml_sub
  
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
    integer                              :: ij, im, ijm, ijml
    real(kind=dbl)                       :: cj1, cj2
    complex(kind=dbl),       allocatable :: flux1(:)
    
    allocate( flux1(this%jmv) )
      
      call this%sol%flux_jml_many1_sub( ir, flux1 )
      
      ij = 0
        ijm  = 1
        ijml = 1
          qr_r_ijm(1) = -flux1(1)
      
      do ij = 1, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) )
        cj2 = -sqrt( (ij+1) / (2*ij+one) )
        
        do im = 0, ij
          ijm  = ij*(ij+1)/2+im+1
          ijml = 3*(ijm-1)-1
          
          qr_r_ijm(ijm) = cj1 * flux1(ijml  ) + &
                        & cj2 * flux1(ijml+2)
        end do
      end do
    
    deallocate( flux1 )
    
  end subroutine qr_r_ijm_sub
  
  module pure subroutine qr_ir_jm_sub(this, ijm, qr_ir_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: qr_ir_jm(*)
    integer                              :: ir
    real(kind=dbl)                       :: j, cj1, cj2
    complex(kind=dbl),       allocatable :: flux1(:,:)
    
    j = i2r_fn( this%j_indx(ijm) )
    
    cj1 = +sqrt((j  )/(2*j+1))
    cj2 = -sqrt((j+1)/(2*j+1))
    
    allocate( flux1(3,this%nd) )
    
      call this%sol%flux_r_many1_sub( ijm, flux1 )
      
      do concurrent ( ir = 1:this%nd )
        qr_ir_jm(ijm) = cj1 * flux1(1,ijm) + &
                      & cj2 * flux1(3,ijm)
      end do
      
    deallocate( flux1 )
    
  end subroutine qr_ir_jm_sub
  
  module pure complex(kind=dbl) function q_rr_fn(this, ir, il, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, il, ijm
    real(kind=dbl)                      :: cr1, cr2
    
    cr1 = this%rad_grid%cc(ir,-1)
    cr2 = this%rad_grid%cc(ir,+1)
    
    q_rr_fn = cr1 * this%sol%flux_fn(ir-1, il, ijm) + &
            & cr2 * this%sol%flux_fn(ir  , il, ijm)
    
  end function q_rr_fn
  
  module pure subroutine q_rr_ijml_sub(this, ir, q_rr_ijml)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: q_rr_ijml(:)
    integer                              :: ijml
    real(kind=dbl)                       :: cr1, cr2
    complex(kind=dbl),       allocatable :: flux1(:), flux2(:)
    
    cr1 = this%rad_grid%cc(ir,-1)
    cr2 = this%rad_grid%cc(ir,+1)
    
    allocate( flux1(this%jmv), flux2(this%jmv) )
      
      call this%sol%flux_jml_many2_sub( ir-1, flux1, flux2 )
      
      do concurrent ( ijml = 1:this%jmv )
        q_rr_ijml(ijml  ) = cr1 * flux1(ijml) + &
                          & cr2 * flux2(ijml)
      end do
    
    deallocate( flux1, flux2 )
    
  end subroutine q_rr_ijml_sub
  
  module pure complex(kind=dbl) function dq_dr_rr_fn(this, ir, il, ijm)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: ir, il, ijm
    real(kind=dbl)                      :: fac1, fac2, fac3, fac4
    
    select case (this%grid_type)
      case ('homog')
        fac1 = zero
        fac2 = this%rad_grid%hdd(ir,-1)
        fac3 = this%rad_grid%hdd(ir,+1)
        fac4 = zero
        
      case ('chebv')
        fac1 = this%rad_grid%dd(ir,-2)
        fac2 = this%rad_grid%dd(ir,-1)
        fac3 = this%rad_grid%dd(ir,+1)
        fac4 = this%rad_grid%dd(ir,+2)
        
    end select
    
    if ( (ir > 3) .and. (ir < this%nd) ) then
      dq_dr_rr_fn = fac1 * this%sol%flux_fn(ir-2,il,ijm) + &
                  & fac2 * this%sol%flux_fn(ir-1,il,ijm) + &
                  & fac3 * this%sol%flux_fn(ir  ,il,ijm) + &
                  & fac4 * this%sol%flux_fn(ir+1,il,ijm)
    
    else if ( ir == 2) then
      dq_dr_rr_fn = fac2 * this%sol%flux_fn(ir-1,il,ijm) + &
                  & fac3 * this%sol%flux_fn(ir  ,il,ijm) + &
                  & fac4 * this%sol%flux_fn(ir+1,il,ijm)
    
    else
      dq_dr_rr_fn = fac1 * this%sol%flux_fn(ir-2,il,ijm) + &
                  & fac2 * this%sol%flux_fn(ir-1,il,ijm) + &
                  & fac3 * this%sol%flux_fn(ir  ,il,ijm)
    
    end if
    
  end function dq_dr_rr_fn
  
  module pure subroutine dq_dr_rr_ijml_sub(this, ir, dq_dr_rr_ijml)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: dq_dr_rr_ijml(:)
    integer                              :: ijml
    real(kind=dbl)                       :: fac1, fac2, fac3, fac4
    complex(kind=dbl),       allocatable :: flux1(:), flux2(:), flux3(:), flux4(:) 
    
    select case (this%grid_type)
      case ('homog')
        fac1 = zero
        fac2 = this%rad_grid%hdd(ir,-1)
        fac3 = this%rad_grid%hdd(ir,+1)
        fac4 = zero
        
      case ('chebv')
        fac1 = this%rad_grid%dd(ir,-2)
        fac2 = this%rad_grid%dd(ir,-1)
        fac3 = this%rad_grid%dd(ir,+1)
        fac4 = this%rad_grid%dd(ir,+2)
        
    end select
    
    if ( (ir > 3) .and. (ir < this%nd) ) then
      
      allocate( flux1(this%jmv), flux2(this%jmv), flux3(this%jmv), flux4(this%jmv) )
        
        call this%sol%flux_jml_many4_sub( ir-2, flux1, flux2, flux3, flux4 )
        
        do concurrent ( ijml = 1:this%jmv )
          dq_dr_rr_ijml = fac1 * flux1(ijml) + &
                        & fac2 * flux2(ijml) + &
                        & fac3 * flux3(ijml) + &
                        & fac4 * flux4(ijml)
        end do
        
      deallocate( flux1, flux2, flux3, flux4 )
      
    else if ( ir == 2) then
      
      allocate( flux2(this%jmv), flux3(this%jmv), flux4(this%jmv) )
        
        call this%sol%flux_jml_many3_sub( ir-1, flux2, flux3, flux4 )
        
        do concurrent ( ijml = 1:this%jmv )
          dq_dr_rr_ijml = fac2 * flux2(ijml) + &
                        & fac3 * flux3(ijml) + &
                        & fac4 * flux4(ijml)
        end do
        
      deallocate( flux2, flux3, flux4 )
      
    else
      
      allocate( flux1(this%jmv), flux2(this%jmv), flux3(this%jmv) )
        
        call this%sol%flux_jml_many3_sub( ir-2, flux1, flux2, flux3 )
        
        do concurrent ( ijml = 1:this%jmv )
          dq_dr_rr_ijml = fac1 * flux1(ijml) + &
                        & fac2 * flux2(ijml) + &
                        & fac3 * flux3(ijml)
        end do
        
      deallocate( flux1, flux2, flux3 )
      
    end if
    
  end subroutine dq_dr_rr_ijml_sub
  
  module pure subroutine divq_rr_ijm_sub(this, ir, divq_rr_ijm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ir
    complex(kind=dbl),       intent(out) :: divq_rr_ijm(:)
    integer                              :: ij, im, ijm, ijml
    real(kind=dbl)                       :: cj1, cj2, cjr1, cjr2
    complex(kind=dbl),       allocatable :: dq_dr(:), q(:)
    
    allocate( dq_dr(this%jmv), q(this%jmv) )
      
      call this%dq_dr_rr_ijml_sub(ir, dq_dr)
      call this%q_rr_ijml_sub(ir, q)
      
      ij = 0
        ijm  = 1
        ijml = 1
          divq_rr_ijm(1) = -( dq_dr(1) + 2 * q(1) / this%rad_grid%rr(ir) )
      
      do ij = 1, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) )
        cj2 = -sqrt( (ij+1) / (2*ij+one) )
        
        cjr1 = -(ij-1) / this%rad_grid%rr(ir)
        cjr2 = +(ij+2) / this%rad_grid%rr(ir)
        
        do im = 0, ij
          ijm  = ij*(ij+1)/2+im+1
          ijml = 3*(ijm-1)-1
          
          divq_rr_ijm(ijm) = cj1 * ( dq_dr(ijml  ) + cjr1 * q(ijml  ) ) + &
                           & cj2 * ( dq_dr(ijml+2) + cjr2 * q(ijml+2) )
        end do
      end do
      
    deallocate( dq_dr, q )
    
  end subroutine divq_rr_ijm_sub
  
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
    integer                              :: ij, im, ijm, ijml
    real(kind=dbl)                       :: cr1, cr2, cj1, cj2
    complex(kind=dbl),       allocatable :: flux1(:), flux2(:)
    
    cr1 = this%rad_grid%cc(ir,-1)
    cr2 = this%rad_grid%cc(ir,+1)
    
    allocate( flux1(this%jmv), flux2(this%jmv) )
    
      call this%sol%flux_jml_many2_sub( ir-1, flux1, flux2 )
      
      do ij = 0, this%jmax
        cj1 = +sqrt( (ij  ) / (2*ij+one) )
        cj2 = -sqrt( (ij+1) / (2*ij+one) )
        
        do im = 0, ij
          ijm  = ij*(ij+1)/2+im+1
          ijml = 3*(ijm-1)-1
          
          qr_rr_ijm(ijm) = cr1 * ( cj1 * flux1(ijml) + cj2 * flux1(ijml+2) ) + &
                         & cr2 * ( cj1 * flux2(ijml) + cj2 * flux2(ijml+2) )
        end do
      end do
    
    deallocate( flux1, flux2 )
    
  end subroutine qr_rr_ijm_sub
  
  module pure subroutine qr_irr_jm_sub(this, ijm, qr_irr_jm)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: ijm
    complex(kind=dbl),       intent(out) :: qr_irr_jm(2:this%nd)
    integer                              :: ir
    real(kind=dbl)                       :: j, cj1, cj2, cr1, cr2
    complex(kind=dbl),       allocatable :: flux1(:,:)
    
    j = i2r_fn( this%j_indx(ijm) )
    
    cj1 = +sqrt((j  )/(2*j+1))
    cj2 = -sqrt((j+1)/(2*j+1))
    
    allocate( flux1(3,this%nd) )
      
      call this%sol%flux_r_many1_sub( ijm, flux1 )
      
      do concurrent ( ir = 2:this%nd )
        cr1 = this%rad_grid%cc(ir,-1)
        cr2 = this%rad_grid%cc(ir,+1)
        
        qr_irr_jm(ijm) = cr1 * ( cj1 * flux1(1,ir-1) + cj2 * flux1(3,ir-1) ) + &
                       & cr2 * ( cj1 * flux1(1,ir  ) + cj2 * flux1(3,ir  ) )
      end do
      
    deallocate( flux1 )
    
  end subroutine qr_irr_jm_sub
  
end submodule Variables_heatflux