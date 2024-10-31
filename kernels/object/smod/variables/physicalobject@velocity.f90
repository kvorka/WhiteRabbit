submodule (physicalobject) velocity
  implicit none ; contains
  
  module procedure v_r_fn
    
    v_r_fn = this%rad_grid%c(ir,-1) * this%sol%velocity_fn(ir  ,il,ijm) + &
           & this%rad_grid%c(ir,+1) * this%sol%velocity_fn(ir+1,il,ijm)
    
  end procedure v_r_fn
  
  module procedure v_r_ijml_sub
    integer                        :: ijml
    real(kind=dbl)                 :: cr1, cr2
    complex(kind=dbl), allocatable :: velocity1(:), velocity2(:)
    
    cr1 = this%rad_grid%c(ir,-1)
    cr2 = this%rad_grid%c(ir,+1)
    
    allocate( velocity1(this%jmv), velocity2(this%jmv) )
      
      call this%sol%velocity_jml_many2_sub( ir, velocity1, velocity2 )
      
      do concurrent ( ijml = 1:this%jmv )
        v_r_ijml(ijml) = cr1 * velocity1(ijml) + &
                       & cr2 * velocity2(ijml)
      end do
      
    deallocate( velocity1, velocity2 )
    
  end procedure v_r_ijml_sub
  
  module procedure v_rr_fn
    
    v_rr_fn = this%sol%velocity_fn(ir,il,ijm)
    
  end procedure v_rr_fn
  
  module procedure v_rr_ijml_sub
    
    call this%sol%velocity_jml_many1_sub( ir, v_rr_ijml )
    
  end procedure v_rr_ijml_sub
  
  module procedure vr_r_fn
    real(kind=dbl) :: j, cj1, cj2, cr1, cr2
    
    j = i2r_fn( this%j_indx(ijm) )
    
    cj1 = sqrt((j  )/(2*j+1)) ; cr1 = this%rad_grid%c(ir,-1)
    cj2 = sqrt((j+1)/(2*j+1)) ; cr2 = this%rad_grid%c(ir,+1)
    
    vr_r_fn = cr1 * ( cj1 * this%sol%velocity_fn(ir  ,-1,ijm) - cj2 * this%sol%velocity_fn(ir  ,+1,ijm) ) + &
            & cr2 * ( cj1 * this%sol%velocity_fn(ir+1,-1,ijm) - cj2 * this%sol%velocity_fn(ir+1,+1,ijm) )
    
  end procedure vr_r_fn
  
  module procedure vr_rr_fn
    integer        :: j
    real(kind=dbl) :: cj1, cj2
    
    j = this%j_indx(ijm)
    
    cj1 = sqrt( (j  ) / (2*j+one) )
    cj2 = sqrt( (j+1) / (2*j+one) )
    
    vr_rr_fn = cj1 * this%sol%velocity_fn(ir,-1,ijm) - cj2 * this%sol%velocity_fn(ir,+1,ijm)
    
  end procedure vr_rr_fn
  
  module procedure vr_r_jm_sub
    integer                        :: ij, ijm
    real(kind=dbl)                 :: cj1, cj2, cr1, cr2
    complex(kind=dbl), allocatable :: v1(:), v2(:)
    
    cr1 = this%rad_grid%c(ir,-1)
    cr2 = this%rad_grid%c(ir,+1)
    
    allocate( v1(this%jmv), v2(this%jmv) )
      
      call this%sol%velocity_jml_many2_sub( ir, v1, v2 )
      
      !ij = 0
        vr_jm(1) = czero
      
      do ij = 1, this%jmax
        cj1 = sqrt( (ij  ) / (2*ij+one) )
        cj2 = sqrt( (ij+1) / (2*ij+one) )
        
        do concurrent ( ijm = jm(ij,0):jm(ij,ij) )
          vr_jm(ijm) = cr1 * ( cj1 * v1(3*ijm-4) - cj2 * v1(3*ijm-2) ) + &
                     & cr2 * ( cj1 * v2(3*ijm-4) - cj2 * v2(3*ijm-2) )
        end do
      end do
      
    deallocate( v1, v2 )
    
  end procedure vr_r_jm_sub
  
  module procedure vr_rr_jm_sub
    integer                        :: ij, ijm
    real(kind=dbl)                 :: cj1, cj2
    complex(kind=dbl), allocatable :: v(:)
    
    allocate( v(this%jmv) )
      
      call this%sol%velocity_jml_many1_sub( ir, v )
      
      !ij = 0
        vr_jm(1) = czero
        
      do ij = 1, this%jmax
        cj1 = sqrt( (ij  ) / (2*ij+one) )
        cj2 = sqrt( (ij+1) / (2*ij+one) )
        
        do concurrent ( ijm = jm(ij,0):jm(ij,ij) )
          vr_jm(ijm) = cj1 * v(3*ijm-4) - cj2 * v(3*ijm-2)
        end do
      end do
      
    deallocate( v )
    
  end procedure vr_rr_jm_sub
  
  module procedure dv_dr_rr_jml_sub
    integer                        :: ijml
    real(kind=dbl)                 :: fac1, fac2, fac3
    complex(kind=dbl), allocatable :: v1(:), v2(:)
    
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
    
    allocate( v1(this%jmv), v2(this%jmv) )
      
      call this%sol%velocity_jml_many3_sub(ir-1, v1, v, v2)
      
      do concurrent ( ijml = 1:this%jmv )
        dv(ijml) = fac1 * v1(ijml) + fac2 * v(ijml) + fac3 * v2(ijml)
      end do
      
    deallocate( v1, v2 )
    
  end procedure dv_dr_rr_jml_sub
  
  module procedure curlv_rr_jml_sub
    integer                        :: ij, im, ijml
    real(kind=dbl)                 :: cjr1, cjr2, cjr3, cjr4, crr
    complex(kind=dbl)              :: cj1, cj2
    complex(kind=dbl), allocatable :: dv(:)
    
    crr = 1 / this%rad_grid%rr(ir)
    
    allocate( dv(this%jmv) )
    
    call this%dv_dr_rr_jml_sub( ir, v, dv )
    
    !ij = 0
      !im = 0
        curlv(1) = czero
    
    do ij = 1, this%jmax
      cj1 = sqrt( (ij+1) / (2*ij+one) ) * cunit
      cj2 = sqrt( (ij  ) / (2*ij+one) ) * cunit
      
      cjr1 = (ij-1) * crr
      cjr2 = (ij  ) * crr
      cjr3 = (ij+1) * crr
      cjr4 = (ij+2) * crr
      
      do im = 0, this%jmax
        ijml = 3*(ij*(ij+1)/2+im)-1
        
        curlv(ijml  ) = cj1 * ( dv(ijml+1) + cjr3 * v(ijml+1) )
        curlv(ijml+1) = cj1 * ( dv(ijml  ) - cjr1 * v(ijml  ) ) + &
                      & cj2 * ( dv(ijml+2) + cjr4 * v(ijml+2) )
        curlv(ijml+2) = cj2 * ( dv(ijml+1) - cjr2 * v(ijml+1) )
      end do
    end do
    
    deallocate( dv )
    
  end procedure curlv_rr_jml_sub
  
end submodule velocity