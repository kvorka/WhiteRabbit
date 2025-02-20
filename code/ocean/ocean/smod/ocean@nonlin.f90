submodule (ocean) nonlin
  implicit none; contains
  
  module procedure coriolis_ocean_sub
    integer                        :: ir, ijm
    complex(kind=dbl), allocatable :: v(:), nlm(:,:)
    
    !$omp parallel private (v, nlm, ijm)
    allocate( v(this%jmv), nlm(3,this%jms) )
    
    !$omp do
    do ir = 2, this%nd
      call this%v_rr_ijml_sub(ir, v)
      
      call zero_carray_sub( 3*this%jms, nlm(1,1) )
      call this%coriolis_rr_jml_sub(v, nlm)
      
      do concurrent ( ijm = 1:this%jms )
        this%nsph1(ijm,ir) = nlm(1,ijm)
        this%ntorr(ijm,ir) = nlm(2,ijm)
        this%nsph2(ijm,ir) = nlm(3,ijm)
      end do
    end do
    !$omp end do
    
    deallocate( v, nlm )
    !$omp end parallel
    
  end procedure coriolis_ocean_sub
  
  module procedure coriolis_vgradv_ocean_sub
    integer                        :: ir, ijm
    real(kind=dbl)                 :: fac
    complex(kind=dbl), allocatable :: v(:), dv(:), nlm(:,:), coriolis(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv), nlm(3,this%jms), coriolis(3,this%jms) )
    
    !$omp parallel private (v, dv, ijm, fac, nlm, coriolis)
    allocate( v(this%jmv) , dv(this%jmv), nlm(3,this%jms), coriolis(3,this%jms) )

    !$omp do
    do ir = 2, this%nd
      call this%dv_dr_rr_jml_sub(ir, v, dv)
      
      call this%lat_grid%vcvgv_sub(this%rad_grid%rr(ir), dv, v, nlm)
      
      select case (this%scaling)
        case ('basics')
          fac = 1 / this%Pr
          
          do concurrent ( ijm = 1:this%jms )
            nlm(1,ijm) = nlm(1,ijm) * fac
            nlm(2,ijm) = nlm(2,ijm) * fac
            nlm(3,ijm) = nlm(3,ijm) * fac
          end do
      end select
      
      call this%coriolis_rr_jml_sub(v, coriolis)
      
      do concurrent ( ijm = 1:this%jms )
        this%nsph1(ijm,ir) = nlm(1,ijm)+coriolis(1,ijm)
        this%ntorr(ijm,ir) = nlm(2,ijm)+coriolis(2,ijm)
        this%nsph2(ijm,ir) = nlm(3,ijm)+coriolis(3,ijm)
      end do
    end do
    !$omp end do

    deallocate( v , dv, nlm, coriolis )
    !$omp end parallel
    
  end procedure coriolis_vgradv_ocean_sub
  
  module procedure fullnl_ocean_sub
    integer                        :: ir, ijm
    real(kind=dbl)                 :: fac
    complex(kind=dbl), allocatable :: v(:), dv(:), T(:), gradT(:), nlm(:,:), buoy(:,:), coriolis(:,:)
    
    !$omp parallel private (v, dv, T, gradT, ijm, fac, nlm, coriolis, buoy)
    allocate( v(this%jmv), dv(this%jmv), T(this%jms), gradT(this%jmv), &
            & nlm(4,this%jms), buoy(2,this%jms), coriolis(3,this%jms)  )
    
    !$omp do
    do ir = 2, this%nd
      call this%dv_dr_rr_jml_sub(ir, v, dv)
      call this%gradT_rr_ijml_sub(ir, T, gradT, -1)
      
      call this%lat_grid%vcvv_vcvgv_sub(this%rad_grid%rr(ir), gradT, dv, v, nlm)
      
      select case (this%scaling)
        case ('basics')
          fac = 1 / this%Pr
          
          do concurrent ( ijm = 1:this%jms )
            nlm(2,ijm) = nlm(2,ijm) * fac
            nlm(3,ijm) = nlm(3,ijm) * fac
            nlm(4,ijm) = nlm(4,ijm) * fac
          end do
      end select
      
      call this%coriolis_rr_jml_sub(v, coriolis)
      call this%buoy_rr_jml_sub(ir, T, buoy)
      
      do concurrent ( ijm = 1:this%jms )
        this%ntemp(ijm,ir) = nlm(1,ijm)
        this%nsph1(ijm,ir) = nlm(2,ijm)+buoy(1,ijm)+coriolis(1,ijm)
        this%ntorr(ijm,ir) = nlm(3,ijm)+            coriolis(2,ijm)
        this%nsph2(ijm,ir) = nlm(4,ijm)+buoy(2,ijm)+coriolis(3,ijm)
      end do
    end do
    !$omp end do
    
    deallocate( v , dv, T , gradT, nlm, buoy, coriolis )
    !$omp end parallel
    
  end procedure fullnl_ocean_sub
  
end submodule nonlin
