submodule (ocean) nonlin
  implicit none; contains
  
  module subroutine coriolis_ocean_sub(this)
    class(T_ocean),    intent(inout) :: this
    integer                          :: ir, ijm
    complex(kind=dbl), allocatable   :: v(:), nlm(:,:)
    
    allocate( v(this%jmv), nlm(3,this%jms) )
    
    !$omp parallel do private (v, nlm, ijm)
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
    !$omp end parallel do
    
    deallocate( v, nlm )
    
  end subroutine coriolis_ocean_sub
  
  module subroutine coriolis_vgradv_ocean_sub(this)
    class(T_ocean),    intent(inout) :: this
    integer                          :: ir, ijm
    real(kind=dbl)                   :: fac
    complex(kind=dbl), allocatable   :: v(:), dv(:), nlm(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv), nlm(3,this%jms) )
    
    !$omp parallel do private (v, dv, ijm, fac, nlm)
    do ir = 2, this%nd
      call this%dv_dr_rr_jml_sub(ir, v, dv)
      
      call this%lat_grid%vcvgv_sub(this%rad_grid%rr(ir), dv, v, nlm)
      
      select case (this%scaling)
        case ('basics')
          fac = 1 / this%Pr
          
          do concurrent ( ijm = 1:this%jms )
            nlm(2,ijm) = nlm(2,ijm) * fac
            nlm(3,ijm) = nlm(3,ijm) * fac
            nlm(4,ijm) = nlm(4,ijm) * fac
          end do
      end select
      
      call this%coriolis_rr_jml_sub(v, nlm)
      
      do concurrent ( ijm = 1:this%jms )
        this%nsph1(ijm,ir) = nlm(1,ijm)
        this%ntorr(ijm,ir) = nlm(2,ijm)
        this%nsph2(ijm,ir) = nlm(3,ijm)
      end do
    end do
    !$omp end parallel do
    
    deallocate( v , dv, nlm )
    
  end subroutine coriolis_vgradv_ocean_sub
  
  module subroutine fullnl_ocean_sub(this)
    class(T_ocean),    intent(inout) :: this
    integer                          :: ir, ijm
    real(kind=dbl)                   :: fac
    complex(kind=dbl), allocatable   :: v(:), dv(:), T(:), gradT(:), nlm(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv),    &
            & T(this%jms) , gradT(this%jmv), &
            & nlm(4,this%jms) )
    
    !$omp parallel do private (v, dv, T, gradT, ijm, fac, nlm)
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
      
      call this%coriolis_rr_jml_sub(v, nlm(2:4,:))
      call this%buoy_rr_jml_sub(ir, T, nlm(2:4,:))
      
      do concurrent ( ijm = 1:this%jms )
        this%ntemp(ijm,ir) = nlm(1,ijm)
        this%nsph1(ijm,ir) = nlm(2,ijm)
        this%ntorr(ijm,ir) = nlm(3,ijm)
        this%nsph2(ijm,ir) = nlm(4,ijm)
      end do
    end do
    !$omp end parallel do
    
    deallocate( v , dv, T , gradT, nlm )
    
  end subroutine fullnl_ocean_sub
  
end submodule nonlin