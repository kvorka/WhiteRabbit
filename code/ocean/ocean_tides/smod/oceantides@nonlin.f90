submodule (oceantides) nonlin
  implicit none; contains
  
  module procedure coriolis_oceanTides_sub
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
    
  end procedure coriolis_oceanTides_sub
  
  module procedure coriolis_vgradv_oceanTides_sub
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
    
  end procedure coriolis_vgradv_oceanTides_sub
  
end submodule nonlin