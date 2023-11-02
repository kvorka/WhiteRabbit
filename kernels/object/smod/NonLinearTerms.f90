submodule (PhysicalObject) NonLinearTerms
  contains
  
  function vgradT_fn(this, i) result(vgradT)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl),      allocatable :: v(:), gradT(:), vgradT(:)
    
    allocate( v(this%jmv), gradT(this%jmv), vgradT(this%jms) )
      
      v = this%sol%velocity_jml_fn(i) ; gradT = -this%mgradT_rrjml_fn(i) ; call this%lat_grid%vcvv_sub( v, gradT, vgradT )
      
    deallocate( v, gradT )
    
  end function vgradT_fn
  
  subroutine coriolis_sub(this, i)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: i
    integer                                :: ijm
    complex(kind=dbl),       allocatable   :: v(:), nlm(:,:)
    
    allocate( v(this%jmv) ) ; call this%sol%velocity_jml_sub(i, v)
    
    allocate( nlm(3,this%jms) ) ; nlm = czero
    
    call this%coriolis_rr_jml_sub(v, nlm)
    
    deallocate( v )
    
    do concurrent ( ijm = 1:this%jms )
      this%nsph1(ijm,i) = nlm(1,ijm)
      this%ntorr(ijm,i) = nlm(2,ijm)
      this%nsph2(ijm,i) = nlm(3,ijm)
    end do
    
    deallocate( nlm )
    
  end subroutine coriolis_sub
  
  subroutine coriolis_vgradv_sub(this, i)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: i
    integer                                :: ijm, i1
    complex(kind=dbl),       allocatable   :: v(:), dv(:), nlm(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv) ) ; call this%dv_dr_rr_jml_sub(i, v, dv)
    
    allocate( nlm(3,this%jms) ) ; call this%lat_grid%vcsv_vcvgv_sub(this%rad_grid%rr(i), dv, v, nlm)
    
    deallocate( dv )
    
    select case (this%scaling)
      case ('basics')
        do concurrent ( ijm = 1:this%jms, i1 = 2:4 )
          nlm(i1,ijm) = nlm(i1,ijm) / this%Pr
        end do
    end select
    
    call this%coriolis_rr_jml_sub(v, nlm)
    
    deallocate( v )
    
    do concurrent ( ijm = 1:this%jms )
      this%nsph1(ijm,i) = nlm(1,ijm)
      this%ntorr(ijm,i) = nlm(2,ijm)
      this%nsph2(ijm,i) = nlm(3,ijm)
    end do
    
    deallocate( nlm )
    
  end subroutine coriolis_vgradv_sub
  
  subroutine fullnl_sub(this, i)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: i
    integer                                :: ijm, i1
    complex(kind=dbl),       allocatable   :: v(:), dv(:), T(:), dT(:), nlm(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv) ) ; call this%dv_dr_rr_jml_sub(i, v, dv)
    allocate( T(this%jms) , dT(this%jmv) ) ; call this%mgradT_rr_jml_sub(i, T, dT)
    
    allocate( nlm(4,this%jms) ) ; call this%lat_grid%vcsv_vcvv_vcvgv_sub(this%rad_grid%rr(i), dT, dv, v, nlm)
    
    deallocate( dv, dT )
    
    select case (this%scaling)
      case ('basics')
        do concurrent ( ijm = 1:this%jms, i1 = 2:4 )
          nlm(i1,ijm) = nlm(i1,ijm) / this%Pr
        end do
    end select
    
    call this%coriolis_rr_jml_sub(v, nlm(2:4,:))
    call this%buoy_rr_jml_sub(i, T, nlm(2:4,:))
    
    deallocate( v, T )
    
    do concurrent ( ijm = 1:this%jms )
      this%ntemp(ijm,i) = nlm(1,ijm)
      this%nsph1(ijm,i) = nlm(2,ijm)
      this%ntorr(ijm,i) = nlm(3,ijm)
      this%nsph2(ijm,i) = nlm(4,ijm)
    end do
    
    deallocate( nlm )
    
  end subroutine fullnl_sub

end submodule NonLinearTerms