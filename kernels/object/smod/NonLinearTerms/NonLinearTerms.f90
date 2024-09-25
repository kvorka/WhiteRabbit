submodule (PhysicalObject) NonLinearTerms
  implicit none ; contains
  
  module pure subroutine coriolis_sub(this, i)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: i
    integer                                :: ijm
    complex(kind=dbl),       allocatable   :: v(:), nlm(:,:)
    
    allocate( v(this%jmv) ) ; call this%v_rr_ijml_sub(i, v)
    
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
  
  module subroutine coriolis_vgradv_sub(this, i)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: i
    integer                                :: ijm, i1
    complex(kind=dbl),       allocatable   :: v(:), dv(:), nlm(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv) ) ; call this%dv_dr_rr_jml_sub(i, v, dv)
    
    allocate( nlm(3,this%jms) ) ; call this%lat_grid%vcvgv_sub(this%rad_grid%rr(i), dv, v, nlm)
    
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
  
  module subroutine fullnl_sub(this, i)
    class(T_physicalObject), intent(inout) :: this
    integer,                 intent(in)    :: i
    integer                                :: ijm, i1
    real(kind=dbl)                         :: fac
    complex(kind=dbl),       allocatable   :: v(:), dv(:), T(:), gradT(:), nlm(:,:)
    
    allocate( v(this%jmv) , dv(this%jmv) ) ; call this%dv_dr_rr_jml_sub(i, v, dv)
    allocate( T(this%jms) , gradT(this%jmv) ) ; call this%gradT_rr_ijml_sub(i, T, gradT, -1)
    
    allocate( nlm(4,this%jms) ) ; call this%lat_grid%vcvv_vcvgv_sub(this%rad_grid%rr(i), gradT, dv, v, nlm)
    
    deallocate( dv, gradT )
    
    select case (this%scaling)
      case ('basics')
        fac = 1 / this%Pr
        
        do concurrent ( ijm = 1:this%jms, i1 = 2:4 )
          nlm(i1,ijm) = nlm(i1,ijm) * fac
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
  
  module subroutine tidal_heating_4_sub(this)
    class(T_physicalObject), intent(inout) :: this
    integer                                :: ir
    complex(kind=dbl),       allocatable   :: Dstrss(:), H(:)
    
    allocate( Dstrss(jml2(this%jmax,this%jmax,this%jmax+2)), H(this%tdheat%jms) )
    
    !$omp parallel do private (Dstrss, H)
    do ir = 1, this%nd
      H      = czero
      Dstrss = this%sol%deviatoric_stress_jml2_fn(ir)
      
      H(1) =     Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0,-2))  / s4pi + &
           & 2 * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2,-2))) / s4pi + &
           &     Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  / s4pi + &
           & 2 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) / s4pi + &
           &     Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  / s4pi + &
           & 2 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2))) / s4pi   ; H(1)%im = zero
      
      H(4) = - 2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0, 0))  + &
           &   2 / sqrt(14*pi)        * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2, 0))) + &
           &   2 / sqrt(14*pi)        * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,-2))) - &
           &   3 * sqrt(5/pi/49) / 14 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  + &
           &   3 * sqrt(5/pi/49) /  7 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) - &
           &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0,+2))  + &
           &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,+2))) + &
           &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2, 0))) + &
           &   5 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  - &
           &  10 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2)))   ; H(4)%im = zero
      
      H(6) =  2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) * Dstrss(jml2(2,2, 0)) + &
           &  2 / sqrt(14*pi)        * Dstrss(jml2(2,2,-2)) * Dstrss(jml2(2,0, 0)) + &
           &  3 * sqrt(5/pi/49) /  7 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2, 0)) + &
           & 12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2,+2)) + &
           & 12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,0,+2)) - &
           & 10 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,0,+2)) * Dstrss(jml2(2,2,+2))
      
      H(11) =  2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0,+2))  + &
            &  1 / sqrt(14*pi)   /  3 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,-2))) + &
            &  1 / sqrt(14*pi)   /  3 * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2,+2))) + &
            &  6 / sqrt(   pi)   / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  + &
            &  2 / sqrt(   pi)   / 49 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) - &
            & 50 / sqrt( 5*pi)   / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0,+2))  - &
            & 25 / sqrt( 5*pi)   /147 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,+2))) - &
            & 25 / sqrt( 5*pi)   /147 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2, 0))) + &
            &  9 / sqrt(   pi)   / 98 * Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  + &
            &  6 / sqrt(   pi)   /196 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2)))   ; H(11)%im = zero
      
      H(13) =      sqrt( 5/pi/42)      * Dstrss(jml2(2,0,-2)) * Dstrss(jml2(2,2,+2)) + &
            &      sqrt( 5/pi/42)      * Dstrss(jml2(2,2,-2)) * Dstrss(jml2(2,0,+2)) + &
            &  2 * sqrt(15/pi)    / 49 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2, 0)) - &
            & 25 * sqrt( 3/pi)    /147 * Dstrss(jml2(2,0, 0)) * Dstrss(jml2(2,2,+2)) - &
            & 25 * sqrt( 3/pi)    /147 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,0,+2)) + &
            &  3 * sqrt(15/pi)    / 98 * Dstrss(jml2(2,0,+2)) * Dstrss(jml2(2,2,+2))
      
      H(15) =      sqrt( 5/pi)    /  3 * Dstrss(jml2(2,2,-2)) * Dstrss(jml2(2,2,+2)) + &
            &      sqrt(10/pi/7)  /  7 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,2, 0)) - &
            & 50 / sqrt(14*pi)    / 21 * Dstrss(jml2(2,2, 0)) * Dstrss(jml2(2,2,+2)) + &
            &  3 * sqrt( 5/pi/14) / 14 * Dstrss(jml2(2,2,+2)) * Dstrss(jml2(2,2,+2))
      
      this%tdheat%htide(ir,:) = this%tdheat%htide(ir,:) + H(:) / this%visc_r_fn(ir) / 2 / this%n_iter
    end do
    !$omp end parallel do
    
    deallocate( Dstrss, H )
    
  end subroutine tidal_heating_4_sub

end submodule NonLinearTerms