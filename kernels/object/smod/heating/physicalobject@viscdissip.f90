submodule (physicalobject) viscdissip
  implicit none ; contains
  
  module procedure viscdissip_power_fn
    integer                     :: ir
    real(kind=dbl), allocatable :: power_ir(:)
    
    allocate( power_ir(this%nd) )
    
      !$omp parallel do
      do ir = 1, this%nd
        power_ir(ir)  = tensnorm2_fn( this%jmax, this%sol%deviatoric_stress_jml2_fn(ir) ) / this%visc_r_fn(ir) / 2
      end do
      !$omp end parallel do
      
      viscdissip_power_fn = this%rad_grid%intV_fn( power_ir )
    
    deallocate( power_ir )
    
  end procedure viscdissip_power_fn
  
  module procedure tidal_heating_4_sub
    integer                        :: ir
    complex(kind=dbl), allocatable :: Dstrss(:), H(:)
    
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
    
  end procedure tidal_heating_4_sub
  
end submodule viscdissip