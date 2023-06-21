submodule(IceMod) TidalHeatingIce
  implicit none

  contains

  pure complex(kind=dbl) function htide_ice_4_fn(this, ir, ijm)
    class(T_ice),      intent(in) :: this
    integer,           intent(in) :: ir, ijm
    complex(kind=dbl)             :: HI
    
    if ( ijm <= jms4 ) then
      HI = this%rad_grid%cc(ir,-1) * this%htide(ir-1,ijm) + this%rad_grid%cc(ir,+1) * this%htide(ir,ijm)
    else
      HI = czero
    end if
    
    htide_ice_4_fn = this%Ds/this%Ra * HI / this%cp_fn(ir)
    
  end function htide_ice_4_fn

  subroutine tidal_heating_ice_4_sub(this)
    class(T_ice),      intent(inout) :: this
    integer                          :: ir
    complex(kind=dbl), allocatable   :: Dstrss(:), H(:)
    
    allocate( Dstrss(jml2(this%jmax,this%jmax,this%jmax+2)), H(jms4) )
    
    do ir = 1, this%nd
      H      = czero
      Dstrss = this%sol%deviatoric_stress_jml2_fn(ir)
      
      H(1) =     Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0,-2))  / sqrt(4*pi) + &
           & 2 * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2,-2))) / sqrt(4*pi) + &
           &     Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  / sqrt(4*pi) + &
           & 2 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) / sqrt(4*pi) + &
           &     Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  / sqrt(4*pi) + &
           & 2 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2))) / sqrt(4*pi)   ; H(1)%im = 0._dbl
      
      H(4) = - 2 / sqrt(14*pi)        * Dstrss(jml2(2,0,-2)) *       Dstrss(jml2(2,0, 0))  + &
           &   2 / sqrt(14*pi)        * Dstrss(jml2(2,2,-2)) * conjg(Dstrss(jml2(2,2, 0))) + &
           &   2 / sqrt(14*pi)        * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,-2))) - &
           &   3 * sqrt(5/pi/49) / 14 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0, 0))  + &
           &   3 * sqrt(5/pi/49) /  7 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2, 0))) - &
           &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,0, 0)) *       Dstrss(jml2(2,0,+2))  + &
           &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2, 0)) * conjg(Dstrss(jml2(2,2,+2))) + &
           &  12 * sqrt(1/pi)    / 49 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2, 0))) + &
           &   5 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,0,+2)) *       Dstrss(jml2(2,0,+2))  - &
           &  10 * sqrt(5/pi)    / 49 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2)))   ; H(4)%im = 0._dbl
      
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
            &  6 / sqrt(   pi)   /196 * Dstrss(jml2(2,2,+2)) * conjg(Dstrss(jml2(2,2,+2)))   ; H(11)%im = 0._dbl
      
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
      
      this%htide(ir,:) = this%htide(ir,:) + H(:) / this%visc_fn(ir) / 2 / this%n_iter
    end do
        
    deallocate( Dstrss, H )
    
  end subroutine tidal_heating_ice_4_sub

end submodule TidalHeatingIce