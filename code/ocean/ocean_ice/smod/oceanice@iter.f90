submodule (oceanice) iter
  implicit none; contains
  
  module procedure iter_oceanIce_sub
    integer        :: k, ijm
    real(kind=dbl) :: avrg_flux
    
    call zero_carray_sub( this%jms, this%bnd%flux_up )
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
    end do
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
        
        do concurrent ( ijm = 1:this%jms )
          this%bnd%flux_up(ijm) = this%bnd%flux_up(ijm) + this%qr_r_fn(this%nd,ijm)
        end do
    end do
    
    avrg_flux = this%bnd%flux_up(1)%re / s4pi
      do concurrent ( ijm = 1:this%jms )
        this%bnd%flux_up(ijm) = this%bnd%flux_up(ijm) / avrg_flux
      end do
    
    call this%vypis_ocean_sub()
    
  end procedure iter_oceanIce_sub
  
end submodule iter