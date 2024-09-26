submodule (OceanMod) Iter_ocean
  implicit none; contains
  
  module subroutine iter_ocean_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: k, ijm
    real(kind=dbl)                :: avrg_flux
    
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
    
  end subroutine iter_ocean_sub
  
  module subroutine speed_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: k
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
    end do
    
  end subroutine speed_sub
  
end submodule Iter_ocean