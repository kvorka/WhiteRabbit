submodule (oceantides) iter
  implicit none; contains
  
  module subroutine iter_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this
    integer                            :: k
    
    this%heating = zero; this%number_of_periods = this%number_of_periods + 1

    do k = 1, this%n_iter
      this%k_of_period = k ; call this%time_scheme_sub()
      this%heating = this%heating + this%viscdissip_power_fn()
    end do

    this%heating = this%heating / this%n_iter

    call vypis_oceanTides_sub(this)

  end subroutine iter_oceanTides_sub
  
end submodule iter