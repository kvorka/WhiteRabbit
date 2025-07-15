submodule (oceanconv) iter
  implicit none; contains
  
  module procedure iter_oceanConv_sub
    integer :: k, ijm
    
    do k = 1, this%n_iter
      this%t = this%t + this%dt
        call this%time_scheme_sub()
    end do
    
    call this%vypis_ocean_sub()
    
  end procedure iter_oceanConv_sub
  
end submodule iter