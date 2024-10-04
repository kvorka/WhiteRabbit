submodule (icetides) solvers
  implicit none; contains
  
  module subroutine compute_iceTides_sub(this, visc_prof)
    class(T_iceTides), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: visc_prof(:,:)
    integer                          :: ijm, ir, n
    real(kind=dbl)                   :: P, Pglobal
    
    call this%sol%nulify_sub() ; Pglobal = zero
    
    if ( this%mparams%initvisc ) then
      if ( this%mparams%visc_radial) then
        do concurrent ( ir = 1:this%nd )
          this%mparams%visc(1,ir) = visc_prof(1,ir)
        end do
      else
        do concurrent ( ijm = 1:this%jms, ir = 1:this%nd )
          this%mparams%visc(ijm,ir) = visc_prof(ijm,ir)
        end do
      end if
    end if
    
    call this%set_visc_sub()
    
    this%t            = zero
    this%dt           = this%period / this%n_iter
    this%tdheat%htide = czero
    
    do
      do n = 1, this%n_iter
        this%t = this%t + this%dt
        
        call this%EE_mech_sub()
        call this%tidal_heating_4_sub()
      end do
          
      P = this%rad_grid%intV_fn( c2r_fn( this%tdheat%htide(:,1) ) )
        if ( abs(P-Pglobal) / P < 1.0d-4 ) then
          exit
        else
          Pglobal = P
          this%tdheat%htide = czero
        end if
    end do
    
  end subroutine compute_iceTides_sub
  
end submodule solvers