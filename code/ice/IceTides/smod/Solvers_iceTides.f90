submodule (IceTidesMod) Solvers_iceTides
  implicit none; contains
  
  module subroutine compute_iceTides_sub(this, visc_prof)
    class(T_iceTides), intent(inout) :: this
    complex(kind=dbl), intent(in)    :: visc_prof(this%nd,*)
    integer                          :: ijm, ir, n
    real(kind=dbl)                   :: P, Pglobal
    
    call this%sol%nulify_sub() ; Pglobal = zero
    
    if ( this%mparams%initvisc ) then
      do concurrent ( ijm = 1:this%jms, ir = 1:this%nd )
        this%mparams%visc(ir,ijm) = visc_prof(ir,ijm)
      end do
    else
      do concurrent ( ir = 1:this%nd )
        this%mparams%visc_radial(ir) = visc_prof(ir,1)
      end do
    end if
    
    this%t = zero
    this%dt = this%period / this%n_iter ; this%htide = czero
    
    do
      do n = 1, this%n_iter
        this%t = this%t + this%dt
        
        call this%EE_mech_sub()
        call this%tidal_heating_sub()
      end do
          
      P = this%rad_grid%intV_fn( real(this%htide(:,1), kind=dbl) )
        if ( abs(P-Pglobal) / P < 1.0d-3 ) then
          exit
        else
          Pglobal    = P
          this%htide = czero
        end if
    end do
    
  end subroutine compute_iceTides_sub
  
end submodule Solvers_iceTides