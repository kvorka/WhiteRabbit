submodule (IceTidesMod) Parameters_iceTides
  implicit none; contains
  
  module subroutine visc_iceTides_jm_sub(this)
    class(T_iceTides), intent(inout) :: this
    integer                          :: ir, i1, i2, i3
    real(kind=dbl)                   :: visc
    real(kind=dbl),    allocatable   :: grid(:,:,:)
    complex(kind=dbl), allocatable   :: visc_jm(:), cc_mj(:)
    
    if ( this%mparams%initvisc ) then
      if ( this%mparams%visc_radial) then
        !$omp parallel do private(visc)
        do ir = 1, this%nd
          visc = s4pi * this%viscU / c2r_fn( this%mparams%visc(1,ir) )
          this%mparams%visc(1,ir) = r2c_fn( s4pi * this%viscU / andrade_visc_fn(this%mu, this%omega, visc) )
          !this%mparams%visc(1,ir) = r2c_fn( s4pi * this%viscU / visc )
        end do
        !$omp end parallel do
        
      else
        !!************!!
        !! TO DO      !!
        !!************!!
      end if
    end if
    
  end subroutine visc_iceTides_jm_sub
  
end submodule Parameters_iceTides