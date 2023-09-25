submodule(OutputOceanMod) Flux
  implicit none
  
  contains
  
  subroutine harm_analysis_flux_sub()
    integer                        :: ijm, error
    real(kind=dbl),    allocatable :: data_flux(:,:)
    complex(kind=dbl)              :: fluxjm
    complex(kind=dbl), allocatable :: spectra(:)
    
    allocate( spectra(jms) ) ; spectra = czero
    
    open(unit=1, file='flux-averaged.spec', status='old', action='read')
      do
        read(1,*,iostat=error) ijm, fluxjm
        
        if (error /= 0) then
          exit
        else
          spectra(ijm) = fluxjm
        end if
      end do
    close(1)
    
    allocate( data_flux(2*nth,nth) ); data_flux = 0._dbl
    
    call harmsy_sub(jmax_ocean, spectra, data_flux)
    
    deallocate( spectra )
    
    open(unit=8, file='ocean_flux_grid', status='new', action='write')
      write(8,*) 'max: ' , maxval(data_flux)
      write(8,*) 'min: ' , minval(data_flux)
    close(8)
    
    call out_data_2d_sub('inp-f.dat', data_flux)
    call out_data_1d_sub('inp-f-zon.dat', data_flux)
    
    deallocate( data_flux )
    
  end subroutine harm_analysis_flux_sub

  subroutine save_spectra_flux_sub()
    complex(kind=dbl), allocatable :: flux(:)
    
    allocate( flux(jms) ) ; flux = czero
    
    call avrg_spectra_2d_sub(path_ocean_flux, flux)
    call out_spectra_2d_sub('flux-averaged.spec', flux)
    
    deallocate( flux )
    
  end subroutine save_spectra_flux_sub
  
end submodule Flux