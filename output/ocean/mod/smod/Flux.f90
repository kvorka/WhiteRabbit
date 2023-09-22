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
    
    call toGrid_sub(jmax_ocean, spectra, data_flux)
    
    deallocate( spectra )
    
    open(unit=8, file='ocean_flux_grid', status='new', action='write')
      write(8,*) 'max: ' , maxval(data_flux)
      write(8,*) 'min: ' , minval(data_flux)
    close(8)
    
    call out_data_sub('inp-f.dat', data_flux)
    
    deallocate( data_flux )
    
  end subroutine harm_analysis_flux_sub

  subroutine save_spectra_flux_sub()
    integer                        :: ij, im, in, error
    complex(kind=dbl)              :: flux_t
    complex(kind=dbl), allocatable :: flux(:)
    
    allocate( flux(jms) )
    
    do in = avrg_start, avrg_end
      open( unit=7, file=path_ocean_flux//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read')
      
      do
        read(7,*,iostat=error) ij, im, flux_t
        
        if (error /= 0) then 
          exit
        else
          flux(jm(ij,im)) = flux(jm(ij,im)) + flux_t / (avrg_end-avrg_start)
        end if
      end do
      
      close(7)
    end do
    
    call out_spectra1_sub('flux-averaged.spec', flux)
      
    deallocate( flux )
    
  end subroutine save_spectra_flux_sub
  
end submodule Flux