submodule(OutputOceanMod) Flux
  implicit none

  contains

  subroutine harm_analysis_flux_sub(path)
    character(len=*),   intent(in) :: path
    integer                        :: jm, error
    real(kind=dbl),    allocatable :: data_flux(:,:)
    complex(kind=dbl)              :: fluxjm
    complex(kind=dbl), allocatable :: spectra(:)
    
    allocate( spectra(jmax_ocean*(jmax_ocean+1)/2+jmax_ocean+1) )
      spectra = cmplx(0._dbl, 0._dbl, kind=dbl)
        open(unit=1, file='flux-averaged.spec', status='old', action='read')
          do
            read(1,*,iostat=error) jm, fluxjm ; if (error /= 0) exit
            spectra(jm) = fluxjm
          end do
        close(1)

    allocate( data_flux(2*nth,nth) ); data_flux = 0._dbl
      call toGrid_sub(jmax_ocean, spectra, data_flux)
    
    deallocate( spectra )
      open(unit=8, file='ocean_flux_range_max', status='new', action='write')
        write(8,'(1f9.2)') maxval(data_flux)
      close(8)
      open(unit=8, file='ocean_flux_range_min', status='new', action='write')
        write(8,'(1f9.2)') minval(data_flux)
      close(8)

      data_flux = data_flux / max( maxval(data_flux), abs(minval(data_flux)) )

      open(unit=8, file='ocean_flux_range_max_grid', status='new', action='write')
        if (maxval(data_flux) > treshold) then
          write(8,'(1f4.2)') treshold
        else
          write(8,'(1f4.2)') maxval(data_flux)
        end if
      close(8)
      open(unit=8, file='ocean_flux_range_min_grid', status='new', action='write')
        if (minval(data_flux) < -treshold) then
          write(8,'(1f5.2)') -treshold
        else
          write(8,'(1f5.2)') minval(data_flux)
        end if
      close(8)

      call out_data_sub('inp-f.dat', data_flux)
    
    deallocate( data_flux )
    
  end subroutine harm_analysis_flux_sub

  subroutine save_spectra_flux_sub(path)
    character(len=*),  intent(in)  :: path
    integer                        :: j, m, n, error
    complex(kind=dbl)              :: flux_t
    complex(kind=dbl), allocatable :: flux(:)

    allocate( flux(jmax_ocean*(jmax_ocean+1)/2+jmax_ocean+1) )

      do n = avrg_start, avrg_end
        open( unit=7, file=path//trim(adjustl(int2str_fn(n)))//'.dat', status='old', action='read')
          do
            read(7,*,iostat=error) j, m, flux_t; if (error /= 0) exit
            flux(j*(j+1)/2+m+1) = flux(j*(j+1)/2+m+1) + flux_t/(avrg_end-avrg_start)
          end do
        close(7)
      end do

      flux(1) = cmplx(0._dbl, 0._dbl, kind=dbl) ; call out_spectra1_sub('flux-averaged.spec', flux)

    deallocate( flux )

  end subroutine save_spectra_flux_sub

end submodule Flux
