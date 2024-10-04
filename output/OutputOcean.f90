program OutputOcean
  use ocean_constants
  use OutputMod
  implicit none
  
  integer, parameter :: n_out = 120
  integer, parameter :: jms   =    (jmax_ocean  )*(jmax_ocean+1)/2 + (jmax_ocean  )  + 1
  integer, parameter :: jms1  =    (jmax_ocean+1)*(jmax_ocean+2)/2 + (jmax_ocean+1)  + 1
  integer, parameter :: jmv   = 3*((jmax_ocean  )*(jmax_ocean+1)/2 + (jmax_ocean  )) + 1
  
  call nuss_curve_sub()
  
  call save_spectra_flux_sub()
  call save_spectra_temp_sub()
  call save_spectra_velc_sub()
  
  contains
  
  subroutine nuss_curve_sub()
    integer         :: n, error
    real(kind=dbl)  :: t, dt, Nuss, Re, sumNuss, sumRe
    
    n = 0
      sumNuss = zero
      sumRe   = zero
    
    open(unit=1, file=path_nuss, status='old', action='read')
    
    do
      read(1,*,iostat=error) t, dt, Nuss, Re
      
      if ( error /= 0 ) then
        exit
      else if ( t > tNuss ) then
        n = n + 1
          sumNuss = sumNuss + Nuss
          sumRe   = sumRe + Re
      end if
    end do
    
    close(1)
    
    open(unit=8, file='nuss', status='new', action='write')
      write(8,*) sumNuss / n , sumRe / n
    close(8)
    
  end subroutine nuss_curve_sub
  
  subroutine save_spectra_flux_sub()
    complex(kind=dbl), allocatable :: flux(:)
    
    allocate( flux(jms) ) ; flux = czero
    
    call avrg_spectra_2d_sub(path_ocean_flux, flux)
    call out_spectra_2d_sub('flux-averaged.spec', flux)
    
    deallocate( flux )
    
  end subroutine save_spectra_flux_sub
  
  subroutine save_spectra_temp_sub()
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: temp(:,:)
    
    allocate( r(nd_ocean+1), temp(jms,nd_ocean+1) ) ; temp = czero
    
    call avrg_spectra_3d_sub(path_ocean_temp, r, temp)
    call out_spectra_3d_sub('temp-averaged.spec', r, temp)
    
    deallocate( r, temp )
    
  end subroutine save_spectra_temp_sub
  
  subroutine save_spectra_velc_sub()
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: velc(:,:)
    
    allocate( r(nd_ocean+1), velc(jmv,nd_ocean+1) ) ; velc = czero
    
    call avrg_spectra_3d_sub(path_ocean_velc, r, velc)
    call out_spectra_3d_sub('velc-averaged.spec', r, velc)
    
    deallocate( r, velc )
    
  end subroutine save_spectra_velc_sub
  
end program OutputOcean
