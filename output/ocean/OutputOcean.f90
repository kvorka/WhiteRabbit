program OutputOcean
  use OutputOceanMod
  implicit none
  
  call nuss_curve_sub()
  
  call save_spectra_flux_sub()
  call save_spectra_temp_sub()
  call save_spectra_velc_sub()
  
end program OutputOcean
