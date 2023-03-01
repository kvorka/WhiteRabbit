program OutputOcean
  use OutputOceanMod
  implicit none

  call nuss_curve_sub(path_nuss)

  call save_spectra_flux_sub(path_ocean_flux)
  call save_spectra_temp_sub(path_ocean_temp)
  call save_spectra_velc_sub(path_ocean_velc)
  
  call harm_analysis_flux_sub('')
  call harm_analysis_temp_sub('')
  call harm_analysis_rad_velc_sub('')
  call harm_analysis_zon_velc_sub('')
  
end program OutputOcean
