program OutputOcean
  use OutputOceanMod
  implicit none

  call nuss_curve_sub('data/Nuss.dat')

  call save_spectra_flux_sub('data/data_ocean_flux/Flux-')
  call save_spectra_temp_sub('data/data_ocean_temp/Temp-')
  call save_spectra_velc_sub('data/data_ocean_veloc/Velc-')
  
  call harm_analysis_flux_sub('')
  call harm_analysis_temp_sub('')
  call harm_analysis_rad_velc_sub('')
  call harm_analysis_zon_velc_sub('')
  
end program OutputOcean
