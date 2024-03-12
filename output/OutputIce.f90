program OutputIce
  use IceConstants
  use OutputMod
  implicit none
  
  call save_spectra_ice_sub(path_ice_shape_dn, 'shape_dn.spec')
  call save_spectra_ice_sub(path_ice_shape_up, 'shape_up.spec')
  call save_spectra_ice_sub(path_ice_topo_dn , 'topo_dn.spec' )
  call save_spectra_ice_sub(path_ice_topo_up , 'topo_up.spec' )
  
  contains
  
  subroutine save_spectra_ice_sub(path_in, path_out)
    character(len=*),  intent(in)  :: path_in, path_out
    complex(kind=dbl), allocatable :: spectra(:)
    
    allocate( spectra(jmax_ice*(jmax_ice+1)/2+jmax_ice+1) )
    spectra = czero

    call avrg_spectra_2d_sub(path_in, spectra)
    call out_spectra_2d_sub(path_out, spectra)
    
    deallocate( spectra )
    
  end subroutine save_spectra_ice_sub
  
end program OutputIce
