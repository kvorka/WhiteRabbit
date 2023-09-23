submodule(OutputOceanMod) Temp
  implicit none
  
  contains
  
  subroutine harm_analysis_temp_sub()
    integer,             parameter :: jmax = 100
    integer                        :: ir, ij, im
    real(kind=dbl),    allocatable :: map(:,:), temp(:,:)
    complex(kind=dbl), allocatable :: temp120(:,:)
    
    allocate( temp120(jms, n_out), map(2*nth,nth), temp(0:nth,n_out) )
      
    call load_data_sub(jms, 'temp-averaged.spec', temp120)
    
    do ij = 0, jmax_ocean
      do im = 0, ij
        if ( (mod(ij,2) /= 0) .or. (im /= 0) ) temp120(jm(ij,im),:) = czero
      end do
    end do
    
    do ir = 1, n_out
      call harmsy_sub(jmax, temp120(:,ir), map) ; call get_zondata_sub(map, temp(:,ir))
    end do
    
    deallocate( temp120, map )
    
    open(unit=8, file='ocean_temp_grid', status='new', action='write')
      write(8,*) 'max: ' , maxval( temp )
      write(8,*) 'min: ' , minval( temp )
    close(8)
    
    call save_data_sub('inp-t.dat', temp)
    
    deallocate( temp )
    
  end subroutine harm_analysis_temp_sub
  
  subroutine save_spectra_temp_sub()
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: temp(:,:)
    
    allocate( r(nd_ocean+1), temp(jms,nd_ocean+1) ) ; temp = czero
    
    call avrg_spectra_3d_sub(path_ocean_temp, r, temp)
    call out_spectra_3d_sub('temp-averaged.spec', r, temp)
    
    deallocate( r, temp )
    
  end subroutine save_spectra_temp_sub
  
end submodule Temp
