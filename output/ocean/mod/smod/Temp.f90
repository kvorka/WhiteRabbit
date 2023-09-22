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
      call toGrid_sub(jmax, temp120(:,ir), map) ; call get_zonal_sub(map, temp(:,ir))
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
    integer                        :: in, ir
    real(kind=dbl),    allocatable :: r_init(:)
    complex(kind=dbl), allocatable :: temp(:,:), temp_i(:)
    
    allocate( r_init(nd_ocean+1), temp(jms,nd_ocean+1), temp_i(jms) ) ; temp = czero
    
    do in = avrg_start, avrg_end
      open(unit=7, file=path_ocean_temp//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read')
      
      do ir = 1, nd_ocean+1
        read(7,*) r_init(ir) , temp_i(:)
        temp(:,ir) = temp(:,ir) + temp_i(:) / (avrg_end-avrg_start)
      end do
      
      close(7)
    end do
    
    call out_spectra_sub('temp-averaged.spec', r_init, temp)
    
    deallocate( r_init, temp, temp_i )
    
  end subroutine save_spectra_temp_sub
  
end submodule Temp
