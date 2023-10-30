submodule(OutputOceanMod) Temp
  implicit none
  
  contains
  
  subroutine harm_analysis_temp_sub()
    integer,           parameter   :: jmax = 100
    integer                        :: ir, ij
    real(kind=dbl),    allocatable :: temp_j0(:), temp_grid(:,:)
    complex(kind=dbl), allocatable :: temp120(:,:)
    
    allocate( temp120(jms,n_out), temp_j0(0:jmax), temp_grid(nth,n_out) )
    
      call load_data_sub(jms, 'temp-averaged.spec', temp120) ; temp_j0 = zero
      
      do ir = 1, n_out
        do ij = 0, jmax, 2
          temp_j0(ij) = c2r_fn( temp120( jm(ij,0) , ir ) )
        end do
        
        call harmsy_Pj0_sub(jmax, temp_j0, temp_grid(:,ir))
      end do
    
    deallocate( temp120, temp_j0 )
    
    open(unit=8, file='ocean_temp_extrms', status='new', action='write')
      write(8,*) 'max: ' , maxval( temp_grid )
      write(8,*) 'min: ' , minval( temp_grid )
    close(8)
    
    call save_data_sub('inp-t.dat', temp_grid)
    
    deallocate( temp_grid )
    
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
