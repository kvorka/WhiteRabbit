submodule(OutputOceanMod) Temp
  implicit none

  contains

  subroutine harm_analysis_temp_sub()
    integer,             parameter :: jmax = 100
    integer                        :: i, k, j, m
    real(kind=dbl),    allocatable :: map(:,:), t(:,:)
    complex(kind=dbl), allocatable :: temp120(:,:)

    allocate( temp120(jms, n_out), map(2*nth,nth), t(0:nth,n_out) )

      call load_data_sub(jms, 'temp-averaged.spec', temp120)

      do j = 0, jmax_ocean
        do m = 0, j
          if ( (mod(j,2) /= 0) .or. (m /= 0) ) temp120(j*(j+1)/2+m+1,:) = cmplx(0._dbl, 0._dbl, kind=dbl)
        end do
      end do

      do i = 1, n_out
        call toGrid_sub(jmax, temp120(:,i), map) ; call get_zonal_sub(map, t(:,i))
      end do

      open(unit=8, file='ocean_t_max', status='new', action='write') ; write(8,'(1g9.2)') maxval(t) ; close(8)
      open(unit=8, file='ocean_t_min', status='new', action='write') ; write(8,'(1g9.2)') minval(t) ; close(8)

      t = t / max( maxval(t), abs(minval(t)) )
        open(unit=8, file='ocean_t_max_g', status='new', action='write'); write(8,'(1f4.2)') min(maxval(t), +treshold); close(8)
        open(unit=8, file='ocean_t_min_g', status='new', action='write'); write(8,'(1f5.2)') max(minval(t), -treshold); close(8)

      call save_data_sub('inp-t.dat', t)
    
    deallocate(t, temp120, map)

  end subroutine harm_analysis_temp_sub
  
  subroutine save_spectra_temp_sub()
    integer                        :: n, i
    real(kind=dbl),    allocatable :: r_init(:)
    complex(kind=dbl), allocatable :: temp(:,:), temp_i(:)

    allocate( r_init(nd_ocean+1), temp(jms,nd_ocean+1), temp_i(jms) )
      
      temp = cmplx(0._dbl, 0._dbl, kind=dbl)

      do n = avrg_start, avrg_end
        open(unit=7, file=path_ocean_temp//trim(adjustl(int2str_fn(n)))//'.dat', status='old', action='read')
          do i = 1, nd_ocean+1
            read(7,*) r_init(i) , temp_i(:)
            temp(:,i) = temp(:,i) + temp_i(:) / (avrg_end-avrg_start)
          end do
        close(7)
      end do

      call out_spectra_sub('temp-averaged.spec', r_init, temp)

    deallocate( r_init, temp, temp_i )

  end subroutine save_spectra_temp_sub

end submodule Temp
