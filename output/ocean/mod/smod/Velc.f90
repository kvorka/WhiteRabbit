submodule(OutputOceanMod) Velc
  implicit none
  
  contains
  
  subroutine harm_analysis_rad_velc_sub()
    integer,              parameter :: jmax = 100
    integer                         :: ir, iph, ith, ij, im, il
    real(kind=dbl),    allocatable  :: comp_x(:,:), comp_y(:,:), comp_z(:,:), comp_r(:,:), comp_t(:,:), comp_p(:,:), rvelc(:,:)
    complex(kind=dbl), allocatable  :: x(:), y(:), z(:), velc120(:,:)
    
    allocate( velc120(jmv,n_out), x(jms1), y(jms1), z(jms1), comp_x(2*nth,nth), comp_y(2*nth,nth), comp_z(2*nth,nth), &
            & comp_r(2*nth,nth), comp_t(2*nth,nth), comp_p(2*nth,nth), rvelc(0:nth,n_out)                             )
    
    call load_data_sub(jmv, 'velc-averaged.spec', velc120)
    
    do ij = 0, jmax_ocean
      do im = 0, ij
        do il = abs(ij-1)-ij, +1
          if( (mod(ij,2) /= 0) .or. (im /= 0) .or. (il == 0) ) velc120(jml(ij,im,il),:) = czero
        end do
      end do
    end do
    
    do ir = 1, n_out
      call vec2scals_sub(jmax_ocean, velc120(:,ir), x, y, z); y = -y
      
      call toGrid_sub(jmax, x, comp_x)
      call toGrid_sub(jmax, y, comp_y)
      call toGrid_sub(jmax, z, comp_z)
      
      do ith = 1, nth
        do iph = 1, 2*nth
          call vecxyz2rtp_sub( 180.5_dbl-ith, iph-1._dbl, comp_x(iph,ith), comp_y(iph,ith), comp_z(iph,ith), &
                             &                            comp_r(iph,ith), comp_t(iph,ith), comp_p(iph,ith)  )
        end do
      end do
      
      call get_zonal_sub(comp_r, rvelc(:,ir))
    end do
    
    deallocate( velc120, x, y, z, comp_x, comp_y, comp_z, comp_r, comp_t, comp_p )
    
    rvelc = Ek_ocean * rvelc / Pr_ocean
      open(unit=8, file='ocean_rvelc_grid', status='new', action='write')
        write(8,*) 'max: ' , maxval( rvelc )
        write(8,*) 'min: ' , minval( rvelc )
      close(8)
    
    call save_data_sub('inp-r.dat', rvelc)
    
    deallocate( rvelc )
    
  end subroutine harm_analysis_rad_velc_sub
  
  subroutine harm_analysis_zon_velc_sub()
    integer,              parameter :: jmax = 100
    integer                         :: ir, iph, ith, ij, im, il
    real(kind=dbl),    allocatable  :: comp_x(:,:), comp_y(:,:), comp_z(:,:), comp_r(:,:), comp_t(:,:), comp_p(:,:), zvelc(:,:)
    complex(kind=dbl), allocatable  :: x(:), y(:), z(:), velc120(:,:)
    
    allocate( velc120(jmv,n_out), x(jms1), y(jms1), z(jms1), comp_x(2*nth,nth), comp_y(2*nth,nth), comp_z(2*nth,nth), &
            & comp_r(2*nth,nth), comp_t(2*nth,nth), comp_p(2*nth,nth), zvelc(0:nth,n_out)                             )
    
    call load_data_sub(jmv, 'velc-averaged.spec', velc120)
    
    do ij = 0, jmax_ocean
      do im = 0, ij
        do il = abs(ij-1)-ij, +1
          if( (mod(ij,2) == 0) .or. (im /= 0) .or. (il /= 0) ) velc120(jml(ij,im,il),:) = czero
        end do
      end do
    end do
    
    do ir = 1, n_out
      call vec2scals_sub(jmax_ocean, velc120(:,ir), x, y, z); y = -y
      
      call toGrid_sub(jmax, x, comp_x)
      call toGrid_sub(jmax, y, comp_y)
      call toGrid_sub(jmax, z, comp_z)
      
      do ith = 1, nth
        do iph = 1, 2*nth
          call vecxyz2rtp_sub( 180.5_dbl-ith, iph-1._dbl, comp_x(iph,ith), comp_y(iph,ith), comp_z(iph,ith), &
                             &                            comp_r(iph,ith), comp_t(iph,ith), comp_p(iph,ith)  )
        end do
      end do
      
      call get_zonal_sub(comp_p, zvelc(:,ir))
    end do
    
    deallocate( velc120, x, y, z, comp_x, comp_y, comp_z, comp_r, comp_t, comp_p )
    
    zvelc = Ek_ocean * zvelc / Pr_ocean
      open(unit=8, file='ocean_zvelc_grid', status='new', action='write')
        write(8,*) 'max: ' , maxval( zvelc )
        write(8,*) 'min: ' , minval( zvelc )
      close(8)
    
    call save_data_sub('inp-p.dat', zvelc)
    
    deallocate( zvelc )
    
  end subroutine harm_analysis_zon_velc_sub
  
  subroutine save_spectra_velc_sub()
    integer                        :: in, ir
    real(kind=dbl),    allocatable :: r_init(:), r_out(:)
    complex(kind=dbl), allocatable :: velc(:,:), velc_i(:), velc_out(:,:)
    
    allocate( r_init(nd_ocean+1), velc(jmv,nd_ocean+1), velc_i(jmv) ) ; velc = czero
    
    do in = avrg_start, avrg_end
      open( unit=7, file=path_ocean_velc//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read' )
      
      do ir = 1, nd_ocean+1
        read(7,*) r_init(ir) , velc_i(:)
        velc(:,ir) = velc(:,ir) + velc_i(:) / (avrg_end-avrg_start)
      end do
      
      close(7)
    end do
    
    call out_spectra_sub('velc-averaged.spec', r_init, velc)
    
    deallocate( r_init, velc, velc_i )
    
  end subroutine save_spectra_velc_sub
  
end submodule Velc
