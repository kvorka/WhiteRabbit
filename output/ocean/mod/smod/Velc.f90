submodule(OutputOceanMod) Velc
  implicit none
  
  contains
  
  subroutine harm_analysis_rad_velc_sub()
    integer,              parameter :: jmax = 100
    integer                         :: ir, iph, ith, ij, im, il
    real(kind=dbl),    allocatable  :: comp_x(:,:), comp_y(:,:), comp_z(:,:), comp_r(:,:), comp_t(:,:), comp_p(:,:), rvelc(:,:)
    complex(kind=dbl), allocatable  :: x(:), y(:), z(:), velc120(:,:)
    
    allocate( velc120(jmv,n_out), x(jms1), y(jms1), z(jms1), comp_x(2*nth,nth), comp_y(2*nth,nth), comp_z(2*nth,nth), &
            & comp_r(2*nth,nth), comp_t(2*nth,nth), comp_p(2*nth,nth), rvelc(nth,n_out)                               )
    
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
      
      call harmsy_sub(jmax, x, comp_x)
      call harmsy_sub(jmax, y, comp_y)
      call harmsy_sub(jmax, z, comp_z)
      
      do ith = 1, nth
        do iph = 1, 2*nth
          call vecxyz2rtp_sub( ith-0.5_dbl, iph-1._dbl, comp_x(iph,ith), comp_y(iph,ith), comp_z(iph,ith), &
                             &                          comp_r(iph,ith), comp_t(iph,ith), comp_p(iph,ith)  )
        end do
      end do
      
      call get_zondata_sub(comp_r, rvelc(:,ir))
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
            & comp_r(2*nth,nth), comp_t(2*nth,nth), comp_p(2*nth,nth), zvelc(nth,n_out)                               )
    
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
      
      call harmsy_sub(jmax, x, comp_x)
      call harmsy_sub(jmax, y, comp_y)
      call harmsy_sub(jmax, z, comp_z)
      
      do ith = 1, nth
        do iph = 1, 2*nth
          call vecxyz2rtp_sub( ith-0.5_dbl, iph-1._dbl, comp_x(iph,ith), comp_y(iph,ith), comp_z(iph,ith), &
                             &                          comp_r(iph,ith), comp_t(iph,ith), comp_p(iph,ith)  )
        end do
      end do
      
      call get_zondata_sub(comp_p, zvelc(:,ir))
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
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: velc(:,:)
    
    allocate( r(nd_ocean+1), velc(jmv,nd_ocean+1) ) ; velc = czero
    
    call avrg_spectra_3d_sub(path_ocean_velc, r, velc)
    call out_spectra_3d_sub('velc-averaged.spec', r, velc)
    
    deallocate( r, velc )
    
  end subroutine save_spectra_velc_sub
  
end submodule Velc
