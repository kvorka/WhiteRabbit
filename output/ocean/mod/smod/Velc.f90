submodule(OutputOceanMod) Velc
  implicit none

  contains

  subroutine harm_analysis_rad_velc_sub(path)
    integer,              parameter :: jmax = 100
    character(len=*),    intent(in) :: path
    integer                         :: i, k, ii, j, m, l
    real(kind=dbl),    allocatable  :: comp_x(:,:), comp_y(:,:), comp_z(:,:), comp_r(:,:), comp_t(:,:), comp_p(:,:), t(:,:)
    complex(kind=dbl), allocatable  :: x(:), y(:), z(:), velc120(:,:)

    allocate( velc120(jmv,n_out), x(jms1), y(jms1), z(jms1), comp_x(2*nth,nth), comp_y(2*nth,nth), comp_z(2*nth,nth), &
            & comp_r(2*nth,nth), comp_t(2*nth,nth), comp_p(2*nth,nth), t(0:nth,n_out)                                )

      call load_data_sub(jmv, 'velc-averaged.spec', velc120)

      do j = 0, jmax_ocean
        do m = 0, j
          do l = abs(j-1), j+1
            if( (mod(j,2) /= 0) .or. (m /= 0) .or. (l == j) ) velc120(3*(j*(j+1)/2+m)+l-j,:) = cmplx(0._dbl, 0._dbl, kind=dbl)
          end do
        end do
      end do

      do i = 1, n_out
        call vec2scals_sub(jmax_ocean, velc120(:,i), x, y, z); y = -y

        call toGrid_sub(jmax, x, comp_x)
        call toGrid_sub(jmax, y, comp_y)
        call toGrid_sub(jmax, z, comp_z)

        do ii = 1, nth
          do k = 1, 2*nth
            call vecxyz2rtp_sub( 180.5_dbl-ii, k-1._dbl, comp_x(k,ii), comp_y(k,ii), comp_z(k,ii), &
                               &                         comp_r(k,ii), comp_t(k,ii), comp_p(k,ii)  )
          end do
        end do

        call get_zonal_sub(comp_r, t(:,i))
      end do

    deallocate( velc120, x, y, z, comp_x, comp_y, comp_z, comp_r, comp_t, comp_p )

      t = Ek_ocean * t / Pr_ocean
        open(unit=8, file='ocean_rv_max', status='new', action='write'); write(8,'(1g9.2e1)') maxval(t); close(8)
        open(unit=8, file='ocean_rv_min', status='new', action='write'); write(8,'(1g9.2e1)') minval(t); close(8)

      t = t / max( maxval(t), abs(minval(t)) )
        open(unit=8, file='ocean_rv_max_g', status='new', action='write'); write(8,'(1f4.2)') min(maxval(t), +treshold); close(8)
        open(unit=8, file='ocean_rv_min_g', status='new', action='write'); write(8,'(1f5.2)') max(minval(t), -treshold); close(8)

        call save_data_sub('inp-r.dat', t)

    deallocate( t )
  
  end subroutine harm_analysis_rad_velc_sub

  subroutine harm_analysis_zon_velc_sub(path)
    integer,              parameter :: jmax = 100
    character(len=*),    intent(in) :: path
    integer                         :: i, k, ii, j, m, l
    real(kind=dbl),    allocatable  :: comp_x(:,:), comp_y(:,:), comp_z(:,:), comp_r(:,:), comp_t(:,:), comp_p(:,:), t(:,:)
    complex(kind=dbl), allocatable  :: x(:), y(:), z(:), velc120(:,:)

    allocate( velc120(jmv,n_out), x(jms1), y(jms1), z(jms1), comp_x(2*nth,nth), comp_y(2*nth,nth), comp_z(2*nth,nth), &
            & comp_r(2*nth,nth), comp_t(2*nth,nth), comp_p(2*nth,nth), t(0:nth,n_out)                                 )

      call load_data_sub(jmv, 'velc-averaged.spec', velc120)

      do j = 0, jmax_ocean
        do m = 0, j
          do l = abs(j-1), j+1
            if( (mod(j,2) == 0) .or. (m /= 0) .or. (l /= j) ) velc120(3*(j*(j+1)/2+m)+l-j,:) = cmplx(0._dbl, 0._dbl, kind=dbl)
          end do
        end do
      end do
  
      do i = 1, n_out
        call vec2scals_sub(jmax_ocean, velc120(:,i), x, y, z); y = -y
    
        call toGrid_sub(jmax, x, comp_x)
        call toGrid_sub(jmax, y, comp_y)
        call toGrid_sub(jmax, z, comp_z)
    
        do ii = 1, nth
          do k = 1, 2*nth
            call vecxyz2rtp_sub( 180.5_dbl-ii, k-1._dbl, comp_x(k,ii), comp_y(k,ii), comp_z(k,ii), &
                               &                         comp_r(k,ii), comp_t(k,ii), comp_p(k,ii)  )
          end do
        end do
    
        call get_zonal_sub(comp_p, t(:,i))
      end do
  
    deallocate( velc120, x, y, z, comp_x, comp_y, comp_z, comp_r, comp_t, comp_p )

      t = Ek_ocean * t / Pr_ocean
        open(unit=8, file='ocean_zv_max', status='new', action='write'); write(8,'(1g9.2e1)') maxval(t); close(8)
        open(unit=8, file='ocean_zv_min', status='new', action='write'); write(8,'(1g9.2e1)') minval(t); close(8)

      t = t / max( maxval(t), abs(minval(t)) )
        open(unit=8, file='ocean_zv_max_g', status='new', action='write'); write(8,'(1f4.2)') min(maxval(t), +treshold); close(8)
        open(unit=8, file='ocean_zv_min_g', status='new', action='write'); write(8,'(1f5.2)') max(minval(t), -treshold); close(8)

      call save_data_sub('inp-p.dat', t)

    deallocate( t )

  end subroutine harm_analysis_zon_velc_sub

  subroutine save_spectra_velc_sub(path)
    character(len=*), intent(in)   :: path
    character(len=10)              :: subor
    integer                        :: n, i
    real(kind=dbl),    allocatable :: r_init(:), r_out(:)
    complex(kind=dbl), allocatable :: velc(:,:), velc_i(:), velc_out(:,:)

    allocate( r_init(nd_ocean+1), velc(jmv,nd_ocean+1), velc_i(jmv) )
      velc = cmplx(0._dbl, 0._dbl, kind=dbl)

      do n = i1, i2
        write(subor,'(1I4)') n
    
        open( unit=7, file=path//trim(adjustl(subor))//'.dat', status='old', action='read' )
          do i = 1, nd_ocean+1
            read(7,*) r_init(i) , velc_i(:)
            velc(:,i) = velc(:,i) + velc_i(:) / (i2-i1)
          end do
        close(7)
      end do

      call out_spectra_sub('velc-averaged.spec', r_init, velc)
    deallocate( r_init, velc, velc_i )

  end subroutine save_spectra_velc_sub

end submodule Velc
