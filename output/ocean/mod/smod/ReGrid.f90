submodule(OutputOceanMod) ReGrid
  implicit none

  contains

  subroutine get_zonal_sub(data_in, data_out)
    real(kind=dbl), intent(in)  :: data_in(:,:)
    real(kind=dbl), intent(out) :: data_out(:)
    integer                     :: k

    data_out(1) = data_in(1,1)

    do k = 1, nth-1
      data_out(k+1) = ( data_in(1,k) + data_in(1,k+1) ) / 2
    end do

    data_out(nth+1) = data_in(1,nth)
    
  end subroutine get_zonal_sub

  subroutine load_data_sub(dim_in, file_in, data_out)
    character(len=*),  intent(in)  :: file_in
    integer,           intent(in)  :: dim_in
    complex(kind=dbl), intent(out) :: data_out(:,:)
    real(kind=dbl),    allocatable :: r_in(:)
    complex(kind=dbl), allocatable :: data_in(:,:)
    integer                        :: i, ii
    real(kind=dbl)                 :: r

    allocate( r_in(nd_ocean+1), data_in(dim_in, nd_ocean+1) )

      open(unit=1, file=file_in, status='old', action='read')
        do i = 1, nd_ocean+1
          read(1,*) r_in(i), data_in(:,i)
        end do
      close(1)

      do i = 1, n_out
        r = r_ud_ocean/(1-r_ud_ocean) + (i-1._dbl)/(n_out-1._dbl)

        do ii = 1, nd_ocean
          if ( (r >= r_in(ii)) .and. (r <= r_in(ii+1)) ) then
            data_out(:,i) = ( (r-r_in(ii)) * data_in(:,ii+1) + (r_in(ii+1)-r) * data_in(:,ii) ) / ( r_in(ii+1)-r_in(ii) )
            exit
          end if
        end do
      end do
    
    deallocate( r_in, data_in )
    
    data_out(1,:) = cmplx(0._dbl, 0._dbl, kind=dbl)

  end subroutine load_data_sub

  subroutine save_data_sub(file_in, data_in)
    character(len=*), intent(in) :: file_in
    real(kind=dbl),   intent(in) :: data_in(:,:)
    integer                      :: i, k

    open(unit=2, file=file_in, status='new', action='write')
    
      do k = 0, nth
        do i = 1, n_out
          write(2,*) k+90, 480+i, data_in(k+1,i)
        end do
      end do

    close(2)

  end subroutine save_data_sub

end submodule ReGrid