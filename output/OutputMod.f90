module OutputMod
  use Paths
  use Harmsy
  use Vector_analysis
  implicit none
  
  public :: get_zondata_sub
  public :: out_data_2d_sub
  public :: out_data_1d_sub
  public :: out_spectra_2d_sub
  public :: out_spectra_3d_sub
  public :: avrg_spectra_2d_sub
  public :: avrg_spectra_3d_sub
  
  contains

  subroutine get_zondata_sub(data_in, data_out)
    real(kind=dbl), intent(in)  :: data_in(:,:)
    real(kind=dbl), intent(out) :: data_out(:)
    integer                     :: ith
    
    data_out(1) = sum( data_in(:,1) ) / ( 2 * nth )
    
    do ith = 1, nth-1
      data_out(ith+1) = sum( data_in(:,ith) + data_in(:,ith+1) ) / ( 2 * nth ) / 2
    end do
    
    data_out(nth+1) = sum ( data_in(:,nth) ) / ( 2 * nth )
    
  end subroutine get_zondata_sub
  
  subroutine out_data_1d_sub(opt, data_in)
    character(len=*), intent(in) :: opt
    real(kind=dbl),   intent(in) :: data_in(:,:)
    integer                      :: ith
      
    open(unit=1, file=opt, status='new', action='write')
      
    ith = 1
      write(1,*) -90, sum( data_in(:,ith) ) / ( 2 * nth )
    
    do ith = 2, nth
      write(1,*) ith-91, sum( data_in(:,ith) + data_in(:,ith-1) ) / ( 2 * nth ) / 2
    end do

    ith = nth
      write(1,*) +90, sum( data_in(:,ith) ) / ( 2 * nth )
      
    close(1)

  end subroutine out_data_1d_sub
  
  subroutine out_data_2d_sub(opt, data_in)
    character(len=*), intent(in) :: opt
    real(kind=dbl),   intent(in) :: data_in(:,:)
    integer                      :: ith, iph
    
    open(unit=1, file=opt, status='new', action='write')
    
    ith = nth
      do iph = 2, 2*nth
        write(1,*) iph-1, +90, data_in(iph,ith)
      end do
    
      iph = 1
        write(1,*) 360, +90, data_in(iph,ith)
    
    do ith = 2, nth
      do iph = 2, 2*nth
        write(1,*) iph-1, ith-91, ( data_in(iph,ith) + data_in(iph,ith-1) ) / 2
      end do
      
      iph = 1
        write(1,*) 360, ith-91, ( data_in(iph,ith) + data_in(iph,ith-1) ) / 2
    end do
    
    ith = 1
      do iph = 2, 2*nth
        write(1,*) iph-1, -90, data_in(iph,ith)
      end do
      
      iph = 1
        write(1,*) 360, -90, data_in(iph,ith)
      
    close(1)
    
  end subroutine out_data_2d_sub
  
  subroutine avrg_spectra_2d_sub(opt, spectra_out)
    character(len=*),  intent(in)  :: opt
    complex(kind=dbl), intent(out) :: spectra_out(:)
    integer                        :: in, ij, im, error
    complex(kind=dbl)              :: spectrajm
    
    do in = avrg_start, avrg_end
      open( unit=7, file=opt//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read' )
      
      do
        read(7,*,iostat=error) ij, im, spectrajm
        
        if (error /= 0) then
          exit
        else
          spectra_out(jm(ij,im)) = spectra_out(jm(ij,im)) + spectrajm / (avrg_end-avrg_start)
        end if
      end do
      
      close(7)
    end do
    
  end subroutine avrg_spectra_2d_sub
  
  subroutine out_spectra_2d_sub(opt, spectra_in)
    character(len=*),     intent(in) :: opt
    complex(kind=dbl), intent(inout) :: spectra_in(:)
    integer                          :: ijm
    
    open(unit=1, file=opt, status='new', action='write')
    
    do ijm = 1, size(spectra_in)
      write(1,*) ijm, spectra_in(ijm)
    end do
    
    close(1)
    
  end subroutine out_spectra_2d_sub
  
  subroutine avrg_spectra_3d_sub(opt, r_out, spectra_out)
    character(len=*),  intent(in)  :: opt
    real(kind=dbl),    intent(out) :: r_out(:)
    complex(kind=dbl), intent(out) :: spectra_out(:,:)
    integer                        :: in, ir
    complex(kind=dbl), allocatable :: spectra(:)
    
    allocate( spectra(size(spectra_out, dim=1)) )
    
    do in = avrg_start, avrg_end
      open(unit=7, file=opt//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read')
      
      do ir = 1, size(r_out)
        read(7,*) r_out(ir) , spectra(:)
        spectra_out(:,ir) = spectra_out(:,ir) + spectra(:) / (avrg_end-avrg_start)
      end do
      
      close(7)
    end do
    
    deallocate( spectra )
    
  end subroutine avrg_spectra_3d_sub
  
  subroutine out_spectra_3d_sub(opt, r_in, data_in)
    character(len=*),  intent(in) :: opt
    real(kind=dbl),    intent(in) :: r_in(:)
    complex(kind=dbl), intent(in) :: data_in(:,:)
    integer                       :: ir
    
    open(unit=1, file=opt, status='new', action='write')
    
    do ir = 1, size(r_in)
      write(1,*) r_in(ir), data_in(:,ir)
    end do
    
    close(1)
    
  end subroutine out_spectra_3d_sub

end module OutputMod