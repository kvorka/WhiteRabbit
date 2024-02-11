module OutputMod
  use Paths
  use Conversions
  implicit none; public; contains
  
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
  
  subroutine avrg_spectra_2d_sub(opt, spectra_out)
    character(len=*),  intent(in)  :: opt
    complex(kind=dbl), intent(out) :: spectra_out(:)
    integer                        :: in, ijm, error
    complex(kind=dbl)              :: spectrajm
    
    do in = avrg_start, avrg_end
      open( unit=7, file=opt//trim(adjustl(int2str_fn(in)))//'.dat', status='old', action='read' )
      
      do
        read(7,*,iostat=error) ijm, spectrajm
        
        if (error /= 0) then
          exit
        else
          spectra_out(ijm) = spectra_out(ijm) + spectrajm / (avrg_end-avrg_start)
        end if
      end do
      
      close(7)
    end do
    
  end subroutine avrg_spectra_2d_sub
  
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
  
end module OutputMod