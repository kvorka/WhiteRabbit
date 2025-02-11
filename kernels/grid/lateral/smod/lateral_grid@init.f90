submodule (lateral_grid) init
  implicit none ; contains
  
  module procedure init_harmonics_sub
    integer :: nL, nF
    
    if ( .not. ( any( addmissible_jmax == jmax ) ) ) then
      write(*,*) 'Due to FFT, this value of jmax is prohibited.'
      stop
    end if
    
    !Prepare indexing object
    call this%reindexing%init_sub( jmax )
    
    !FFT
    nF = 3*(jmax+3)
    call this%fourtrans%init_sub( nF )
    
    !Sums of associated Legendre polynomials
    nL = (3*(jmax+2)/2+1)/2+step+1-mod((3*(jmax+2)/2+1)/2+1,step)
    call this%lgp%init_sub( jmax+2, nL, i2r_fn( nF ) )
    
  end procedure init_harmonics_sub
  
  module procedure deallocate_harmonics_sub
    
    call this%fourtrans%deallocate_sub()
    call this%lgp%deallocate_sub()
    
  end procedure deallocate_harmonics_sub
  
end submodule init