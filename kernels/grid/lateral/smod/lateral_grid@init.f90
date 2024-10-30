submodule (lateral_grid) init
  implicit none ; contains
  
  module procedure init_harmonics_sub
    
    if ( .not. ( any( addmissible_jmax == jmax ) ) ) then
      write(*,*) 'Due to FFT, this value of jmax is prohibited.'
      stop
    end if
    
    call this%reindexing%init_sub( jmax )
    call this%fourtrans%init_sub( 3*(jmax+3) )
    call this%lgp%init_sub( jmax+2, (3*(jmax+3)/2+1)/2+5-mod((3*(jmax+3)/2+1)/2+1,4), 3._dbl*(jmax+3) )
    
  end procedure init_harmonics_sub
  
  module procedure deallocate_harmonics_sub
    
    call this%fourtrans%deallocate_sub()
    call this%lgp%deallocate_sub()
    
  end procedure deallocate_harmonics_sub
  
end submodule init
