submodule (lateral_grid) init
  implicit none ; contains
  
  module subroutine init_harmonics_sub(this, jmax)
    class(T_lateralGrid), intent(inout) :: this
    integer,              intent(in)    :: jmax
    
    if ( .not. ( any( addmissible_jmax == jmax ) ) ) then
      write(*,*) 'Due to FFT, this value of jmax is prohibited.'
      stop
    end if
    
    call this%reindexing%init_sub( jmax )
    call this%fourtrans%init_sub( 3*(jmax+3) )
    call this%lgp%init_sub( jmax+2, (3*(jmax+3)/2+1)/2+5-mod((3*(jmax+3)/2+1)/2+1,4) )
    
  end subroutine init_harmonics_sub
  
  module pure subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call this%fourtrans%deallocate_sub()
    call this%lgp%deallocate_sub()
    
  end subroutine deallocate_harmonics_sub
  
end submodule init
