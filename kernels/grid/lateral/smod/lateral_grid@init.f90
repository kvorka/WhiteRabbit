submodule (lateral_grid) init
  implicit none ; contains
  
  module subroutine init_harmonics_sub(this, jmax)
    class(T_lateralGrid), intent(inout) :: this
    integer,              intent(in)    :: jmax
    
    if ( .not. ( any( addmissible_jmax == jmax ) ) ) then
      write(*,*) 'Due to FFT, this value of jmax is prohibited.'
      stop
    end if
    
    this%jmax = jmax+2
    
    this%nFourier  = 3 * ( this%jmax+1 )
    this%nLegendre = (3*(this%jmax+1)/2+1)/2+5-mod((3*(this%jmax+1)/2+1)/2+1,4)
    
    call legef_roots_sub( this%nLegendre, this%cosx )
    call legef_weights_sub( this%nLegendre, this%cosx, this%weight )
    
    call legep_factors_sub( this%jmax, this%amj, this%bmj, this%cmm )
    
    call this%reindexing%init_sub( this%jmax-2 )
    call this%fourtrans%init_sub( this%nFourier )
    
  end subroutine init_harmonics_sub
  
  module pure subroutine deallocate_harmonics_sub(this)
    class(T_lateralGrid), intent(inout) :: this
    
    call this%fourtrans%deallocate_sub()
    
    if ( allocated(this%cosx)   ) deallocate( this%cosx   )
    if ( allocated(this%weight) ) deallocate( this%weight )
    if ( allocated(this%amj)    ) deallocate( this%amj   )
    if ( allocated(this%bmj)    ) deallocate( this%bmj   )
    if ( allocated(this%cmm)    ) deallocate( this%cmm     )
    
  end subroutine deallocate_harmonics_sub
  
end submodule init
