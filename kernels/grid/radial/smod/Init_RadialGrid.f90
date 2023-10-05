submodule (RadialGrid) Init_RadialGrid
  implicit none
  
  contains
  
  subroutine init_grid_sub(this, nr, rd, ru, grid_type)
    class(T_radialGrid), intent(inout) :: this
    integer,             intent(in)    :: nr
    real(kind=dbl),      intent(in)    :: rd, ru
    character(len=5),    intent(in)    :: grid_type
    integer                            :: i
    real(kind=dbl)                     :: dr
    
    this%nd = nr; allocate( this%r(this%nd), this%rr(this%nd+1) )
    
    select case (grid_type)
      case('chebv')
        forall ( i=1:(this%nd  ) ) this%r(i)  = ( rd + ru ) / 2 - cos( (2*i-1) * pi / (2*this%nd) ) / cos( pi/(2*this%nd) ) / 2
        forall ( i=1:(this%nd+1) ) this%rr(i) = ( rd + ru ) / 2 - cos( (  i-1) * pi / (  this%nd) ) / cos( pi/(2*this%nd) ) / 2
      
      case('homog')
        dr = (ru-rd)/(this%nd-1)
          forall ( i=1:(this%nd  ) ) this%r(i)  = rd        + (i-1) * dr
          forall ( i=1:(this%nd+1) ) this%rr(i) = rd - dr/2 + (i-1) * dr
    end select
    
    this%volume = 4 * pi * ( ru**3 - rd**3 ) / 3
    
  end subroutine init_grid_sub
  
  subroutine deallocate_grid_sub(this)
    class(T_radialGrid), intent(inout) :: this
    
    if ( allocated(this%r)  ) deallocate( this%r  )
    if ( allocated(this%rr) ) deallocate( this%rr )
    
  end subroutine deallocate_grid_sub
  
end submodule Init_RadialGrid