submodule (radial_grid) init
  implicit none ; contains
  
  module procedure init_grid_sub
    integer        :: i
    real(kind=dbl) :: dr
    
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
    
  end procedure init_grid_sub
  
  module procedure deallocate_grid_sub
    
    if ( allocated(this%r)  ) deallocate( this%r  )
    if ( allocated(this%rr) ) deallocate( this%rr )
    
  end procedure deallocate_grid_sub
  
end submodule init