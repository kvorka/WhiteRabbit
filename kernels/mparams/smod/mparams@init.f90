submodule (mparams) init
  implicit none; contains
  
  module subroutine init_mparams_sub(this, nd, jmax)
    class(T_Mparams), intent(inout) :: this
    integer,          intent(in)    :: nd, jmax
    
    this%nd   = nd
    this%jmax = jmax
    this%jms  = jmax*(jmax+1)/2+jmax+1
    
    this%initvisc   = .false.
    this%initlambda = .false.
    this%initcp     = .false.
    this%initalpha  = .false.
    
  end subroutine init_mparams_sub
  
  module subroutine deallocate_mparams_sub(this)
    class(T_Mparams), intent(inout) :: this
    
    if ( allocated(this%visc)   ) deallocate( this%visc   )
    if ( allocated(this%lambda) ) deallocate( this%lambda )
    if ( allocated(this%cp)     ) deallocate( this%cp     )
    if ( allocated(this%alpha)  ) deallocate( this%alpha  )
    
  end subroutine deallocate_mparams_sub

end submodule init