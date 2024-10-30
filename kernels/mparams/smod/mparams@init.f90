submodule (mparams) init
  implicit none; contains
  
  module procedure init_mparams_sub
    
    this%nd   = nd
    this%jmax = jmax
    this%jms  = jmax*(jmax+1)/2+jmax+1
    
    this%initvisc   = .false.
    this%initlambda = .false.
    this%initcp     = .false.
    this%initalpha  = .false.
    
  end procedure init_mparams_sub
  
  module procedure deallocate_mparams_sub
    
    if ( allocated(this%visc)   ) deallocate( this%visc   )
    if ( allocated(this%lambda) ) deallocate( this%lambda )
    if ( allocated(this%cp)     ) deallocate( this%cp     )
    if ( allocated(this%alpha)  ) deallocate( this%alpha  )
    
  end procedure deallocate_mparams_sub

end submodule init