submodule (ocean) init
  implicit none; contains
  
  module procedure init_ocean_sub
    
    call this%init_objects_sub( nd = nd_ocean, jmax = jmax_ocean, r_ud = r_ud_ocean, rgrid = grid_type_ocean, &
                              & gmod = gravity_ocean , g = (1-r_ud_ocean)**2 , noharm = noharm_ocean          )
    
    this%n_iter = n_iter_ocean
    this%cf     = 0.6_dbl
    this%ab     = 1.5_dbl
    
    this%Pr = Pr_ocean
    this%Ra = Ra_ocean
    this%Ek = Ek_ocean
    this%Cl = Kl_ocean
    
    this%D_ud         = D_ud_ocean
    this%rheology     = rheology_ocean
    this%mechanic_bnd = mechanic_bnd_ocean
    this%thermal_bnd  = thermal_bnd_ocean
    this%scaling      = scaling_ocean
    
    open(unit=11, file='data/Nuss.dat', status='new', action='write')
    open(unit=12, file='data/Laws.dat', status='new', action='write')
    
  end procedure init_ocean_sub
  
  module procedure deallocate_ocean_sub
    
    if ( allocated(this%nsph1) ) deallocate( this%nsph1 )
    if ( allocated(this%nsph2) ) deallocate( this%nsph2 )
    if ( allocated(this%ntorr) ) deallocate( this%ntorr )
    if ( allocated(this%ntemp) ) deallocate( this%ntemp )
    
    close(11); close(12)
    call this%deallocate_objects_sub()

  end procedure deallocate_ocean_sub
  
end submodule init