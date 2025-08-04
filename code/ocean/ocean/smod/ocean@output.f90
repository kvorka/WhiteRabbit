submodule (ocean) output
  implicit none; contains
  
  module procedure vypis_ocean_sub
    
    !!write(11,*) this%t, this%dt, this%nuss_fn(), this%reynolds_fn(choice='convective')
    write(11,*) this%t, this%dt, this%nuss_fn(), this%reynolds_fn(), this%reynolds_poloidal_fn(), &
              & this%reynolds_torroidal_fn(), this%temperature_fn()
    write(12,*) this%t, this%dt, this%laws_temp_fn(), this%laws_mech_fn()
    
    call this%vypis_sub(8, 'data/data_ocean_temp' , 'temperature')
    call this%vypis_sub(8, 'data/data_ocean_veloc', 'velocity'   )
    call this%vypis_sub(8, 'data/data_ocean_fluxu', 'fluxu'      )
    call this%vypis_sub(8, 'data/data_ocean_fluxd', 'fluxd'      )
    
    this%poc = this%poc + 1
    
  end procedure vypis_ocean_sub
  
end submodule output