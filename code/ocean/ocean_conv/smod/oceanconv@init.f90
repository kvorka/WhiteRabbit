submodule (oceanconv) init
  implicit none; contains
  
  module procedure init_oceanConv_sub
    integer           :: j, m, error
    complex(kind=dbl) :: flux
    
    call this%init_ocean_sub()
    
    if ( this%init_thermal_bnd ) then
      
      select case ( this%thermal_bnd )
        case ('fluxd')
          call this%bnd%init_flux_dn_sub()
          
          open(unit=35, file='code/ocean/heat_flux.cmplx', status='old', action='read')
            do
              read(35,*,iostat=error) j, m, flux
              
              if ( error == 0 ) then
                this%bnd%flux_dn(jm(j,m)) = flux
              else
                exit
              end if
            end do
          close(35)
          
          this%bnd%flux_dn = this%bnd%flux_dn / ( this%bnd%flux_dn(1)%re / s4pi )
      end select
    end if
    
    call this%init_eq_temp_sub()
      allocate( this%ntemp(this%jms,2:this%nd) )
      this%ntemp = czero
      
      call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    
    call this%init_eq_torr_sub()
      allocate( this%ntorr(this%jms,2:this%nd) )
      this%ntorr = czero
      
      call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_eq_mech_sub()
      allocate( this%nsph1(this%jms,2:this%nd), this%nsph2(this%jms,2:this%nd) )
      this%nsph1 = czero
      this%nsph2 = czero
      
      call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_state_sub()
    
  end procedure init_oceanConv_sub
  
end submodule init