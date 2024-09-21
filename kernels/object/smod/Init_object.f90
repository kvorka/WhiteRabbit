submodule (PhysicalObject) Init_object
  implicit none ; contains
  
  module subroutine init_objects_sub( this, nd, jmax, r_ud, rgrid, gmod, g, noobj, noharm )
    class(T_physicalObject),    intent(inout) :: this
    integer,                    intent(in)    :: nd, jmax
    real(kind=dbl),             intent(in)    :: r_ud, g
    character(len=*),           intent(in)    :: rgrid, gmod
    logical,          optional, intent(in)    :: noobj, noharm
    integer                                   :: j, m
    
    this%nd        = nd
    this%jmax      = jmax
    this%r_ud      = r_ud
    this%rd        = r_ud / ( 1-r_ud )
    this%ru        = 1    / ( 1-r_ud )
    this%grid_type = rgrid
    
    if (present(noobj)) then
      this%noobj = noobj
    else
      this%noobj = .false.
    end if
    
    if (present(noharm)) then
      this%noharm = noharm
    else if (this%noobj .eqv. .true.) then
      this%noharm = .true.
    else
      this%noharm = .false.
    end if
    
    this%jms  =   ( this%jmax*(this%jmax+1)/2+this%jmax ) + 1
    this%jms2 =   ( (this%jmax+2)*(this%jmax+3)/2+this%jmax+2 ) + 1
    this%jmv  = 3*( this%jmax*(this%jmax+1)/2+this%jmax ) + 1
    this%jmt  = 5*( this%jmax*(this%jmax+1)/2+this%jmax ) + 1
    
    if (.not. this%noobj) then
      call this%rad_grid%init_sub(this%nd, this%rd, this%ru, this%grid_type)
      if (.not. this%noharm) call this%lat_grid%init_sub(this%jmax)
      call this%sol%init_sub(this%nd, this%jmax)
      call this%mat%init_sub(this%nd, this%jmax, this%grid_type)
      call this%gravity%init_sub(gmod, g)
      call this%mparams%init_sub(this%nd, this%jmax)
      call this%bnd%init_sub(this%jms)

      allocate( this%j_indx(this%jms) )
        do j = 0, this%jmax
          do m = 0, j
            this%j_indx(j*(j+1)/2+m+1) = j
          end do
        end do
      
      this%gd = this%gravity%g_fn(this%rd)
      this%gu = this%gravity%g_fn(this%ru)
    end if
    
    this%poc = 0
    this%t   = zero
    
    call this%set_dt_sub()
    
  end subroutine init_objects_sub
  
  module subroutine set_dt_sub(this)
    class(T_physicalObject), intent(inout) :: this
    
    if (this%noobj) then
      call this%rad_grid%init_sub(this%nd, this%rd, this%ru, this%grid_type)
      this%dt  = 0.49_dbl * ( this%rad_grid%r(2)-this%rad_grid%r(1) )**2
      call this%rad_grid%deallocate_sub()
    else
      this%dt  = 0.49_dbl * ( this%rad_grid%r(2)-this%rad_grid%r(1) )**2
    end if
    
  end subroutine set_dt_sub
  
  module subroutine deallocate_objects_sub(this)
    class(T_physicalObject), intent(inout) :: this
    
    if ( allocated(this%j_indx)  ) deallocate( this%j_indx  )
    if ( allocated(this%flux_up) ) deallocate( this%flux_up )
    
    if ( allocated(this%rsph1) ) deallocate( this%rsph1 )
    if ( allocated(this%rsph2) ) deallocate( this%rsph2 )
    if ( allocated(this%rtorr) ) deallocate( this%rtorr )
    if ( allocated(this%rtemp) ) deallocate( this%rtemp )
    
    if ( allocated(this%nsph1) ) deallocate( this%nsph1 )
    if ( allocated(this%nsph2) ) deallocate( this%nsph2 )
    if ( allocated(this%ntorr) ) deallocate( this%ntorr )
    if ( allocated(this%ntemp) ) deallocate( this%ntemp )

    if (.not. this%noobj) then
      call this%sol%deallocate_sub()
      call this%mat%deallocate_sub()
      call this%rad_grid%deallocate_sub()
      call this%mparams%deallocate_sub()
      call this%bnd%deallocate_sub()
      call this%tdheat%deallocate_sub()
      
      if (.not. this%noharm) call this%lat_grid%deallocate_sub()
    end if
    
  end subroutine deallocate_objects_sub

end submodule Init_object