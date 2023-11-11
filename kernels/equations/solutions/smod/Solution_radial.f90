submodule(Solution) Solution_radial
  implicit none ; contains
  
  module pure function temp_i_fn(this, ijm) result(temp)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ijm
    integer                        :: ir, is
    complex(kind=dbl), allocatable :: temp(:)
    
    allocate( temp(this%nd+1) ) ; temp = czero
    
    if ( this%inittemp ) then
      do concurrent ( ir = 1:this%nd+1 )
        is = 3*(ir-1)+1 ; temp(ir) = this%temp(is,ijm)
      end do
    end if
    
  end function temp_i_fn
  
  module pure function flux_i_fn(this, il, ijm) result(flux)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: il, ijm
    integer                        :: ir, is
    complex(kind=dbl), allocatable :: flux(:)
    
    allocate( flux(this%nd) ) ; flux = czero
    
    if ( this%inittemp ) then
      select case (il)
        case (-1)
          do concurrent ( ir = 1:this%nd )
            is = 3*(ir-1)+1 ; flux(ir) = this%temp(is+1,ijm)
          end do
        
        case ( 0)
          flux = czero
        
        case (+1)
          do concurrent ( ir = 1:this%nd )
            is = 3*(ir-1)+1 ; flux(ir) = this%temp(is+2,ijm)
          end do
      end select
    end if
    
  end function flux_i_fn
  
  module pure function velocity_i_fn(this, il, ijm) result(velocity)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: il, ijm
    integer                        :: ir, is
    complex(kind=dbl), allocatable :: velocity(:)
    
    allocate( velocity(this%nd+1) ) ; velocity = czero
    
    select case (il)
      case (-1)
        if ( this%initsfer ) then
          do concurrent ( ir = 1:this%nd+1 )
            is = 6*(ir-1)+1 ; velocity(ir) = this%mech(is,ijm)
          end do
        end if
      
      case (0)
        if ( this%inittorr ) then
          do concurrent ( ir = 1:this%nd+1 )
            is = 3*(ir-1)+1 ; velocity(ir) = this%torr(is,ijm)
          end do
        end if
      
      case (+1)
        if ( this%initsfer ) then
          do concurrent ( ir = 1:this%nd+1 )
            is = 6*(ir-1)+1 ; velocity(ir) = this%mech(is+1,ijm)
          end do
        end if
    end select
    
  end function velocity_i_fn
  
  module pure function deviatoric_stress_i_fn(this, il, ijm) result(dstress)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: il, ijm
    integer                        :: ir, is
    complex(kind=dbl), allocatable :: dstress(:)
    
    allocate( dstress(this%nd) ) ; dstress = czero
    
    if ( ijm >= 2 ) then
      select case (il)
        case (-2)
          if ( this%initsfer ) then
            do concurrent ( ir = 1:this%nd )
              is = 6*(ir-1)+1 ; dstress(ir) = this%mech(is+2,ijm)
            end do
          end if
        
        case (-1)
          if ( this%inittorr ) then
            do concurrent ( ir = 1:this%nd )
              is = 3*(ir-1)+1 ; dstress(ir) = this%torr(is+1,ijm)
            end do
          end if
        
        case ( 0)
          if ( this%initsfer ) then
            do concurrent ( ir = 1:this%nd )
              is = 6*(ir-1)+1 ; dstress(ir) = this%mech(is+4,ijm)
            end do
          end if
        
        case (+1)
          if ( this%inittorr ) then
            do concurrent ( ir = 1:this%nd )
              is = 3*(ir-1)+1 ; dstress(ir) = this%torr(is+2,ijm)
            end do
          end if
        
        case (+2)
          if ( this%initsfer ) then
            do concurrent ( ir = 1:this%nd )
              is = 6*(ir-1)+1 ; dstress(ir) = this%mech(is+5,ijm)
            end do
          end if
      end select
    end if
    
  end function deviatoric_stress_i_fn

end submodule Solution_radial
