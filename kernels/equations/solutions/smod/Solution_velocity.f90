submodule (Solution) Solution_velocity
  implicit none; contains
  
  module pure complex(kind=dbl) function velocity_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: isp, ist
    
    velocity_fn = czero
    
    if ( ijm >= 2 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      select case (il)
        case (-1)
          if ( this%initsfer ) velocity_fn = this%mech(isp,ijm)
        case ( 0)
          if ( this%inittorr ) velocity_fn = this%torr(ist,ijm)
        case (+1)
          if ( this%initsfer ) velocity_fn = this%mech(isp+1,ijm)
      end select
    end if
    
  end function velocity_fn
  
  module pure subroutine velocity_rr_il_many1_sub(this, il, ijm, velocity1)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: il, ijm
    complex(kind=dbl), intent(out) :: velocity1(:)
    integer                        :: ir
    
    do concurrent ( ir = 1:this%nd+1 )
      velocity1(ir) = czero
    end do
    
    select case (il)
      case (-1)
        if ( this%initsfer ) then
          do concurrent ( ir = 1:this%nd+1 )
            velocity1(ir) = this%mech(6*(ir-1)+1,ijm)
          end do
        end if
      
      case (0)
        if ( this%inittorr ) then
          do concurrent ( ir = 1:this%nd+1 )
            velocity1(ir) = this%torr(3*(ir-1)+1,ijm)
          end do
        end if
      
      case (+1)
        if ( this%initsfer ) then
          do concurrent ( ir = 1:this%nd+1 )
            velocity1(ir) = this%mech(6*(ir-1)+2,ijm)
          end do
        end if
      
    end select
    
  end subroutine velocity_rr_il_many1_sub
  
  module pure subroutine velocity_rr_many1_sub(this, ijm, velocity1)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ijm
    complex(kind=dbl), intent(out) :: velocity1(:,:)
    integer                        :: ir, isp, ist
    
    if ( (this%initsfer) .and. (this%inittorr) ) then
      do concurrent ( ir = 1:this%nd+1 )
        isp = 6*(ir-1)+1
        ist = 3*(ir-1)+1
        
        velocity1(1,ir) = this%mech(isp  ,ijm)
        velocity1(2,ir) = this%torr(ist  ,ijm)
        velocity1(3,ir) = this%mech(isp+1,ijm)
      end do
    
    else if (this%initsfer) then
      do concurrent ( ir = 1:this%nd+1 )
        isp = 6*(ir-1)+1
        
        velocity1(1,ir) = this%mech(isp  ,ijm)
        velocity1(2,ir) = czero
        velocity1(3,ir) = this%mech(isp+1,ijm)
      end do
    
    else if (this%inittorr) then
      do concurrent ( ir = 1:this%nd+1 )
        ist = 3*(ir-1)+1
        
        velocity1(1,ir) = czero
        velocity1(2,ir) = this%torr(ist  ,ijm)
        velocity1(3,ir) = czero
      end do
    
    end if
    
  end subroutine velocity_rr_many1_sub
  
end submodule Solution_velocity