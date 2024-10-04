submodule (solution) velocity
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
  
  module pure subroutine conv_velocity_jml_sub(this, ir, velocity)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: velocity(:)
    integer                        :: ij, im, il, ijm
      
    call zero_carray_sub( this%jmv, velocity )
    
    do ij = 0, this%jmax
      do im = 1, ij
        ijm = ij*(ij+1)/2+im+1
        
        do il = abs(ij-1)-ij, +1
          velocity(3*(ijm-1)+il) = this%velocity_fn(ir, il, ijm)
        end do
      end do
    end do
    
  end subroutine conv_velocity_jml_sub
  
  module pure subroutine velocity_jml_many1_sub(this, ir, velocity1)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: velocity1(:)
    integer                        :: ijm, ijml, isp, ist
    
    call zero_carray_sub( this%jmv, velocity1 )
    
    if ( this%initsfer .and. this%inittorr ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = this%mech(isp,ijm)
        velocity1(ijml+1) = this%torr(ist,ijm)
        velocity1(ijml+2) = this%mech(isp+1,ijm)
      end do
    
    else if ( this%initsfer ) then
      isp = 6*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = this%mech(isp,ijm)
        velocity1(ijml+1) = czero
        velocity1(ijml+2) = this%mech(isp+1,ijm)
      end do
      
    else if ( this%inittorr ) then
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = czero
        velocity1(ijml+1) = this%torr(ist,ijm)
        velocity1(ijml+2) = czero
      end do
    end if
    
  end subroutine velocity_jml_many1_sub
  
  module pure subroutine velocity_jml_many2_sub(this, ir, velocity1, velocity2)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:)
    integer                        :: ijm, ijml, isp, ist
    
    call zero_carray_sub( this%jmv, velocity1 )
    call zero_carray_sub( this%jmv, velocity2 )
    
    if ( this%initsfer .and. this%inittorr ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = this%mech(isp ,ijm)
        velocity1(ijml+1) = this%torr(ist ,ijm)
        velocity1(ijml+2) = this%mech(isp+1,ijm)
        
        velocity2(ijml  ) = this%mech(isp+6,ijm)
        velocity2(ijml+1) = this%torr(ist+3,ijm)
        velocity2(ijml+2) = this%mech(isp+7,ijm)
      end do
    
    else if ( this%initsfer ) then
      isp = 6*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = this%mech(isp ,ijm)
        velocity1(ijml+2) = this%mech(isp+1,ijm)
        
        velocity2(ijml  ) = this%mech(isp+6,ijm)
        velocity2(ijml+2) = this%mech(isp+7,ijm)
      end do
      
    else if ( this%inittorr ) then
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml+1) = this%torr(ist  ,ijm)
        
        velocity2(ijml+1) = this%torr(ist+3,ijm)
      end do
    end if
    
  end subroutine velocity_jml_many2_sub
  
  module pure subroutine velocity_jml_many3_sub(this, ir, velocity1, velocity2, velocity3)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:), velocity3(:)
    integer                        :: ijm, ijml, isp, ist
    
    call zero_carray_sub( this%jmv, velocity1 )
    call zero_carray_sub( this%jmv, velocity2 )
    call zero_carray_sub( this%jmv, velocity3 )
    
    if ( this%initsfer .and. this%inittorr ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = this%mech(isp ,ijm)
        velocity1(ijml+1) = this%torr(ist ,ijm)
        velocity1(ijml+2) = this%mech(isp+1,ijm)
        
        velocity2(ijml  ) = this%mech(isp+6,ijm)
        velocity2(ijml+1) = this%torr(ist+3,ijm)
        velocity2(ijml+2) = this%mech(isp+7,ijm)
        
        velocity3(ijml  ) = this%mech(isp+12,ijm)
        velocity3(ijml+1) = this%torr(ist+ 6,ijm)
        velocity3(ijml+2) = this%mech(isp+13,ijm)
      end do
    
    else if ( this%initsfer ) then
      isp = 6*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml  ) = this%mech(isp ,ijm)
        velocity1(ijml+2) = this%mech(isp+1,ijm)
        
        velocity2(ijml  ) = this%mech(isp+6,ijm)
        velocity2(ijml+2) = this%mech(isp+7,ijm)
        
        velocity3(ijml  ) = this%mech(isp+12,ijm)
        velocity3(ijml+2) = this%mech(isp+13,ijm)
      end do
      
    else if ( this%inittorr ) then
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity1(ijml+1) = this%torr(ist  ,ijm)
        
        velocity2(ijml+1) = this%torr(ist+3,ijm)
        
        velocity3(ijml+1) = this%torr(ist+6,ijm)
      end do
    end if
    
  end subroutine velocity_jml_many3_sub
  
end submodule velocity