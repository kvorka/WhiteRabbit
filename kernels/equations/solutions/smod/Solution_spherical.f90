submodule(Solution) Solution_spherical
  implicit none ; contains
  
  module pure function velocity_jml_fn(this, ir) result(velocity)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), allocatable :: velocity(:)
    integer                        :: ijm, ijml, isp, ist
    
    allocate( velocity(this%jmv) ) ; velocity = czero
    
    if ( this%inittorr .and. this%initsfer ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity(ijml  ) = this%mech(isp,ijm)
        velocity(ijml+1) = this%torr(ist,ijm)
        velocity(ijml+2) = this%mech(isp+1,ijm)
      end do
      
    else if ( this%initsfer ) then
      isp = 6*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity(ijml  ) = this%mech(isp,ijm)
        velocity(ijml+1) = czero
        velocity(ijml+2) = this%mech(isp+1,ijm)
      end do
      
    else if ( this%inittorr ) then
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity(ijml  ) = czero
        velocity(ijml+1) = this%torr(ist,ijm)
        velocity(ijml+2) = czero
      end do
    end if
    
  end function velocity_jml_fn
  
  module pure function conv_velocity_jml_fn(this, ir) result(velocity)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), allocatable :: velocity(:)
    integer                        :: ij, im, il, ijm
      
    allocate( velocity(this%jmv) ) ; velocity = czero
    
    do ij = 0, this%jmax
      do im = 1, ij
        ijm = ij*(ij+1)/2+im+1
        
        do il = abs(ij-1)-ij, +1
          velocity(3*(ijm-1)+il) = this%velocity_fn(ir, il, ijm)
        end do
      end do
    end do
    
  end function conv_velocity_jml_fn
  
  module pure function deviatoric_stress_jml2_fn(this, ir) result(dstress)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), allocatable :: dstress(:)
    integer                        :: ij, im, il, ijm
    
    allocate( dstress(this%jmt) ) ; dstress = czero
    
    do ij = 0, this%jmax
      do im = 0, ij
        ijm = ij*(ij+1)/2+im+1
        
        do il = abs(ij-2)-ij, +2
          dstress(5*(ijm-1)+il-1) = this%deviatoric_stress_fn(ir, il, ijm)
        end do
      end do
    end do
    
  end function deviatoric_stress_jml2_fn
  
  module pure subroutine velocity_jml_sub(this, ir, velocity)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: velocity(:)
    integer                        :: ijm, ijml, isp, ist
    
    velocity = czero
    
    if ( this%initsfer .and. this%inittorr ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity(ijml  ) = this%mech(isp,ijm)
        velocity(ijml+1) = this%torr(ist,ijm)
        velocity(ijml+2) = this%mech(isp+1,ijm)
      end do
    
    else if ( this%initsfer ) then
      isp = 6*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity(ijml  ) = this%mech(isp,ijm)
        velocity(ijml+1) = czero
        velocity(ijml+2) = this%mech(isp+1,ijm)
      end do
      
    else if ( this%inittorr ) then
      ist = 3*(ir-1)+1
      
      do concurrent ( ijm = 2:this%jms )
        ijml = 3*(ijm-1)-1
        
        velocity(ijml  ) = czero
        velocity(ijml+1) = this%torr(ist,ijm)
        velocity(ijml+2) = czero
      end do
    end if
    
  end subroutine velocity_jml_sub
  
  module pure subroutine velocity_jml_many_sub(this, ir, velocity1, velocity2, velocity3)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: velocity1(:), velocity2(:), velocity3(:)
    integer                        :: ijm, ijml, isp, ist
    
    isp = 6*(ir-1)+1
    ist = 3*(ir-1)+1
    
    ijml = 1
      velocity1(ijml) = czero
      velocity2(ijml) = czero
      velocity3(ijml) = czero
    
    do concurrent ( ijm = 2:this%jms )
      ijml = 3*(ijm-1)-1
      
      velocity1(ijml  ) = this%mech(isp,ijm)
      velocity1(ijml+1) = this%torr(ist,ijm)
      velocity1(ijml+2) = this%mech(isp+1,ijm)
      
      velocity2(ijml  ) = this%mech(isp+6,ijm)
      velocity2(ijml+1) = this%torr(ist+3,ijm)
      velocity2(ijml+2) = this%mech(isp+7,ijm)
      
      velocity3(ijml  ) = this%mech(isp+12,ijm)
      velocity3(ijml+1) = this%torr(ist+ 6,ijm)
      velocity3(ijml+2) = this%mech(isp+13,ijm)
    end do
    
  end subroutine velocity_jml_many_sub
  
end submodule Solution_spherical