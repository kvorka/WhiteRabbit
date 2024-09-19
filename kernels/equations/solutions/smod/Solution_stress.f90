submodule (Solution) Solution_stress
  implicit none; contains
  
  module pure complex(kind=dbl) function deviatoric_stress_fn(this, ir, il, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm, il
    integer                       :: isp, ist
    
    deviatoric_stress_fn = czero
    
    if ( ijm >= 2 ) then
      isp = 6*(ir-1)+1
      ist = 3*(ir-1)+1
      
      select case (il)
        case (-2)
          if ( this%initsfer ) deviatoric_stress_fn = this%mech(isp+2,ijm)
        case (-1)
          if ( this%inittorr ) deviatoric_stress_fn = this%torr(ist+1,ijm)
        case ( 0)
          if ( this%initsfer ) deviatoric_stress_fn = this%mech(isp+4,ijm)
        case (+1)
          if ( this%inittorr ) deviatoric_stress_fn = this%torr(ist+2,ijm)
        case (+2)
          if ( this%initsfer ) deviatoric_stress_fn = this%mech(isp+5,ijm)
      end select
    end if
    
  end function deviatoric_stress_fn
  
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
  
end submodule Solution_stress