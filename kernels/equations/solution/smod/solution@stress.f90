submodule (solution) stress
  implicit none; contains
  
  module procedure deviatoric_stress_fn
    integer :: isp, ist
    
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
    
  end procedure deviatoric_stress_fn
  
  module procedure deviatoric_stress_jml2_fn
    integer :: ij, im, il, ijm
    
    allocate( dstress(this%jmt) )
      call zero_carray_sub( this%jmt, dstress )
    
    do ij = 0, this%jmax
      do im = 0, ij
        ijm = ij*(ij+1)/2+im+1
        
        do il = abs(ij-2)-ij, +2
          dstress(5*(ijm-1)+il-1) = this%deviatoric_stress_fn(ir, il, ijm)
        end do
      end do
    end do
    
  end procedure deviatoric_stress_jml2_fn
  
end submodule stress