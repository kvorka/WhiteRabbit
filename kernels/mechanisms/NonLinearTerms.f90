module NonLinearTerms
  use PhysicalObject
  implicit none
  
  public :: vgradT_fn
  public :: vgradv_fn
  public :: fullnl_sub
  
  contains
  
  function vgradT_fn(this, i) result(vgradT)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: vgradT(this%jms)
    complex(kind=dbl),      allocatable :: v(:), mgradT(:)
    
    allocate( v(this%jmv), mgradT(this%jmv) )

      v = this%sol%velocity_jml_fn(i) ; mgradT = this%mgradT_rrjml_fn(i)

      vgradT = -this%lat_grid%vcvv_fn( v, mgradT )
  
    deallocate( v, mgradT )

  end function vgradT_fn
   
  function vgradv_fn(this, i) result(vgradv)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: vgradv(this%jmv)
    complex(kind=dbl),      allocatable :: v(:), dv(:)
  
    allocate( v(this%jmv), dv(this%jmv) )
      
      v  = this%sol%velocity_jml_fn(i)
      dv = this%dv_dr_rrjml_fn(i,v)
      
      vgradv = this%lat_grid%vcsv_vcvgv_fn(this%rad_grid%rr(i), dv, v)
  
    deallocate( v, dv )
  
  end function vgradv_fn
  
  subroutine fullnl_sub(this, i, ntemp, nsph1, ntorr, nsph2)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: i
    complex(kind=dbl),       intent(out) :: ntemp(:), nsph1(:), ntorr(:), nsph2(:)
    integer                              :: ijm, ijml
    complex(kind=dbl),       allocatable :: v(:), nlm(:), dv(:), dT(:)

    allocate( v(this%jmv), dv(this%jmv), dT(this%jmv) )
    
      v  = this%sol%velocity_jml_fn(i)
      dv = this%dv_dr_rrjml_fn(i,v)
      dT = this%mgradT_rrjml_fn(i)
    
    allocate( nlm(this%jmv) )
      
      call this%lat_grid%vcsv_vcvv_vcvgv_sub(this%rad_grid%rr(i), dT, dv, v, ntemp, nlm)

    deallocate( dv, dT )

      nlm = nlm / this%Pr - this%buoy_rr_jml_fn(i) + this%coriolis_rr_jml_fn(v=v)
    
    deallocate( v )
      
      do ijm = 2, this%jms
        ijml = 3*(ijm-1)-1

        nsph1(ijm) = nlm(ijml  )
        ntorr(ijm) = nlm(ijml+1)
        nsph2(ijm) = nlm(ijml+2)
      end do
      
    deallocate( nlm )
    
  end subroutine fullnl_sub

end module NonLinearTerms