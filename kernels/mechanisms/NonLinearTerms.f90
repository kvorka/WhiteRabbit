module NonLinearTerms
  use PhysicalObject
  implicit none
  
  public :: vgradT_fn
  public :: vgradv_fn
  public :: fullnl_sub
  public :: erT_fn
  public :: coriolis_fn
  
  contains
  
  function vgradT_fn(this, i) result(vgradT)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: vgradT(this%jms)
    
    vgradT = -this%lat_grid%vcvv_fn( this%sol%velocity_jml_fn(i), this%mgradT_rrjml_fn(i) )
  
  end function vgradT_fn
   
  function vgradv_fn(this, i) result(vgradv)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: vgradv(this%jmv)
    complex(kind=dbl),      allocatable :: v(:)
  
    allocate( v(this%jmv) ) ; v = this%sol%velocity_jml_fn(i)
      
      vgradv = this%lat_grid%vcsv_vcvgv_fn(this%rad_grid%rr(i), this%dv_dr_rrjml_fn(i,v), v)
  
    deallocate(v)
  
  end function vgradv_fn
  
  function erT_fn(this, i) result(erT)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: i
    complex(kind=dbl)                   :: erT(this%jmv)
    
    erT = ersv_fn( this%jmax, this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i)) * this%sol%temp_jm_fn(i) )
  
  end function erT_fn
  
  function coriolis_fn(this, i, v) result(coriolis)
    class(T_physicalObject),     intent(in) :: this
    integer,           optional, intent(in) :: i                   !priestorova diskretizacia
    complex(kind=dbl), optional, intent(in) :: v(:)
    complex(kind=dbl)                       :: coriolis(this%jmv)  !Coriolisova sila
    
    if (present(i)) coriolis = ezvv_fn( this%jmax, this%sol%velocity_jml_fn(i) )
    if (present(v)) coriolis = ezvv_fn( this%jmax, v )
  
  end function coriolis_fn
  
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

      nlm = nlm / this%Pr - this%Ra * erT_fn(this,i) + coriolis_fn(this=this, v=v) * 2 / this%Ek
    
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