module NonLinearTerms
  use Math
  use PhysicalObject
  implicit none
  
  public :: vgradT_fn
  public :: vgradv_fn
  public :: fullnl_sub
  public :: erT_fn
  public :: coriolis_fn
  
  contains
  
  function vgradT_fn(this, i) result(vgradT)
    class(T_physicalObject),    intent(in) :: this
    integer,                    intent(in) :: i
    complex(kind=dbl), dimension(this%jms) :: vgradT
    
    vgradT = -this%lat_grid%vcvv_fn( this%sol%velocity_jml_fn(i),                                               &
                                   & this%rad_grid%cc(i,-1) * this%sol%flux_jml_fn(i-1) / this%lambda_fn(i-1) + &
                                   & this%rad_grid%cc(i,+1) * this%sol%flux_jml_fn(i  ) / this%lambda_fn(i)     )
  
  end function vgradT_fn
   
  function vgradv_fn(this, i) result(vgradv)
    class(T_physicalObject),          intent(in) :: this
    integer,                          intent(in) :: i       !priestorova diskretizacia
    complex(kind=dbl), dimension(this%jmv)       :: vgradv  !v \cdot grad(v)
    complex(kind=dbl), dimension(:), allocatable :: dv, v   !radialna derivacia rychlosti, rychlost
  
    associate( jmax     => this%jmax,       &
             & jmv      => this%jmv,        &
             & sol      => this%sol,        &
             & lat_grid => this%lat_grid,   &
             & r        => this%rad_grid%r, &
             & rr       => this%rad_grid%rr )
  
    allocate( dv(jmv), v(jmv) )
      v  = sol%velocity_jml_fn(i)
  
      dv = ( (rr(i+1)-rr(i))/(rr(i-1)-rr(i)) * ( sol%velocity_jml_fn(i-1)-sol%velocity_jml_fn(i) ) - &
           & (rr(i-1)-rr(i))/(rr(i+1)-rr(i)) * ( sol%velocity_jml_fn(i+1)-sol%velocity_jml_fn(i) )   ) / (rr(i+1)-rr(i-1))
      
      vgradv = lat_grid%vcsv_fn(ervs_fn(jmax, v), dv) + lat_grid%vcvgv_fn(v, v) / rr(i)
  
    deallocate(dv, v)
    
    end associate
  
  end function vgradv_fn
  
  function erT_fn(this, i) result(erT)
    class(T_physicalObject),    intent(in) :: this
    integer,                    intent(in) :: i      !priestorova diskretizacia
    complex(kind=dbl), dimension(this%jmv) :: erT    !er*T
    
    erT = ersv_fn( this%jmax, this%alpha_fn(i) * this%gravity%g_fn(this%rad_grid%rr(i)) * this%sol%temp_jm_fn(i) )
  
  end function erT_fn
  
  function coriolis_fn(this, i) result(coriolis)
    class(T_physicalObject),    intent(in) :: this
    integer,                    intent(in) :: i         !priestorova diskretizacia
    complex(kind=dbl), dimension(this%jmv) :: coriolis  !Coriolisova sila
    
    coriolis = ezvv_fn( this%jmax, this%sol%velocity_jml_fn(i) )
  
  end function coriolis_fn
  
  subroutine fullnl_sub(this, i, nlscal, nlvect)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: i
    complex(kind=dbl),       intent(out) :: nlscal(:), nlvect(:)
    complex(kind=dbl),       allocatable :: dv(:), v(:), q(:)
  
    associate( rr => this%rad_grid%rr ) ; allocate( dv(this%jmv), v(this%jmv), q(this%jmv) )

      v  = this%sol%velocity_jml_fn(i)

      q  = this%rad_grid%cc(i,-1) * this%sol%flux_jml_fn(i-1) / this%lambda_fn(i-1) + &
         & this%rad_grid%cc(i,+1) * this%sol%flux_jml_fn(i  ) / this%lambda_fn(i)

      dv = ( (rr(i+1)-rr(i))/(rr(i-1)-rr(i)) * ( this%sol%velocity_jml_fn(i-1) - v ) - &
           & (rr(i-1)-rr(i))/(rr(i+1)-rr(i)) * ( this%sol%velocity_jml_fn(i+1) - v )   ) / (rr(i+1)-rr(i-1))

      call this%lat_grid%vcsv_vcvv_vcvgv_sub(rr(i), q, dv, v, nlscal, nlvect)
  
    deallocate(dv, v, q) ; end associate

    nlvect = nlvect / this%Pr - this%Ra * erT_fn(this,i) + coriolis_fn(this,i) * 2/this%Ek
    
  end subroutine fullnl_sub

end module NonLinearTerms