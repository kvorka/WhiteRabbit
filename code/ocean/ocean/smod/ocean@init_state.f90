submodule (ocean) init_state
  implicit none; contains
  
  subroutine init_conduction_sub(this)
    class(T_ocean), intent(inout) :: this
    integer                       :: is, ir, ij, im
    real(kind=dbl)                :: realp, imagp, dt_help, ab_help, cf_help
    
    !! Set the time-step to infinity and method to fully implicit
    dt_help = this%dt
    ab_help = this%ab
    cf_help = this%cf
    
    this%dt = huge(zero)
    this%ab = one
    this%cf = one
    
    !! Solve for conductive state at degree zero
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%solve_temp_sub( ijmstart=1, ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.false. )
    
    !! Set the time-stepping to the initial choice
    this%dt = dt_help
    this%ab = ab_help
    this%cf = cf_help
    
    !! Add a random perturbation to the conductive state
    do ir = 1, this%nd+1
      is  = 3*(ir-1)+1
      
      do ij = 1, this%jmax
        im = 0
          call random_number( realp )
          
          this%sol%temp(is,jm(ij,im)) = cmplx( realp, zero, kind=dbl ) / 1e4
        
        do im = 1, ij
          call random_number( realp )
          call random_number( imagp )
          
          this%sol%temp(is,jm(ij,im)) = cmplx( realp, imagp, kind=dbl ) / 1e4
        end do
      end do
    end do
    
  end subroutine init_conduction_sub
  
  subroutine init_fromFile_sub(this)
    class(T_ocean),  intent(inout) :: this
    integer                        :: ir, ijm, ndI, jmsI, jmvI
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: velc(:), temp(:,:), spher1(:,:), torr(:,:), spher2(:,:)
    
    ndI  = nd_init_ocean
    jmsI = jm(jmax_init_ocean,jmax_init_ocean)
    jmvI = jml(jmax_init_ocean,jmax_init_ocean,+1)
    
    allocate( r(ndI+1)         )
    allocate( velc(jmvI)       )
    allocate( temp(ndI+1,jmsI) )
    
    allocate( spher1(ndI+1,jmsI) ); spher1 = czero
    allocate( spher2(ndI+1,jmsI) ); spher2 = czero
    allocate( torr(ndI+1,jmsI)   ); torr = czero
    
    open(unit=8, file='code/ocean/inittemp', status='old', action='read')
      do ir = 1, ndI+1
        read(8,*) r(ir), temp(ir,:)
      end do
    close(8)
    
    open(unit=8, file='code/ocean/initvelc', status='old', action='read')
      do ir = 1, ndI+1
        read(8,*) r(ir), velc
        
        do ijm = 2, jmsI
          spher1(ir,ijm) = velc(3*(ijm-1)-1)
          torr(  ir,ijm) = velc(3*(ijm-1)  )
          spher2(ir,ijm) = velc(3*(ijm-1)+1)
        end do
      end do
    close(8)
    
    deallocate(velc)
    
    do ir = 1, this%nd+1
      this%sol%temp(3*(ir-1)+1,:) = this%rad_grid%interpolation_fn( this%jms, ir, r, temp   )
      this%sol%mech(6*(ir-1)+1,:) = this%rad_grid%interpolation_fn( this%jms, ir, r, spher1 )
      this%sol%torr(3*(ir-1)+1,:) = this%rad_grid%interpolation_fn( this%jms, ir, r, torr   )
      this%sol%mech(6*(ir-1)+2,:) = this%rad_grid%interpolation_fn( this%jms, ir, r, spher2 )
    end do
    
    deallocate( r, spher1, spher2, torr, temp )
    
  end subroutine init_fromFile_sub
  
  module procedure init_state_ocean_sub
    real(kind=dbl) :: ab_help, cf_help
    
    if (.not. init_through_file_ocean) then
      call init_conduction_sub(this)
    else
      call init_fromFile_sub(this)
    end if
    
    !! First step with implicit stepping
    ab_help = this%ab
    cf_help = this%cf
    
    this%ab = one
    this%cf = one
    
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
      call this%time_scheme_sub()
      call this%vypis_ocean_sub()
    
    !! Set the time-step to the initial choice
    this%ab = ab_help
    this%cf = cf_help
    
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
  end procedure init_state_ocean_sub
  
end submodule init_state