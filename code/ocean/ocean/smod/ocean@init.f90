submodule (ocean) init
  implicit none; contains
  
  module procedure init_ocean_sub
    
    call this%init_objects_sub( nd = nd_ocean, jmax = jmax_ocean, r_ud = r_ud_ocean, rgrid = grid_type_ocean, &
                              & gmod = gravity_ocean , g = (1-r_ud_ocean)**2 , noharm = noharm_ocean          )
    
    call this%bnd%init_flux_up_sub()
    
    this%n_iter = n_iter_ocean
    this%cf     = 0.6_dbl
    this%ab     = 1.5_dbl
    
    this%Pr = Pr_ocean
    this%Ra = Ra_ocean
    this%Ek = Ek_ocean
    this%Cl = Kl_ocean
    
    this%D_ud         = D_ud_ocean
    this%rheology     = rheology_ocean
    this%mechanic_bnd = mechanic_bnd_ocean
    this%thermal_bnd  = thermal_bnd_ocean
    this%scaling      = scaling_ocean
    
    open(unit=11, file='data/Nuss.dat', status='new', action='write')
    open(unit=12, file='data/Laws.dat', status='new', action='write')
    
  end procedure init_ocean_sub
  
  module procedure deallocate_ocean_sub
    
    if ( allocated(this%nsph1) ) deallocate( this%nsph1 )
    if ( allocated(this%nsph2) ) deallocate( this%nsph2 )
    if ( allocated(this%ntorr) ) deallocate( this%ntorr )
    if ( allocated(this%ntemp) ) deallocate( this%ntemp )
    
    close(11); close(12)
    call this%deallocate_objects_sub()

  end procedure deallocate_ocean_sub
  
  module procedure init_state_ocean_sub
    integer                        :: i, j, m, ijm, ndI1, jmsI, jmvI
    real(kind=dbl)                 :: ab_help, dt_help, cf_help, re, im
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: velc(:), temp(:,:), spher1(:,:), torr(:,:), spher2(:,:)

    if (.not. init_through_file_ocean) then
      !! Solve for conductive state at degree zero
      dt_help = this%dt
      ab_help = this%ab
      cf_help = this%cf
      
      this%dt = huge(zero)
      this%ab = one
      this%cf = one
      
      call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
      call this%solve_temp_sub( ijmstart=1, ijmend=this%jms, ijmstep=1, rematrix=.false., matxsol=.false. )
      
      this%dt = dt_help
      this%ab = ab_help
      this%cf = cf_help
      
      !! Add a random perturbation to the conductive state
      do i = 1, this%nd+1
        do j = 1, this%jmax
          do m = 0, j
            ijm = jm(j,m)
            
            if (m == 0) then
              call random_number( re )
              this%sol%temp(3*(i-1)+1, ijm)%re = re / 1e4
            else
              call random_number( re ); call random_number( im )
              this%sol%temp(3*(i-1)+1, ijm) = cmplx(re, im, kind=dbl) / 1e4
            end if
            
          end do
        end do
      end do
      
    else
      ndI1 = nd_init_ocean+1; jmsI = jm(jmax_init_ocean,jmax_init_ocean); jmvI = jml(jmax_init_ocean,jmax_init_ocean,+1)

      allocate( r(ndI1), velc(jmvI), temp(ndI1,jmsI), spher1(ndI1,jmsI), spher2(ndI1,jmsI), torr(ndI1,jmsI) )
        spher1 = czero; spher2 = czero; torr = czero

        open(unit=8, file='code/ocean/inittemp', status='old', action='read')
          do i = 1, ndI1
            read(8,*) r(i), temp(i,:)
          end do
        close(8)

        open(unit=8, file='code/ocean/initvelc', status='old', action='read')
          do i = 1, ndI1
            read(8,*) r(i), velc

            do ijm = 2, jmsI
              spher1(i,ijm) = velc(3*(ijm-1)-1)
              torr(  i,ijm) = velc(3*(ijm-1)  )
              spher2(i,ijm) = velc(3*(ijm-1)+1)
            end do
          end do
        close(8)

      deallocate(velc)

      do i = 1, this%nd+1
        this%sol%temp(3*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, temp  )
        this%sol%mech(6*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, spher1)
        this%sol%torr(3*(i-1)+1,:) = this%rad_grid%interpolation_fn(this%jms, i, r, torr  )
        this%sol%mech(6*(i-1)+2,:) = this%rad_grid%interpolation_fn(this%jms, i, r, spher2)
      end do

      deallocate(r, spher1, spher2, torr, temp)
    end if
    
    ab_help = this%ab
    cf_help = this%cf
    
    this%ab = one
    this%cf = one
    
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
      call this%time_scheme_sub()
      call this%vypis_ocean_sub()
    
    this%ab = ab_help
    this%cf = cf_help
    
    call this%prepare_mat_temp_sub( ijstart=0 , ijend=this%jmax )
    call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
  end procedure init_state_ocean_sub
  
end submodule init