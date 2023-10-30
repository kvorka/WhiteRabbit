submodule(OutputOceanMod) Velc
  implicit none
  
  contains
  
  subroutine harm_analysis_rad_velc_sub()
    integer,              parameter :: jmax = 100
    integer                         :: ir, ij
    real(kind=dbl),    allocatable  :: velc_r_j0(:), velc_r_grid(:,:)
    complex(kind=dbl), allocatable  :: velc120(:,:)
    
    allocate( velc120(jmv,n_out), velc_r_j0(0:jmax), velc_r_grid(nth,n_out) )
    
      call load_data_sub(jmv, 'velc-averaged.spec', velc120) ; velc_r_j0 = zero
      
      do ir = 1, n_out
        do ij = 2, jmax, 2
          velc_r_j0(ij) = sqrt( (ij  ) / (2*ij+1._dbl) ) * velc120(jml(ij,0,-1),ir) - &
                        & sqrt( (ij+1) / (2*ij+1._dbl) ) * velc120(jml(ij,0,+1),ir)
        end do
        
        call harmsy_Pj0_sub(jmax, velc_r_j0, velc_r_grid(:,ir))
      end do
    
    deallocate( velc120, velc_r_j0 )
    
      velc_r_grid = Ek_ocean / Pr_ocean * velc_r_grid
      
      open(unit=8, file='ocean_rvelc_grid', status='new', action='write')
        write(8,*) 'max: ' , maxval( velc_r_grid )
        write(8,*) 'min: ' , minval( velc_r_grid )
      close(8)
      
      call save_data_sub('inp-r.dat', velc_r_grid)
      
    deallocate( velc_r_grid )
    
  end subroutine harm_analysis_rad_velc_sub
  
  subroutine harm_analysis_zon_velc_sub()
    integer,           parameter   :: jmax = 100
    integer                        :: ir, ij
    real(kind=dbl),    allocatable :: velc_phi_j0(:), velc_phi_grid(:,:)
    complex(kind=dbl), allocatable :: velc120(:,:)
    
    allocate( velc120(jmv,n_out), velc_phi_j0(jmax), velc_phi_grid(nth,n_out) )
    
      call load_data_sub(jmv, 'velc-averaged.spec', velc120) ; velc_phi_j0 = zero
      
      do ir = 1, n_out
        do ij = 1, jmax, 2
          velc_phi_j0(ij) = c2r_fn( -cunit * velc120( jml(ij,0,0) , ir ) )
        end do
        
        call harmsy_Pj1_sub(jmax, velc_phi_j0, velc_phi_grid(:,ir))
      end do
    
    deallocate( velc120, velc_phi_j0 )
    
      velc_phi_grid = Ek_ocean / Pr_ocean * velc_phi_grid
      
      open(unit=8, file='ocean_velc_phi_extrms', status='new', action='write')
        write(8,*) 'max: ' , maxval( velc_phi_grid )
        write(8,*) 'min: ' , minval( velc_phi_grid )
      close(8)
      
      call save_data_sub('inp-p.dat', velc_phi_grid)
      
    deallocate( velc_phi_grid )
    
  end subroutine harm_analysis_zon_velc_sub
  
  subroutine save_spectra_velc_sub()
    real(kind=dbl),    allocatable :: r(:)
    complex(kind=dbl), allocatable :: velc(:,:)
    
    allocate( r(nd_ocean+1), velc(jmv,nd_ocean+1) ) ; velc = czero
    
    call avrg_spectra_3d_sub(path_ocean_velc, r, velc)
    call out_spectra_3d_sub('velc-averaged.spec', r, velc)
    
    deallocate( r, velc )
    
  end subroutine save_spectra_velc_sub
  
end submodule Velc
