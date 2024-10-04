submodule (physicalobject) thermal_matrix
  implicit none ; contains
  
  module pure function matica_temp_hom_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: j_in
    real(kind=dbl),          intent(in)  :: a_in
    real(kind=dbl),          allocatable :: matica(:,:)
    integer                              :: ir, is
    real(kind=dbl)                       :: j, dT_dr
    
    allocate( matica(7,3*this%nd+1) ); associate(grid => this%rad_grid)
    
    call zero_rarray_sub( 7*(3*this%nd+1), matica )
    
    j = i2r_fn(j_in)
    
    ir = 1
      is = 1
        select case (this%thermal_bnd)
          case('phase')
            if ( j_in == 0 ) then
              matica(4,is) = grid%c(ir,-1)
              matica(7,is) = grid%c(ir,+1)
            else
              dT_dr = c2r_fn( this%dT_dr_r_fn(ir,1) ) / s4pi
                matica(4,is) = grid%c(ir,-1) / ( dT_dr - this%Cl )
                matica(5,is) = -sqrt((j  )/(2*j+1)) * this%Raf * this%dt
                matica(6,is) = +sqrt((j+1)/(2*j+1)) * this%Raf * this%dt
                matica(7,is) = grid%c(ir,+1) / ( dT_dr - this%Cl )
            end if
          
          case('basic')
            matica(4,is) = grid%c(ir,-1)
            matica(7,is) = grid%c(ir,+1)
        end select
    
    do ir = 1, this%nd
      is = 3*(ir-1)+1
        
        if (ir > 1) then
          matica(2,is) = +a_in*sqrt((j  )/(2*j+1))*(grid%hdd(ir,-1) - grid%cc(ir,-1)*(j-1)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica(3,is) = -a_in*sqrt((j+1)/(2*j+1))*(grid%hdd(ir,-1) + grid%cc(ir,-1)*(j+2)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica(4,is) = 1/this%dt
          matica(5,is) = +a_in*sqrt((j  )/(2*j+1))*(grid%hdd(ir,+1) - grid%cc(ir,+1)*(j-1)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica(6,is) = -a_in*sqrt((j+1)/(2*j+1))*(grid%hdd(ir,+1) + grid%cc(ir,+1)*(j+2)/grid%rr(ir)) / this%cp_rr_fn(ir)
        end if
        
        matica(3,is+1) = +sqrt((j  )/(2*j+1))*(grid%hd(ir,-1) + grid%c(ir,-1)*(j+1)/grid%r(ir))
        matica(4,is+1) = 1 / this%lambda_r_fn(ir)
        matica(6,is+1) = +sqrt((j  )/(2*j+1))*(grid%hd(ir,+1) + grid%c(ir,+1)*(j+1)/grid%r(ir))
        
        matica(2,is+2) = -sqrt((j+1)/(2*j+1))*(grid%hd(ir,-1) - grid%c(ir,-1)*(j  )/grid%r(ir))
        matica(4,is+2) = 1 / this%lambda_r_fn(ir)
        matica(5,is+2) = -sqrt((j+1)/(2*j+1))*(grid%hd(ir,+1) - grid%c(ir,+1)*(j  )/grid%r(ir))
    end do
    
    ir = this%nd
      is = 3*this%nd+1
        select case (this%thermal_bnd)
          case('phase')
            if ( j_in == 0 ) then
              matica(1,is) = grid%c(ir,-1)
              matica(4,is) = grid%c(ir,+1)
            else
              dT_dr = c2r_fn( this%dT_dr_r_fn(ir,1) ) / s4pi
                matica(1,is) = grid%c(ir,-1) / dT_dr
                matica(4,is) = grid%c(ir,+1) / dT_dr
            end if
          
          case('basic')
            matica(1,is) = grid%c(ir,-1)
            matica(4,is) = grid%c(ir,+1)
        end select
        
    end associate
    
  end function matica_temp_hom_fn
  
  module pure function matica_temp_chb_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: ir, is
    real(kind=dbl)                      :: j, dT_dr

    allocate(matica(11,3*this%nd+1) ); associate( grid => this%rad_grid )
    
    call zero_rarray_sub( 11*(3*this%nd+1), matica )
    
    j = i2r_fn(j_in)
    
    ir = 1
      is = 1
        select case (this%thermal_bnd)
          case('phase')
            if ( j_in == 0 ) then
              matica(6,is) = grid%c(ir,-1)
              matica(9,is) = grid%c(ir,+1)
            else
              dT_dr = c2r_fn( this%dT_dr_r_fn(ir,1) ) / s4pi
                matica(6,is) = grid%c(ir,-1) / dT_dr
                matica(7,is) = -sqrt((j  )/(2*j+1)) * this%Raf * this%dt
                matica(8,is) = +sqrt((j+1)/(2*j+1)) * this%Raf * this%dt
                matica(9,is) = grid%c(ir,+1) / dT_dr
            end if
          
          case('basic')
            matica(6,is) = grid%c(ir,-1)
            matica(9,is) = grid%c(ir,+1)
        end select
    
    do ir = 1, this%nd
      is = 3*(ir-1)+1
      
        if (ir > 1) then
          matica( 1,is) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(ir,-2)                                   ) / this%cp_rr_fn(ir)
          matica( 2,is) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(ir,-2)                                   ) / this%cp_rr_fn(ir)
          matica( 4,is) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(ir,-1) - grid%cc(ir,-1)*(j-1)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica( 5,is) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(ir,-1) + grid%cc(ir,-1)*(j+2)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica( 6,is) = 1 / this%dt
          matica( 7,is) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(ir,+1) - grid%cc(ir,+1)*(j-1)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica( 8,is) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(ir,+1) + grid%cc(ir,+1)*(j+2)/grid%rr(ir)) / this%cp_rr_fn(ir)
          matica(10,is) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(ir,+2)                                   ) / this%cp_rr_fn(ir)
          matica(11,is) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(ir,+2)                                   ) / this%cp_rr_fn(ir)
        end if
        
        matica( 2,is+1) = +sqrt((j  )/(2*j+1))*(grid%d(ir,-2)                                 )
        matica( 5,is+1) = +sqrt((j  )/(2*j+1))*(grid%d(ir,-1) + grid%c(ir,-1)*(j+1)/grid%r(ir))
        matica( 6,is+1) = 1 / this%lambda_r_fn(ir)
        matica( 8,is+1) = +sqrt((j  )/(2*j+1))*(grid%d(ir,+1) + grid%c(ir,+1)*(j+1)/grid%r(ir))
        matica(11,is+1) = +sqrt((j  )/(2*j+1))*(grid%d(ir,+2)                                 )
        
        matica( 1,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,-2)                                 )
        matica( 4,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,-1) - grid%c(ir,-1)*(j  )/grid%r(ir))
        matica( 6,is+2) = 1 / this%lambda_r_fn(ir)
        matica( 7,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,+1) - grid%c(ir,+1)*(j  )/grid%r(ir))
        matica(10,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,+2)                                 )
    end do
    
    ir = this%nd
      is = 3*ir+1
        select case (this%thermal_bnd)
          case('phase')
            if ( j_in == 0 ) then
              matica(3,is) = grid%c(ir,-1)
              matica(6,is) = grid%c(ir,+1)
            else
              dT_dr = c2r_fn( this%dT_dr_r_fn(ir,1) ) / s4pi
                matica(3,is) = grid%c(ir,-1) / dT_dr
                matica(6,is) = grid%c(ir,+1) / dT_dr
            end if
          
          case('basic')
            matica(3,is) = grid%c(ir,-1)
            matica(6,is) = grid%c(ir,+1)
        end select
    
    end associate
    
  end function matica_temp_chb_fn
  
end submodule thermal_matrix