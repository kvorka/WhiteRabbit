submodule (PhysicalObject) ThermalMatrices
  implicit none
  
  contains
  
  pure function matica_temp_hom_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in)  :: this
    integer,                 intent(in)  :: j_in
    real(kind=dbl),          intent(in)  :: a_in
    real(kind=dbl),          allocatable :: matica(:,:)
    integer                              :: i
    real(kind=dbl)                       :: j, q
    
    allocate( matica(7,3*this%nd+1) )
    
    associate(grid => this%rad_grid); j = real(j_in, kind=dbl); matica = 0._dbl
    
    select case (this%thermal_bnd)
      case('phase')
        if ( j_in == 0 ) then
          matica(4, 1) = 0.5_dbl
          matica(7, 1) = 0.5_dbl
        else
          q = -real(-this%sol%flux_fn(1,1,1), kind=dbl) / sqrt(4*pi) / this%lambda_fn(1)
            matica(4, 1) = 0.5_dbl / q
            matica(5, 1) = -sqrt((j  )/(2*j+1)) * this%Raf * this%dt
            matica(6, 1) = +sqrt((j+1)/(2*j+1)) * this%Raf * this%dt
            matica(7, 1) = 0.5_dbl / q      
        end if
      
      case('basic')
        matica(4, 1) = 0.5_dbl
        matica(7, 1) = 0.5_dbl
    end select
    
    do i = 1, this%nd
      if (i > 1) then
        matica(2, 3*(i-1)+1) = +a_in*sqrt((j  )/(2*j+1))*(-1/(grid%r(i)-grid%r(i-1)) - (j-1)/grid%rr(i)/2) / this%cp_fn(i)
        matica(3, 3*(i-1)+1) = -a_in*sqrt((j+1)/(2*j+1))*(-1/(grid%r(i)-grid%r(i-1)) + (j+2)/grid%rr(i)/2) / this%cp_fn(i)
        matica(4, 3*(i-1)+1) = 1/this%dt
        matica(5, 3*(i-1)+1) = +a_in*sqrt((j  )/(2*j+1))*(+1/(grid%r(i)-grid%r(i-1)) - (j-1)/grid%rr(i)/2) / this%cp_fn(i)
        matica(6, 3*(i-1)+1) = -a_in*sqrt((j+1)/(2*j+1))*(+1/(grid%r(i)-grid%r(i-1)) + (j+2)/grid%rr(i)/2) / this%cp_fn(i)
      end if
      
      matica(3, 3*(i-1)+2) = +sqrt((j  )/(2*j+1))*(-1/(grid%rr(i+1)-grid%rr(i)) + (j+1)/grid%r(i)/2) * this%lambda_fn(i)
      matica(4, 3*(i-1)+2) = 1._dbl
      matica(6, 3*(i-1)+2) = +sqrt((j  )/(2*j+1))*(+1/(grid%rr(i+1)-grid%rr(i)) + (j+1)/grid%r(i)/2) * this%lambda_fn(i)
      
      matica(2, 3*(i-1)+3) = -sqrt((j+1)/(2*j+1))*(-1/(grid%rr(i+1)-grid%rr(i)) - (j  )/grid%r(i)/2) * this%lambda_fn(i)
      matica(4, 3*(i-1)+3) = 1._dbl
      matica(5, 3*(i-1)+3) = -sqrt((j+1)/(2*j+1))*(+1/(grid%rr(i+1)-grid%rr(i)) - (j  )/grid%r(i)/2) * this%lambda_fn(i)
    end do
    
    select case (this%thermal_bnd)
      case('phase')
        if ( j_in == 0 ) then
          matica(1, 3*this%nd+1) = 0.5_dbl
          matica(4, 3*this%nd+1) = 0.5_dbl
        else
          q = -real(-this%sol%flux_fn(this%nd,1,1), kind=dbl) / sqrt(4*pi) / this%lambda_fn(this%nd)
            matica(1, 3*this%nd+1) = 0.5_dbl / q
            matica(4, 3*this%nd+1) = 0.5_dbl / q
        end if
      
      case('basic')
        matica(1, 3*this%nd+1) = 0.5_dbl
        matica(4, 3*this%nd+1) = 0.5_dbl
    end select
    
    end associate
    
  end function matica_temp_hom_fn
  
  pure function matica_temp_chb_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: i
    real(kind=dbl)                      :: j, q

    allocate(matica(11,3*this%nd+1) )
    
    associate( grid => this%rad_grid ); j = real(j_in, kind=dbl); matica = 0._dbl
    
    select case (this%thermal_bnd)
      case('phase')
        if ( j_in == 0 ) then
          matica(6,1) = grid%c(1,-1)
          matica(9,1) = grid%c(1,+1)
        else
          q = -real(-this%sol%flux_fn(1,1,1), kind=dbl) / sqrt(4*pi) / this%lambda_fn(1)
            matica(6,1) = grid%c(1,-1) / q
            matica(7,1) = -sqrt((j  )/(2*j+1)) * this%Raf * this%dt
            matica(8,1) = +sqrt((j+1)/(2*j+1)) * this%Raf * this%dt
            matica(9,1) = grid%c(1,+1) / q
        end if
      
      case('basic')
        matica(6,1) = grid%c(1,-1)
        matica(9,1) = grid%c(1,+1)
    end select
    
    do i = 1, this%nd
      if (i > 1) then
        matica( 1,3*(i-1)+1) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(i,-2)                                 ) / this%cp_fn(i)
        matica( 2,3*(i-1)+1) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(i,-2)                                 ) / this%cp_fn(i)
        matica( 4,3*(i-1)+1) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(i,-1) - grid%cc(i,-1)*(j-1)/grid%rr(i)) / this%cp_fn(i)
        matica( 5,3*(i-1)+1) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(i,-1) + grid%cc(i,-1)*(j+2)/grid%rr(i)) / this%cp_fn(i)
        matica( 6,3*(i-1)+1) = 1/this%dt
        matica( 7,3*(i-1)+1) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(i,+1) - grid%cc(i,+1)*(j-1)/grid%rr(i)) / this%cp_fn(i)
        matica( 8,3*(i-1)+1) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(i,+1) + grid%cc(i,+1)*(j+2)/grid%rr(i)) / this%cp_fn(i)
        matica(10,3*(i-1)+1) = +a_in*sqrt((j  )/(2*j+1))*(grid%dd(i,+2)                                 ) / this%cp_fn(i)
        matica(11,3*(i-1)+1) = -a_in*sqrt((j+1)/(2*j+1))*(grid%dd(i,+2)                                 ) / this%cp_fn(i)
      end if
    
      matica( 2,3*(i-1)+2) = +sqrt((j  )/(2*j+1))*(grid%d(i,-2)                               ) * this%lambda_fn(i)
      matica( 5,3*(i-1)+2) = +sqrt((j  )/(2*j+1))*(grid%d(i,-1) + grid%c(i,-1)*(j+1)/grid%r(i)) * this%lambda_fn(i)
      matica( 6,3*(i-1)+2) = 1._dbl
      matica( 8,3*(i-1)+2) = +sqrt((j  )/(2*j+1))*(grid%d(i,+1) + grid%c(i,+1)*(j+1)/grid%r(i)) * this%lambda_fn(i)
      matica(11,3*(i-1)+2) = +sqrt((j  )/(2*j+1))*(grid%d(i,+2)                               ) * this%lambda_fn(i)
      
      matica( 1,3*(i-1)+3) = -sqrt((j+1)/(2*j+1))*(grid%d(i,-2)                               ) * this%lambda_fn(i)
      matica( 4,3*(i-1)+3) = -sqrt((j+1)/(2*j+1))*(grid%d(i,-1) - grid%c(i,-1)*(j  )/grid%r(i)) * this%lambda_fn(i)
      matica( 6,3*(i-1)+3) = 1._dbl
      matica( 7,3*(i-1)+3) = -sqrt((j+1)/(2*j+1))*(grid%d(i,+1) - grid%c(i,+1)*(j  )/grid%r(i)) * this%lambda_fn(i)
      matica(10,3*(i-1)+3) = -sqrt((j+1)/(2*j+1))*(grid%d(i,+2)                               ) * this%lambda_fn(i)
    end do
    
    select case (this%thermal_bnd)
      case('phase')
        if ( j_in == 0 ) then
          matica(3,3*this%nd+1) = grid%c(this%nd,-1)
          matica(6,3*this%nd+1) = grid%c(this%nd,+1)
        else
          q = -real(-this%sol%flux_fn(this%nd,1,1), kind=dbl) / sqrt(4*pi) / this%lambda_fn(this%nd)
            matica(3,3*this%nd+1) = grid%c(this%nd,-1) / q
            matica(6,3*this%nd+1) = grid%c(this%nd,+1) / q    
        end if
      
      case('basic')
        matica(3,3*this%nd+1) = grid%c(this%nd,-1)
        matica(6,3*this%nd+1) = grid%c(this%nd,+1)
    end select
    
    end associate
    
  end function matica_temp_chb_fn
  
end submodule ThermalMatrices