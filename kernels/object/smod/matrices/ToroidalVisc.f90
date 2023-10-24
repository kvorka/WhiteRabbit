submodule (PhysicalObject) ToroidalVisc
  implicit none
  
  contains
  
  pure function matica_torr_chb_viscos_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: i
    real(kind=dbl)                      :: j

    allocate( matica(11, 3*this%nd+1) ) ; matica = 0._dbl
  
    associate( grid => this%rad_grid ); j = real(j_in, kind=dbl)
    
    select case (this%mechanic_bnd)
      case('frees')
        matica(7,1) = +sqrt((j-1)/(2*(2*j+1)))
        matica(8,1) = -sqrt((j+2)/(2*(2*j+1)))
      
      case('noslp')
        matica(6,1) = grid%c(1,-1)
        matica(9,1) = grid%c(1,+1)
    end select
  
    do i = 1, this%nd
      if (i > 1) then
        matica( 1,3*(i-1)+1) = +a_in*sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,-2)                                 )
        matica( 2,3*(i-1)+1) = -a_in*sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,-2)                                 )
        matica( 4,3*(i-1)+1) = +a_in*sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,-1) - grid%cc(i,-1)*(j-1)/grid%rr(i))
        matica( 5,3*(i-1)+1) = -a_in*sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,-1) + grid%cc(i,-1)*(j+2)/grid%rr(i))
        matica( 6,3*(i-1)+1) = -1/(this%Pr * this%dt)
        matica( 7,3*(i-1)+1) = +a_in*sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,+1) - grid%cc(i,+1)*(j-1)/grid%rr(i))
        matica( 8,3*(i-1)+1) = -a_in*sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,+1) + grid%cc(i,+1)*(j+2)/grid%rr(i))
        matica(10,3*(i-1)+1) = +a_in*sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,+2)                                 )
        matica(11,3*(i-1)+1) = -a_in*sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,+2)                                 )
      end if
      
      matica( 2,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,-2)                               )
      matica( 5,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,-1) + grid%c(i,-1)*(j+1)/grid%r(i))
      matica( 6,3*(i-1)+2) = 1._dbl
      matica( 8,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,+1) + grid%c(i,+1)*(j+1)/grid%r(i))
      matica(11,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,+2)                               )
      
      matica( 1,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,-2)                               )
      matica( 4,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,-1) - grid%c(i,-1)*(j  )/grid%r(i))
      matica( 6,3*(i-1)+3) = 1._dbl
      matica( 7,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,+1) - grid%c(i,+1)*(j  )/grid%r(i))
      matica(10,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,+2)                               )
    end do
    
    select case (this%mechanic_bnd)
      case('frees')
        matica(4,3*this%nd+1) = +sqrt((j-1)/(2*(2*j+1)))
        matica(5,3*this%nd+1) = -sqrt((j+2)/(2*(2*j+1)))
        
      case('noslp')
        matica(3,3*this%nd+1) = grid%c(this%nd,-1)
        matica(6,3*this%nd+1) = grid%c(this%nd,+1)
    end select
    
    end associate
  
  end function matica_torr_chb_viscos_fn

  pure function matica_torr_chb_christ_viscos_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: i
    real(kind=dbl)                      :: j

    allocate( matica(11, 3*this%nd+1) ) ; matica = 0._dbl
  
    associate( grid => this%rad_grid ); j = real(j_in, kind=dbl)
    
    select case (this%mechanic_bnd)
      case('frees')
        matica(7,1) = +sqrt((j-1)/(2*(2*j+1)))
        matica(8,1) = -sqrt((j+2)/(2*(2*j+1)))
      
      case('noslp')
        matica(6,1) = grid%c(1,-1)
        matica(9,1) = grid%c(1,+1)
    end select
  
    do i = 1, this%nd
      if (i > 1) then
        matica( 1,3*(i-1)+1) = +this%Ek * a_in * sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,-2)                                 )
        matica( 2,3*(i-1)+1) = -this%Ek * a_in * sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,-2)                                 )
        matica( 4,3*(i-1)+1) = +this%Ek * a_in * sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,-1) - grid%cc(i,-1)*(j-1)/grid%rr(i))
        matica( 5,3*(i-1)+1) = -this%Ek * a_in * sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,-1) + grid%cc(i,-1)*(j+2)/grid%rr(i))
        matica( 6,3*(i-1)+1) = -1/(this%dt)
        matica( 7,3*(i-1)+1) = +this%Ek * a_in * sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,+1) - grid%cc(i,+1)*(j-1)/grid%rr(i))
        matica( 8,3*(i-1)+1) = -this%Ek * a_in * sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,+1) + grid%cc(i,+1)*(j+2)/grid%rr(i))
        matica(10,3*(i-1)+1) = +this%Ek * a_in * sqrt((j-1)/(2*(2*j+1)))*(grid%dd(i,+2)                                 )
        matica(11,3*(i-1)+1) = -this%Ek * a_in * sqrt((j+2)/(2*(2*j+1)))*(grid%dd(i,+2)                                 )
      end if
      
      matica( 2,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,-2)                               )
      matica( 5,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,-1) + grid%c(i,-1)*(j+1)/grid%r(i))
      matica( 6,3*(i-1)+2) = 1._dbl
      matica( 8,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,+1) + grid%c(i,+1)*(j+1)/grid%r(i))
      matica(11,3*(i-1)+2) = -2*sqrt((j-1)/(2*(2*j+1)))*(grid%d(i,+2)                               )
      
      matica( 1,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,-2)                               )
      matica( 4,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,-1) - grid%c(i,-1)*(j  )/grid%r(i))
      matica( 6,3*(i-1)+3) = 1._dbl
      matica( 7,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,+1) - grid%c(i,+1)*(j  )/grid%r(i))
      matica(10,3*(i-1)+3) = +2*sqrt((j+2)/(2*(2*j+1)))*(grid%d(i,+2)                               )
    end do
    
    select case (this%mechanic_bnd)
      case('frees')
        matica(4,3*this%nd+1) = +sqrt((j-1)/(2*(2*j+1)))
        matica(5,3*this%nd+1) = -sqrt((j+2)/(2*(2*j+1)))
        
      case('noslp')
        matica(3,3*this%nd+1) = grid%c(this%nd,-1)
        matica(6,3*this%nd+1) = grid%c(this%nd,+1)
    end select
    
    end associate
  
  end function matica_torr_chb_christ_viscos_fn
  
end submodule ToroidalVisc