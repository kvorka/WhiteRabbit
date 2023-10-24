submodule (PhysicalObject) SpheroidalViscel
  implicit none
  
  contains
  
  pure function matica_mech_chb_viscel_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: i
    real(kind=dbl)                      :: j

    allocate( matica(23,6*this%nd+2) ) ; matica = zero
  
    associate( grid => this%rad_grid ); j = i2r_fn(j_in)
    
    select case (this%mechanic_bnd)
      case('shape')
        matica(12, 1) = +sqrt(j    /(2*j+1)) * grid%c(1,-1)
        matica(13, 1) = -sqrt((j+1)/(2*j+1)) * grid%c(1,-1)
        matica(14, 1) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rad / this%gd
        matica(15, 1) = -sqrt(1._dbl/3._dbl)                 / this%Rad / this%gd
        matica(16, 1) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rad / this%gd
        matica(17, 1) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rad / this%gd
        matica(18, 1) = +sqrt(j    /(2*j+1)) * grid%c(1,+1)
        matica(19, 1) = -sqrt((j+1)/(2*j+1)) * grid%c(1,+1)
                
        matica(13, 2) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
        matica(15, 2) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
        matica(16, 2) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
    end select
    
    do i = 1, this%nd
      if (i > 1) then
        matica( 2, 6*(i-1)+1) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(i,-2)                                 )
        matica( 3, 6*(i-1)+1) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(i,-2)                                 )
        matica( 4, 6*(i-1)+1) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(i,-2)                                 )
        matica( 8, 6*(i-1)+1) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(i,-1) - grid%cc(i,-1)*(j-2)/grid%rr(i))
        matica( 9, 6*(i-1)+1) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(i,-1) + grid%cc(i,-1)*(j+1)/grid%rr(i))
        matica(10, 6*(i-1)+1) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(i,-1) + grid%cc(i,-1)*(j+1)/grid%rr(i))
        matica(12, 6*(i-1)+1) = 0._dbl
        matica(14, 6*(i-1)+1) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(i,+1) - grid%cc(i,+1)*(j-2)/grid%rr(i))
        matica(15, 6*(i-1)+1) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(i,+1) + grid%cc(i,+1)*(j+1)/grid%rr(i))
        matica(16, 6*(i-1)+1) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(i,+1) + grid%cc(i,+1)*(j+1)/grid%rr(i))
        matica(20, 6*(i-1)+1) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(i,+2)                                 )
        matica(21, 6*(i-1)+1) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(i,+2)                                 )
        matica(22, 6*(i-1)+1) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(i,+2)                                 )
        
        matica( 2, 6*(i-1)+2) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(i,-2)                                 )
        matica( 3, 6*(i-1)+2) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(i,-2)                                 )
        matica( 4, 6*(i-1)+2) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(i,-2)                                 )
        matica( 8, 6*(i-1)+2) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(i,-1) - grid%cc(i,-1)*(j  )/grid%rr(i))
        matica( 9, 6*(i-1)+2) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(i,-1) - grid%cc(i,-1)*(j  )/grid%rr(i))
        matica(10, 6*(i-1)+2) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(i,-1) + grid%cc(i,-1)*(j+3)/grid%rr(i))
        matica(12, 6*(i-1)+2) = 0._dbl
        matica(14, 6*(i-1)+2) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(i,+1) - grid%cc(i,+1)*(j  )/grid%rr(i))
        matica(15, 6*(i-1)+2) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(i,+1) - grid%cc(i,+1)*(j  )/grid%rr(i))
        matica(16, 6*(i-1)+2) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(i,+1) + grid%cc(i,+1)*(j+3)/grid%rr(i))
        matica(20, 6*(i-1)+2) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(i,+2)                                 )
        matica(21, 6*(i-1)+2) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(i,+2)                                 )
        matica(22, 6*(i-1)+2) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(i,+2)                                 )
      end if
      
      matica( 4, 6*(i-1)+3) = +sqrt(j  )*(grid%d(i,-2)                               )
      matica( 5, 6*(i-1)+3) = -sqrt(j+1)*(grid%d(i,-2)                               )
      matica(10, 6*(i-1)+3) = +sqrt(j  )*(grid%d(i,-1) - grid%c(i,-1)*(j-1)/grid%r(i))
      matica(11, 6*(i-1)+3) = -sqrt(j+1)*(grid%d(i,-1) + grid%c(i,-1)*(j+2)/grid%r(i))
      matica(16, 6*(i-1)+3) = +sqrt(j  )*(grid%d(i,+1) - grid%c(i,+1)*(j-1)/grid%r(i))
      matica(17, 6*(i-1)+3) = -sqrt(j+1)*(grid%d(i,+1) + grid%c(i,+1)*(j+2)/grid%r(i))
      matica(22, 6*(i-1)+3) = +sqrt(j  )*(grid%d(i,+2)                               )
      matica(23, 6*(i-1)+3) = -sqrt(j+1)*(grid%d(i,+2)                               )
      
      matica( 3, 6*(i-1)+4) = -2*sqrt((j-1)/(2*j-1))*(grid%d(i,-2)                           ) * this%Ramu
      matica( 9, 6*(i-1)+4) = -2*sqrt((j-1)/(2*j-1))*(grid%d(i,-1) + grid%c(i,-1)*j/grid%r(i)) * this%Ramu
      matica(11, 6*(i-1)+4) = 1._dbl + this%Ramu * this%dt / this%visc_fn(i)
      matica(15, 6*(i-1)+4) = -2*sqrt((j-1)/(2*j-1))*(grid%d(i,+1) + grid%c(i,+1)*j/grid%r(i)) * this%Ramu
      matica(21, 6*(i-1)+4) = -2*sqrt((j-1)/(2*j-1))*(grid%d(i,+2)                           ) * this%Ramu
      
      matica( 2, 6*(i-1)+5) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(i,-2)                               ) * this%Ramu
      matica( 3, 6*(i-1)+5) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(i,-2)                               ) * this%Ramu
      matica( 8, 6*(i-1)+5) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(i,-1) - grid%c(i,-1)*(j-1)/grid%r(i)) * this%Ramu
      matica( 9, 6*(i-1)+5) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(i,-1) + grid%c(i,-1)*(j+2)/grid%r(i)) * this%Ramu
      matica(12, 6*(i-1)+5) = 1._dbl + this%Ramu * this%dt / this%visc_fn(i)
      matica(14, 6*(i-1)+5) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(i,+1) - grid%c(i,+1)*(j-1)/grid%r(i)) * this%Ramu
      matica(15, 6*(i-1)+5) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(i,+1) + grid%c(i,+1)*(j+2)/grid%r(i)) * this%Ramu
      matica(20, 6*(i-1)+5) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(i,+2)                               ) * this%Ramu
      matica(21, 6*(i-1)+5) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(i,+2)                               ) * this%Ramu
      
      matica( 2, 6*(i-1)+6) = +2*sqrt((j+2)/(2*j+3))*(grid%d(i,-2)                               ) * this%Ramu
      matica( 8, 6*(i-1)+6) = +2*sqrt((j+2)/(2*j+3))*(grid%d(i,-1) - grid%c(i,-1)*(j+1)/grid%r(i)) * this%Ramu
      matica(12, 6*(i-1)+6) = 1._dbl + this%Ramu * this%dt / this%visc_fn(i)
      matica(14, 6*(i-1)+6) = +2*sqrt((j+2)/(2*j+3))*(grid%d(i,+1) - grid%c(i,+1)*(j+1)/grid%r(i)) * this%Ramu
      matica(20, 6*(i-1)+6) = +2*sqrt((j+2)/(2*j+3))*(grid%d(i,+2)                               ) * this%Ramu
    end do
    
    select case (this%mechanic_bnd)
      case('shape')
        matica( 8, 6*this%nd+1) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
        matica(10, 6*this%nd+1) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
        matica(11, 6*this%nd+1) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        
        matica( 5, 6*this%nd+2) = +sqrt(j    /(2*j+1)) * grid%c(this%nd,-1)
        matica( 6, 6*this%nd+2) = -sqrt((j+1)/(2*j+1)) * grid%c(this%nd,-1)
        matica( 7, 6*this%nd+2) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rau / this%gu
        matica( 8, 6*this%nd+2) = -sqrt(1._dbl/3._dbl)                 / this%Rau / this%gu
        matica( 9, 6*this%nd+2) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rau / this%gu
        matica(10, 6*this%nd+2) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rau / this%gu
        matica(11, 6*this%nd+2) = +sqrt(j    /(2*j+1)) * grid%c(this%nd,+1)
        matica(12, 6*this%nd+2) = -sqrt((j+1)/(2*j+1)) * grid%c(this%nd,+1)
    end select  
      
    end associate
  
  end function matica_mech_chb_viscel_fn
  
  pure function matica_mech_hom_viscel_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: i
    real(kind=dbl)                      :: j, dr
    
    allocate( matica(15,6*this%nd+2) ) ; matica = zero
    
    associate( grid => this%rad_grid ); j = i2r_fn(j_in); dr = grid%r(2)-grid%r(1)
    
    select case (this%mechanic_bnd)
      case('shape')
        matica( 8,1) = +sqrt(j    /(2*j+1)) * this%dt / 2
        matica( 9,1) = -sqrt((j+1)/(2*j+1)) * this%dt / 2
        matica(10,1) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rad / this%gd
        matica(11,1) = -sqrt(1._dbl/3._dbl)                 / this%Rad / this%gd
        matica(12,1) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rad / this%gd
        matica(13,1) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rad / this%gd
        matica(14,1) = +sqrt(j    /(2*j+1)) * this%dt / 2
        matica(15,1) = -sqrt((j+1)/(2*j+1)) * this%dt / 2
        
        matica( 9,2) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
        matica(11,2) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
        matica(12,2) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        
      case('frees')
        matica( 8,1) = +sqrt(j  ) / 2
        matica( 9,1) = -sqrt(j+1) / 2
        matica(14,1) = +sqrt(j  ) / 2
        matica(15,1) = -sqrt(j+1) / 2
        
        matica( 9,2) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
        matica(11,2) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
        matica(12,2) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        
      case('noslp')
        matica( 8,1) = 0.5_dbl
        matica(14,1) = 0.5_dbl
        
        matica( 8,2) = 0.5_dbl
        matica(14,2) = 0.5_dbl
        
    end select
      
    do i = 1, this%nd
      if (i > 1) then
        matica( 4, 6*(i-1)+1) = +a_in*sqrt((j-1)        /   (2*j-1)         )*(-1/dr - (j-2)/grid%rr(i)/2)
        matica( 5, 6*(i-1)+1) = -a_in*sqrt((j  )        /(3*(2*j+1)        ))*(-1/dr + (j+1)/grid%rr(i)/2)
        matica( 6, 6*(i-1)+1) = -a_in*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(-1/dr + (j+1)/grid%rr(i)/2)
        matica( 8, 6*(i-1)+1) = -1/(this%Pr * this%dt)
        matica(10, 6*(i-1)+1) = +a_in*sqrt((j-1)        /   (2*j-1)         )*(+1/dr - (j-2)/grid%rr(i)/2)
        matica(11, 6*(i-1)+1) = -a_in*sqrt((j  )        /(3*(2*j+1)        ))*(+1/dr + (j+1)/grid%rr(i)/2)
        matica(12, 6*(i-1)+1) = -a_in*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(+1/dr + (j+1)/grid%rr(i)/2)
        
        matica( 4, 6*(i-1)+2) = +a_in*sqrt((j+1)        /(3*(2*j+1)        ))*(-1/dr - (j  )/grid%rr(i)/2)
        matica( 5, 6*(i-1)+2) = +a_in*sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(-1/dr - (j  )/grid%rr(i)/2)
        matica( 6, 6*(i-1)+2) = -a_in*sqrt((j+2)        /(2*j+3)            )*(-1/dr + (j+3)/grid%rr(i)/2)
        matica( 8, 6*(i-1)+2) = -1/(this%Pr * this%dt)
        matica(10, 6*(i-1)+2) = +a_in*sqrt((j+1)        /(3*(2*j+1)        ))*(+1/dr - (j  )/grid%rr(i)/2)
        matica(11, 6*(i-1)+2) = +a_in*sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(+1/dr - (j  )/grid%rr(i)/2)
        matica(12, 6*(i-1)+2) = -a_in*sqrt((j+2)        /(2*j+3)            )*(+1/dr + (j+3)/grid%rr(i)/2)
      end if
      
      matica( 6, 6*(i-1)+3) = +sqrt(j  )*(-1/dr - (j-1)/grid%r(i)/2)
      matica( 7, 6*(i-1)+3) = -sqrt(j+1)*(-1/dr + (j+2)/grid%r(i)/2)
      matica(12, 6*(i-1)+3) = +sqrt(j  )*(+1/dr - (j-1)/grid%r(i)/2)
      matica(13, 6*(i-1)+3) = -sqrt(j+1)*(+1/dr + (j+2)/grid%r(i)/2)
      
      matica( 5, 6*(i-1)+4) = -2*sqrt((j-1)/(2*j-1))*(-1/dr + j/grid%r(i)/2) * this%visc_fn(i)
      matica( 7, 6*(i-1)+4) = 1._dbl + this%Ramu * this%visc_fn(i) / this%dt
      matica(11, 6*(i-1)+4) = -2*sqrt((j-1)/(2*j-1))*(+1/dr + j/grid%r(i)/2) * this%visc_fn(i)
      
      matica( 4, 6*(i-1)+5) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(-1/dr - (j-1)/grid%r(i)/2) * this%visc_fn(i)
      matica( 5, 6*(i-1)+5) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(-1/dr + (j+2)/grid%r(i)/2) * this%visc_fn(i)
      matica( 8, 6*(i-1)+5) = 1._dbl + this%Ramu * this%visc_fn(i) / this%dt 
      matica(10, 6*(i-1)+5) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(+1/dr - (j-1)/grid%r(i)/2) * this%visc_fn(i)
      matica(11, 6*(i-1)+5) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(+1/dr + (j+2)/grid%r(i)/2) * this%visc_fn(i)
      
      matica( 4, 6*(i-1)+6) = +2*sqrt((j+2)/(2*j+3))*(-1/dr - (j+1)/grid%r(i)/2) * this%visc_fn(i)
      matica( 8, 6*(i-1)+6) = 1._dbl + this%Ramu * this%visc_fn(i) / this%dt
      matica(10, 6*(i-1)+6) = +2*sqrt((j+2)/(2*j+3))*(+1/dr - (j+1)/grid%r(i)/2) * this%visc_fn(i)
    end do
    
    select case (this%mechanic_bnd)
      case('shape')
        matica(4, 6*this%nd+1) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
        matica(6, 6*this%nd+1) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
        matica(7, 6*this%nd+1) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        
        matica(1, 6*this%nd+2) = +sqrt(j    /(2*j+1)) * this%dt / 2
        matica(2, 6*this%nd+2) = -sqrt((j+1)/(2*j+1)) * this%dt / 2
        matica(3, 6*this%nd+2) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rau / this%gu
        matica(4, 6*this%nd+2) = -sqrt(1._dbl/3._dbl)                 / this%Rau / this%gu
        matica(5, 6*this%nd+2) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rau / this%gu
        matica(6, 6*this%nd+2) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rau / this%gu
        matica(7, 6*this%nd+2) = +sqrt(j    /(2*j+1)) * this%dt / 2
        matica(8, 6*this%nd+2) = -sqrt((j+1)/(2*j+1)) * this%dt / 2
        
      case('frees')
        matica(4, 6*this%nd+1) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
        matica(6, 6*this%nd+1) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
        matica(7, 6*this%nd+1) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        
        matica(1, 6*this%nd+2) = +sqrt(j    ) / 2
        matica(2, 6*this%nd+2) = -sqrt((j+1)) / 2
        matica(7, 6*this%nd+2) = +sqrt(j    ) / 2
        matica(8, 6*this%nd+2) = -sqrt((j+1)) / 2
        
      case('noslp')
        matica(2,6*this%nd+1) = 0.5_dbl
        matica(8,6*this%nd+1) = 0.5_dbl
        
        matica(2,6*this%nd+1) = 0.5_dbl
        matica(8,6*this%nd+1) = 0.5_dbl
    end select
    
    end associate
      
  end function matica_mech_hom_viscel_fn
  
end submodule SpheroidalViscel