submodule (PhysicalObject) SpheroidalViscel
  implicit none ; contains
  
  module pure function matica_mech_hom_viscel_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: ir, is
    real(kind=dbl)                      :: j
    
    allocate( matica(15,6*this%nd+2) ); associate( grid => this%rad_grid )
    
    call zero_rarray_sub( 15*(6*this%nd+2), matica )
    
    j = i2r_fn(j_in)
    
    ir = 1
      is = 1
        select case (this%mechanic_bnd)
          case('shape')
            matica( 8,is) = +sqrt(j    /(2*j+1)) * this%dt * grid%c(ir,-1)
            matica( 9,is) = -sqrt((j+1)/(2*j+1)) * this%dt * grid%c(ir,-1)
            matica(10,is) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rad / this%gd
            matica(11,is) = -sqrt(1._dbl/3._dbl)                 / this%Rad / this%gd
            matica(12,is) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rad / this%gd
            matica(13,is) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rad / this%gd
            matica(14,is) = +sqrt(j    /(2*j+1)) * this%dt * grid%c(ir,+1)
            matica(15,is) = -sqrt((j+1)/(2*j+1)) * this%dt * grid%c(ir,+1)
            
            matica( 9,is+1) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
            matica(11,is+1) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
            matica(12,is+1) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        end select
      
    do ir = 1, this%nd
      is = 6*(ir-1)+1
      
        if (ir > 1) then
          matica( 4,is) = +a_in*sqrt((j-1)        /   (2*j-1)         )*(grid%hdd(ir,-1) - grid%cc(ir,-1)*(j-2)/grid%rr(ir))
          matica( 5,is) = -a_in*sqrt((j  )        /(3*(2*j+1)        ))*(grid%hdd(ir,-1) + grid%cc(ir,-1)*(j+1)/grid%rr(ir))
          matica( 6,is) = -a_in*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%hdd(ir,-1) + grid%cc(ir,-1)*(j+1)/grid%rr(ir))
          matica( 8,is) = -1/(this%Pr * this%dt)
          matica(10,is) = +a_in*sqrt((j-1)        /   (2*j-1)         )*(grid%hdd(ir,+1) - grid%cc(ir,+1)*(j-2)/grid%rr(ir))
          matica(11,is) = -a_in*sqrt((j  )        /(3*(2*j+1)        ))*(grid%hdd(ir,+1) + grid%cc(ir,+1)*(j+1)/grid%rr(ir))
          matica(12,is) = -a_in*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%hdd(ir,+1) + grid%cc(ir,+1)*(j+1)/grid%rr(ir))
          
          matica( 4,is+1) = +a_in*sqrt((j+1)        /(3*(2*j+1)        ))*(grid%hdd(ir,-1) - grid%cc(ir,-1)*(j  )/grid%rr(ir))
          matica( 5,is+1) = +a_in*sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%hdd(ir,-1) - grid%cc(ir,-1)*(j  )/grid%rr(ir))
          matica( 6,is+1) = -a_in*sqrt((j+2)        /(2*j+3)            )*(grid%hdd(ir,-1) + grid%cc(ir,-1)*(j+3)/grid%rr(ir))
          matica( 8,is+1) = -1/(this%Pr * this%dt)
          matica(10,is+1) = +a_in*sqrt((j+1)        /(3*(2*j+1)        ))*(grid%hdd(ir,+1) - grid%cc(ir,+1)*(j  )/grid%rr(ir))
          matica(11,is+1) = +a_in*sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%hdd(ir,+1) - grid%cc(ir,+1)*(j  )/grid%rr(ir))
          matica(12,is+1) = -a_in*sqrt((j+2)        /(2*j+3)            )*(grid%hdd(ir,+1) + grid%cc(ir,+1)*(j+3)/grid%rr(ir))
        end if
        
        matica( 6,is+2) = +sqrt((j  )/(2*j+1))*(grid%hd(ir,-1) - grid%c(ir,-1)*(j-1)/grid%r(ir))
        matica( 7,is+2) = -sqrt((j+1)/(2*j+1))*(grid%hd(ir,-1) + grid%c(ir,-1)*(j+2)/grid%r(ir))
        matica(12,is+2) = +sqrt((j  )/(2*j+1))*(grid%hd(ir,+1) - grid%c(ir,+1)*(j-1)/grid%r(ir))
        matica(13,is+2) = -sqrt((j+1)/(2*j+1))*(grid%hd(ir,+1) + grid%c(ir,+1)*(j+2)/grid%r(ir))
        
        matica( 5,is+3) = -2*sqrt((j-1)/(2*j-1))*(grid%hd(ir,-1) + grid%c(ir,-1)*j/grid%r(ir))
        matica( 7,is+3) = 1 / this%visc_r_fn(ir) + this%Ramu / this%dt
        matica(11,is+3) = -2*sqrt((j-1)/(2*j-1))*(grid%hd(ir,+1) + grid%c(ir,+1)*j/grid%r(ir))
        
        matica( 4,is+4) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%hd(ir,-1) - grid%c(ir,-1)*(j-1)/grid%r(ir))
        matica( 5,is+4) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%hd(ir,-1) + grid%c(ir,-1)*(j+2)/grid%r(ir))
        matica( 8,is+4) = 1 / this%visc_r_fn(ir) + this%Ramu / this%dt 
        matica(10,is+4) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%hd(ir,+1) - grid%c(ir,-1)*(j-1)/grid%r(ir))
        matica(11,is+4) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%hd(ir,+1) + grid%c(ir,-1)*(j+2)/grid%r(ir))
        
        matica( 4,is+5) = +2*sqrt((j+2)/(2*j+3))*(grid%hd(ir,-1) - grid%c(ir,-1)*(j+1)/grid%r(ir))
        matica( 8,is+5) = one / this%visc_r_fn(ir) + this%Ramu / this%dt
        matica(10,is+5) = +2*sqrt((j+2)/(2*j+3))*(grid%hd(ir,+1) - grid%c(ir,+1)*(j+1)/grid%r(ir))
    end do
    
    ir = this%nd
      is = 6*ir+1
        select case (this%mechanic_bnd)
          case('shape')
            matica(4,is) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
            matica(6,is) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
            matica(7,is) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
            
            matica(1,is+1) = +sqrt(j    /(2*j+1)) * this%dt * grid%c(ir,-1)
            matica(2,is+1) = -sqrt((j+1)/(2*j+1)) * this%dt * grid%c(ir,-1)
            matica(3,is+1) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rau / this%gu
            matica(4,is+1) = -sqrt(1._dbl/3._dbl)                 / this%Rau / this%gu
            matica(5,is+1) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rau / this%gu
            matica(6,is+1) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rau / this%gu
            matica(7,is+1) = +sqrt(j    /(2*j+1)) * this%dt * grid%c(ir,+1)
            matica(8,is+1) = -sqrt((j+1)/(2*j+1)) * this%dt * grid%c(ir,+1)
        end select
    
    end associate
    
  end function matica_mech_hom_viscel_fn
  
  module pure function matica_mech_chb_viscel_fn(this, j_in, a_in) result(matica)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: j_in
    real(kind=dbl),          intent(in) :: a_in
    real(kind=dbl),         allocatable :: matica(:,:)
    integer                             :: ir, is
    real(kind=dbl)                      :: j
    
    allocate( matica(23,6*this%nd+2) ) ; associate( grid => this%rad_grid )
    
    call zero_rarray_sub( 23*(6*this%nd+2), matica )
    
    j = i2r_fn(j_in)
    
    ir = 1
      is = 1
        select case (this%mechanic_bnd)
          case('shape')
            matica(12,is) = +sqrt(j    /(2*j+1)) * grid%c(ir,-1)
            matica(13,is) = -sqrt((j+1)/(2*j+1)) * grid%c(ir,-1)
            matica(14,is) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rad / this%gd
            matica(15,is) = -sqrt(1._dbl/3._dbl)                 / this%Rad / this%gd
            matica(16,is) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rad / this%gd
            matica(17,is) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rad / this%gd
            matica(18,is) = +sqrt(j    /(2*j+1)) * grid%c(ir,+1)
            matica(19,is) = -sqrt((j+1)/(2*j+1)) * grid%c(ir,+1)
            
            matica(13,is+1) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
            matica(15,is+1) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
            matica(16,is+1) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
        end select
    
    do ir = 1, this%nd
      is = 6*(ir-1)+1
      
        if (ir > 1) then
          matica( 2,is) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(ir,-2)                                   )
          matica( 3,is) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(ir,-2)                                   )
          matica( 4,is) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(ir,-2)                                   )
          matica( 8,is) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(ir,-1) - grid%cc(ir,-1)*(j-2)/grid%rr(ir))
          matica( 9,is) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(ir,-1) + grid%cc(ir,-1)*(j+1)/grid%rr(ir))
          matica(10,is) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(ir,-1) + grid%cc(ir,-1)*(j+1)/grid%rr(ir))
          matica(12,is) = -1/ ( this%Pr * this%dt )
          matica(14,is) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(ir,+1) - grid%cc(ir,+1)*(j-2)/grid%rr(ir))
          matica(15,is) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(ir,+1) + grid%cc(ir,+1)*(j+1)/grid%rr(ir))
          matica(16,is) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(ir,+1) + grid%cc(ir,+1)*(j+1)/grid%rr(ir))
          matica(20,is) = +sqrt((j-1)        /   (2*j-1)         )*(grid%dd(ir,+2)                                   )
          matica(21,is) = -sqrt((j  )        /(3*(2*j+1)        ))*(grid%dd(ir,+2)                                   )
          matica(22,is) = -sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%dd(ir,+2)                                   )
          
          matica( 2,is+1) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(ir,-2)                                   )
          matica( 3,is+1) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(ir,-2)                                   )
          matica( 4,is+1) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(ir,-2)                                   )
          matica( 8,is+1) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(ir,-1) - grid%cc(ir,-1)*(j  )/grid%rr(ir))
          matica( 9,is+1) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(ir,-1) - grid%cc(ir,-1)*(j  )/grid%rr(ir))
          matica(10,is+1) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(ir,-1) + grid%cc(ir,-1)*(j+3)/grid%rr(ir))
          matica(12,is+1) = -1 / ( this%Pr * this%dt )
          matica(14,is+1) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(ir,+1) - grid%cc(ir,+1)*(j  )/grid%rr(ir))
          matica(15,is+1) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(ir,+1) - grid%cc(ir,+1)*(j  )/grid%rr(ir))
          matica(16,is+1) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(ir,+1) + grid%cc(ir,+1)*(j+3)/grid%rr(ir))
          matica(20,is+1) = +sqrt((j+1)        /(3*(2*j+1)        ))*(grid%dd(ir,+2)                                   )
          matica(21,is+1) = +sqrt((j  )*(2*j-1)/(6*(2*j+3)*(2*j+1)))*(grid%dd(ir,+2)                                   )
          matica(22,is+1) = -sqrt((j+2)        /(2*j+3)            )*(grid%dd(ir,+2)                                   )
        end if
        
        matica( 4,is+2) = +sqrt((j  )/(2*j+1))*(grid%d(ir,-2)                                 )
        matica( 5,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,-2)                                 )
        matica(10,is+2) = +sqrt((j  )/(2*j+1))*(grid%d(ir,-1) - grid%c(ir,-1)*(j-1)/grid%r(ir))
        matica(11,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,-1) + grid%c(ir,-1)*(j+2)/grid%r(ir))
        matica(16,is+2) = +sqrt((j  )/(2*j+1))*(grid%d(ir,+1) - grid%c(ir,+1)*(j-1)/grid%r(ir))
        matica(17,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,+1) + grid%c(ir,+1)*(j+2)/grid%r(ir))
        matica(22,is+2) = +sqrt((j  )/(2*j+1))*(grid%d(ir,+2)                                 )
        matica(23,is+2) = -sqrt((j+1)/(2*j+1))*(grid%d(ir,+2)                                 )
        
        matica( 3,is+3) = -2*sqrt((j-1)/(2*j-1))*(grid%d(ir,-2)                             )
        matica( 9,is+3) = -2*sqrt((j-1)/(2*j-1))*(grid%d(ir,-1) + grid%c(ir,-1)*j/grid%r(ir))
        matica(11,is+3) = 1 / this%visc_r_fn(ir) + this%Ramu / this%dt
        matica(15,is+3) = -2*sqrt((j-1)/(2*j-1))*(grid%d(ir,+1) + grid%c(ir,+1)*j/grid%r(ir))
        matica(21,is+3) = -2*sqrt((j-1)/(2*j-1))*(grid%d(ir,+2)                             )
        
        matica( 2,is+4) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(ir,-2)                                 )
        matica( 3,is+4) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(ir,-2)                                 )
        matica( 8,is+4) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(ir,-1) - grid%c(ir,-1)*(j-1)/grid%r(ir))
        matica( 9,is+4) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(ir,-1) + grid%c(ir,-1)*(j+2)/grid%r(ir))
        matica(12,is+4) = 1 / this%visc_r_fn(ir) + this%Ramu / this%dt
        matica(14,is+4) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(ir,+1) - grid%c(ir,+1)*(j-1)/grid%r(ir))
        matica(15,is+4) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(ir,+1) + grid%c(ir,+1)*(j+2)/grid%r(ir))
        matica(20,is+4) = +2*sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))*(grid%d(ir,+2)                                 )
        matica(21,is+4) = -2*sqrt((j  )*(2*j-1)/(6*(2*j+1)*(2*j+3)))*(grid%d(ir,+2)                                 )
        
        matica( 2,is+5) = +2*sqrt((j+2)/(2*j+3))*(grid%d(ir,-2)                                 )
        matica( 8,is+5) = +2*sqrt((j+2)/(2*j+3))*(grid%d(ir,-1) - grid%c(ir,-1)*(j+1)/grid%r(ir))
        matica(12,is+5) = 1 / this%visc_r_fn(ir) + this%Ramu / this%dt
        matica(14,is+5) = +2*sqrt((j+2)/(2*j+3))*(grid%d(ir,+1) - grid%c(ir,+1)*(j+1)/grid%r(ir))
        matica(20,is+5) = +2*sqrt((j+2)/(2*j+3))*(grid%d(ir,+2)                                 )
    end do
    
    ir = this%nd
      is = 6*this%nd+1
        select case (this%mechanic_bnd)
          case('shape')
            matica( 8,is) = +sqrt((j+1)*(j-1)/(  (2*j+1)*(2*j-1)))
            matica(10,is) = -sqrt((3  )      /(2*(2*j-1)*(2*j+3)))
            matica(11,is) = -sqrt((j  )*(j+2)/(  (2*j+1)*(2*j+3)))
            
            matica( 5,is+1) = +sqrt(j    /(2*j+1)) * grid%c(ir,-1)
            matica( 6,is+1) = -sqrt((j+1)/(2*j+1)) * grid%c(ir,-1)
            matica( 7,is+1) = +sqrt(j*(j-1)/((2*j-1)*(2*j+1)))     / this%Rau / this%gu
            matica( 8,is+1) = -sqrt(1._dbl/3._dbl)                 / this%Rau / this%gu
            matica( 9,is+1) = -sqrt(2*j*(j+1)/(3*(2*j-1)*(2*j+3))) / this%Rau / this%gu
            matica(10,is+1) = +sqrt((j+1)*(j+2)/((2*j+1)*(2*j+3))) / this%Rau / this%gu
            matica(11,is+1) = +sqrt(j    /(2*j+1)) * grid%c(ir,+1)
            matica(12,is+1) = -sqrt((j+1)/(2*j+1)) * grid%c(ir,+1)
        end select  
    
    end associate
    
  end function matica_mech_chb_viscel_fn
  
end submodule SpheroidalViscel