submodule(Solution) Solution_init
  implicit none

  contains

  subroutine init_stemp_sub(this)
    class(T_solution), intent(inout) :: this
      
    allocate( this%temp(3*this%nd+1, this%jms) ) ; this%temp = cmplx(0._dbl, 0._dbl, kind=dbl)
    
  end subroutine init_stemp_sub

  subroutine init_storr_sub(this)
    class(T_solution), intent(inout) :: this
      
    allocate( this%torr(3*this%nd+1, this%jms) ) ; this%torr = cmplx(0._dbl, 0._dbl, kind=dbl)
    
  end subroutine init_storr_sub

  subroutine init_smech_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%mech(6*this%nd+2,this%jms) ) ; this%mech = cmplx(0._dbl, 0._dbl, kind=dbl)
    
  end subroutine init_smech_sub

  subroutine init_layers_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%u_up(this%jms) ); this%u_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( this%u_dn(this%jms) ); this%u_dn = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( this%u_I2(this%jms) ); this%u_I2 = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( this%u_C(this%jms)  ); this%u_C  = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    allocate( this%t_dn(this%jms) ); this%t_dn = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( this%t_up(this%jms) ); this%t_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    
    allocate( this%v_up(this%jms) ); this%v_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( this%v_dn(this%jms) ); this%v_dn = cmplx(0._dbl, 0._dbl, kind=dbl)
    
  end subroutine init_layers_sub

  subroutine init_layer_u_sub(this)
    class(T_solution), intent(inout) :: this
    
    allocate( this%t_up(this%jms) ); this%t_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    allocate( this%u_up(this%jms) ); this%u_up = cmplx(0._dbl, 0._dbl, kind=dbl)
    
  end subroutine init_layer_u_sub

end submodule Solution_init