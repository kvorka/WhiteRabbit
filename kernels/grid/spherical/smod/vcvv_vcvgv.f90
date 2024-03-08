submodule (SphericalHarmonics) vcvv_vcvgv
  implicit none ; contains
  
  module pure subroutine grid_op_vcvv_vcvgv_sub(this, nstep, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    integer,                intent(in)    :: nstep
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:,:), gin(:,:,:,:)
    real(kind=dbl), allocatable           :: tmp(:)
    
    call this%fourtrans%exec_c2r_sub(15*nstep, sumNS, grid)
    
    allocate( tmp(3) )
    
    gin(1:3,1:5,1:nstep,1:this%nFourier) => grid(1:15*nstep*this%nFourier)
    gout(1:4,1:nstep,1:this%nFourier)    => grid(1: 4*nstep*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, nstep
        tmp = gin(1:3,1,i2,i)
        
        gout(1,i2,i) = gin(1,2,i2,i) * tmp(1) + gin(2,2,i2,i) * tmp(2) + gin(3,2,i2,i) * tmp(3)
        gout(2,i2,i) = gin(1,3,i2,i) * tmp(1) + gin(2,3,i2,i) * tmp(2) + gin(3,3,i2,i) * tmp(3)
        gout(3,i2,i) = gin(1,4,i2,i) * tmp(1) + gin(2,4,i2,i) * tmp(2) + gin(3,4,i2,i) * tmp(3)
        gout(4,i2,i) = gin(1,5,i2,i) * tmp(1) + gin(2,5,i2,i) * tmp(2) + gin(3,5,i2,i) * tmp(3)
      end do
    end do
    
    deallocate( tmp )
    
    call this%fourtrans%exec_r2c_sub(4*nstep, grid, sumNS)
    
  end subroutine grid_op_vcvv_vcvgv_sub
  
  module pure subroutine vcvv_vcvgv_sub(this, ri, q, dv_r, v, cjm)
    class(T_lateralGrid), intent(in)  :: this
    real(kind=dbl),       intent(in)  :: ri
    complex(kind=dbl),    intent(in)  :: dv_r(*), q(*), v(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: ca(:), cc(:), cr(:)
    
    !Array preparation
    allocate( ca( 5*this%jmv1) ); call zero_carray_sub(  5*this%jmv1, ca(1) ) 
    allocate( cc(15*this%jms2) ); call zero_carray_sub( 15*this%jms2, cc(1) ) 
    allocate( cr( 4*this%jms2) ); call zero_carray_sub(  4*this%jms2, cr(1) )
    
    call this%vec2vec_jml_to_jml_sub( v(1), ca(1), 5, 1 )
    call this%vec2vec_jml_to_jml_sub( q(1), ca(1), 5, 2 )
    call this%gradvec2vec_jmlk_to_jml_sub( ri, v(1), dv_r(1), ca(1), 5, 3 )
    
    call this%vec2scal_jml_to_mj_sub( ca(1), 5, cc(1) )
    
    deallocate(ca)
    
    !Transform
    call this%lege_transform_sub( 4, 15, cc(1), cr(1), grid_op_vcvv_vcvgv_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 4, 1, cjm(1), 4, 1 )
    call this%scal2vecscal_mj_to_jm_sub( cr(1), 4, 2, cjm(1), 4, 2 )
    
    !Cleaning
    deallocate( cc, cr )
    
    !Rescaling
    call this%rescale_sub( cjm(1), 4*this%jms )
    
  end subroutine vcvv_vcvgv_sub
  
end submodule vcvv_vcvgv