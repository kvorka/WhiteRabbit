submodule (SphericalHarmonics) vcsum
  implicit none ; contains
  
  module pure subroutine grid_op_vcsum_sub(this, nstep, grid, sumNS)
    class(T_lateralGrid),   intent(in)    :: this
    integer,                intent(in)    :: nstep
    real(kind=dbl), target, intent(out)   :: grid(*)
    complex(kind=dbl),      intent(inout) :: sumNS(*)
    integer                               :: i, i2
    real(kind=dbl), pointer               :: gout(:,:), gin(:,:,:)
    
    call this%fourtrans%exec_c2r_sub(2*nstep, sumNS, grid)
    
    gin(1:2,1:nstep,1:this%nFourier) => grid(1:2*nstep*this%nFourier)
    gout(1:nstep,1:this%nFourier)    => grid(1:  nstep*this%nFourier)
    
    do i = 1, this%nFourier
      do i2 = 1, nstep
        gout(i2,i) = gin(1,i2,i) * gin(2,i2,i)
      end do
    end do
    
    call this%fourtrans%exec_r2c_sub(nstep, grid, sumNS)
    
  end subroutine grid_op_vcsum_sub
  
  module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    allocate( cc(2*this%jms2) ); call zero_carray_sub( 2*this%jms2, cc(1) )
    allocate( cr(this%jms2)   ); call zero_carray_sub(   this%jms2, cr(1) )
    
    call this%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 2, 1 )
    call this%scal2scal_jm_to_mj_sub( cbjm(1), cc(1), 2, 2 )
    
    !Transform
    call this%lege_transform_sub( 1, 2, cc(1), cr(1), grid_op_vcsum_sub )
    
    !Rearranging indexing
    call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1)
    
    !Cleaning
    deallocate( cr, cc )
    
    !Rescaling
    call this%rescale_sub( cjm(1), this%jms )
    
  end subroutine vcsum_sub
  
end submodule vcsum