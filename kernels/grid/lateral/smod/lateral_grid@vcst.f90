submodule (lateral_grid) vcst
  implicit none ; contains
  
  pure subroutine grid_op_vcst_sub(step, nfour, grid)
    integer,                intent(in)    :: nfour, step
    real(kind=dbl), target, intent(inout) :: grid(*)
    integer                               :: i1, i2, i3
    real(kind=dbl), pointer               :: gin(:,:,:), gout(:,:,:)
    real(kind=dbl), allocatable           :: tmp(:), tmp1(:)
    
    gin(1:step,1:6,1:nfour)  => grid(1:6*step*nfour)
    gout(1:step,1:5,1:nfour) => grid(1:5*step*nfour)
    
    allocate( tmp(step), tmp1(step) )
    
    do i1 = 1, nfour
      tmp = gin(1:step,1,i1)
      
      do i3 = 1, 5
        tmp1 = gin(1:step,i3+1,i1)
        
        do concurrent ( i2 = 1:step )
          gout(i2,i3,i1) = tmp(i2) * tmp1(i2)
        end do
      end do
    end do
    
    deallocate( tmp, tmp1 )
    
  end subroutine grid_op_vcst_sub
  
  module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
    complex(kind=dbl),    intent(out) :: cjml2(*)
    complex(kind=dbl),    allocatable :: cc(:), cr(:)
    
    !Array preparation
    call this%reindexing%allocate_scalars_sub( 6, cc )
    call this%reindexing%allocate_scalars_sub( 5, cr )
    
    call this%reindexing%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 6, 1 )
    call this%reindexing%devtens2scal_jml2_to_mj_sub( cbjml2(1), cc(1), 6, 2 )
    
    !Transform
    call this%transform_sub( 5, 6, cc(1), cr(1), grid_op_vcst_sub )
    
    !Rearranging indexing
    call this%reindexing%scal2devtens_mj_to_jml2_sub( cr(1), 5, 1, cjml2(1) )
    
    !Cleaning
    deallocate( cc, cr )
    
  end subroutine vcst_sub
  
end submodule vcst