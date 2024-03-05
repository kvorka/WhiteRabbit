submodule (SphericalHarmonics) vcst
  implicit none ; contains
  
  module pure subroutine vcst_sub(this, cajm, cbjml2, cjml2)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjml2(*)
    complex(kind=dbl),    intent(out) :: cjml2(*)
    integer                           :: i
    real(kind=dbl),       allocatable :: pmj(:), pmj1(:), pmj2(:), cosx(:), w(:), grid(:)
    complex(kind=dbl),    allocatable :: ssym(:), asym(:), sumN(:), sumS(:), cc(:), cr(:)
    
    !Array preparation
    allocate( cc(6*this%jms2) ); call zero_carray_sub( 6*this%jms2, cc(1) )
    allocate( cr(5*this%jms2) ); call zero_carray_sub( 5*this%jms2, cr(1) )
      
      call this%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 6, 6 )
      call this%devtens2scal_jml2_to_mj_sub( cbjml2(1), cc(1), 6, 1 )
    
    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmj(16), pmj1(16), pmj2(16), cosx(16), w(16), ssym(96), asym(96), &
            & sumN(96*this%jmax3), sumS(96*this%jmax3), grid(96*this%nFourier)  )
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call this%partial_backward_16_sub( i, 6, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                         & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_16_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( i, 5, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                        & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call this%partial_backward_8_sub( i, 6, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_8_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( i, 5, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call this%partial_backward_4_sub( i, 6, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_4_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( i, 5, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call this%partial_backward_2_sub( i, 6, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
          
        call this%grid_op_2_vcst_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcst_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( i, 5, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
    deallocate( cc, sumN, sumS, grid, pmj, pmj1, pmj2, cosx, w, ssym, asym )
      
      !Rearranging indexing
      call this%scal2devtens_mj_to_jml2_sub( cr(1), 5, 1, cjml2(1) )
      
    deallocate( cr )
    
    !Rescaling
    call this%rescale_sub( cjml2(1), this%jmt )
    
  end subroutine vcst_sub
  
end submodule vcst