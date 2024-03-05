submodule (SphericalHarmonics) vcsum
  implicit none ; contains
  
  module pure subroutine vcsum_sub(this, cajm, cbjm, cjm)
    class(T_lateralGrid), intent(in)  :: this
    complex(kind=dbl),    intent(in)  :: cajm(*), cbjm(*)
    complex(kind=dbl),    intent(out) :: cjm(*)
    integer                           :: i
    real(kind=dbl),       allocatable :: pmm(:), pmj(:), pmj1(:), pmj2(:), cosx(:), w(:), sinx(:), grid(:)
    complex(kind=dbl),    allocatable :: cc(:), cr(:), ssym(:), asym(:), sumN(:), sumS(:)
    
    !Array preparation
    allocate( cc(2*this%jms2) ); call zero_carray_sub( 2*this%jms2, cc(1) )
    allocate( cr(this%jms2)   ); call zero_carray_sub(   this%jms2, cr(1) )
      
      call this%scal2scal_jm_to_mj_sub( cajm(1), cc(1), 2, 1 )
      call this%scal2scal_jm_to_mj_sub( cbjm(1), cc(1), 2, 2 )
      
    !Allocating needed memory :: no reallocate for lower stepping
    allocate( pmm(16), pmj(16), pmj1(16), pmj2(16), cosx(16), w(16), sinx(16), ssym(32), &
            & asym(32), sumN(32*this%jmax3), sumS(32*this%jmax3), grid(32*this%nFourier)      )
      
      !Stepping of the algorithm :: 16
      do i = 1, (this%nLegendre/16)*16, 16
        call zero_carray_sub( 32*this%jmax3, sumN(1) )
        call zero_carray_sub( 32*this%jmax3, sumS(1) )

        call this%lege_init_16_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_16_sub( 2, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                         & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_16_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_16_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_16_sub( 1, w(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                        & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 8
      do i = (this%nLegendre/16)*16+1, (this%nLegendre/8)*8, 8
        call zero_carray_sub( 16*this%jmax3, sumN(1) )
        call zero_carray_sub( 16*this%jmax3, sumS(1) )
        
        call this%lege_init_8_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_8_sub( 2, cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1),    &
                                        & pmj(1), ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_8_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_8_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_8_sub( 1, w(1), cosx(1), sinx(1), pmm(1), pmj2(1), pmj1(1), &
                                       & pmj(1), ssym(1), asym(1), cr(1), sumN(1), sumS(1)         )
      end do
      
      !Stepping of the algorithm :: 4
      do i = (this%nLegendre/8)*8+1, (this%nLegendre/4)*4, 4
        call zero_carray_sub( 8*this%jmax3, sumN(1) )
        call zero_carray_sub( 8*this%jmax3, sumS(1) )
        
        call this%lege_init_4_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_4_sub( i, 2, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_4_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_4_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_4_sub( i, 1, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
      !Stepping of the algorithm :: 2
      do i = (this%nLegendre/4)*4+1, this%nLegendre, 2
        call zero_carray_sub( 4*this%jmax3, sumN(1) )
        call zero_carray_sub( 4*this%jmax3, sumS(1) )
        
        call this%lege_init_2_sub( i-1, cosx(1), sinx(1), w(1) )
        
        call this%partial_backward_2_sub( i, 2, cosx(1), pmj2(1), pmj1(1), pmj(1),  &
                                        & ssym(1), asym(1), cc(1), sumN(1), sumS(1) )
        
        call this%grid_op_2_vcsum_sub( grid(1), sumN(1) )
        call this%grid_op_2_vcsum_sub( grid(1), sumS(1) )
        
        call this%partial_forward_2_sub( i, 1, w(1), cosx(1), pmj2(1), pmj1(1), pmj(1), &
                                       & ssym(1), asym(1), cr(1), sumN(1), sumS(1)      )
      end do
      
    deallocate( cc, sumN, sumS, grid, pmm, pmj, pmj1, pmj2, cosx, sinx, w, ssym, asym )
      
      !Rearranging indexing
      call this%scal2scal_mj_to_jm_sub( cr(1), 1, 1, cjm(1), 1, 1)
      
    deallocate(cr)
    
    !Rescaling
    call this%rescale_sub( cjm(1), this%jms )
    
  end subroutine vcsum_sub
  
end submodule vcsum