submodule (SphericalHarmonics) lege_recomb_16
  implicit none; contains
  
  module pure subroutine pmj_backward_recomb_16_sub(nsum, ssym, asym, sumN, sumS)
    integer,           intent(in)    :: nsum
    complex(kind=dbl), intent(in)    :: ssym(16,*), asym(16,*)
    complex(kind=dbl), intent(inout) :: sumN(nsum,*), sumS(nsum,*)
    integer                          :: i1, i2
    
    do concurrent ( i2=1:16, i1=1:nsum )
      sumN(i1,i2) = ssym(i2,i1) + asym(i2,i1)
      sumS(i1,i2) = ssym(i2,i1) - asym(i2,i1)
    end do
    
  end subroutine pmj_backward_recomb_16_sub
  
  module pure subroutine pmj_forward_recomb_16_sub(nsum, weight, sumN, sumS, ssym, asym)
    integer,           intent(in)  :: nsum
    real(kind=dbl),    intent(in)  :: weight(*)
    complex(kind=dbl), intent(in)  :: sumN(nsum,*), sumS(nsum,*)
    complex(kind=dbl), intent(out) :: ssym(16,*), asym(16,*)
    integer                        :: i1, i2
    
    do concurrent ( i2=1:16, i1=1:nsum )
      ssym(i2,i1) = weight(i2) * ( sumN(i1,i2) + sumS(i1,i2) )
      asym(i2,i1) = weight(i2) * ( sumN(i1,i2) - sumS(i1,i2) )
    end do
    
  end subroutine pmj_forward_recomb_16_sub
  
end submodule lege_recomb_16