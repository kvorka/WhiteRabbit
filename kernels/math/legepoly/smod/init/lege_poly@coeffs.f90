submodule (lege_poly) coeffs
  implicit none; contains
  
  module procedure compute_coeffs_sub
    integer :: im, ij, imj
    
    allocate( this%amj((this%jmax+2)*(this%jmax+1)/2), &
            & this%bmj((this%jmax+2)*(this%jmax+1)/2), &
            & this%cmm(0:this%jmax) )
      
    do im = 0, this%jmax
      if ( im == 0 ) then
        this%cmm(im) = one / s4pi
      else
        this%cmm(im) = -sqrt( (2*im+one) / (2*im) )
      end if
      
      ij = im
        imj = im*(this%jmax+1)-im*(im+1)/2+ij+1
        
        this%amj(imj) = one
        this%bmj(imj) = one
      
      do ij = im+1, this%jmax
        imj = im*(this%jmax+1)-im*(im+1)/2+ij+1
        
        this%amj(imj) = sqrt((2*ij-1)*(2*ij+one)                    /(         (ij-im)*(ij+im)))
        this%bmj(imj) = sqrt(         (2*ij+one)*(ij-im-1)*(ij+im-1)/((2*ij-3)*(ij-im)*(ij+im))) / this%amj(imj-1)
      end do
    end do
    
  end procedure compute_coeffs_sub
  
end submodule coeffs