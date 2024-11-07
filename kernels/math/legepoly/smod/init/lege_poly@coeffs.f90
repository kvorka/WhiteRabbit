submodule (lege_poly) coeffs
  implicit none; contains
  
  module procedure compute_coeffs_sub
    integer :: im, ij, imj
    
    allocate( this%ab(2,(this%jmax+2)*(this%jmax+1)/2) )
      
    do im = 0, this%jmax
      ij = im
        imj = im*(this%jmax+1)-im*(im+1)/2+ij+1
        
        if ( im == 0 ) then
          this%ab(1,imj) = one / s4pi
        else
          this%ab(1,imj) = -sqrt( (2*im+one) / (2*im) )
        end if
      
      if ( im < this%jmax ) then
        ij = im+1
          imj = im*(this%jmax+1)-im*(im+1)/2+ij+1
          
          this%ab(1,imj) = zero
          this%ab(2,imj) = sqrt(2*ij+one)
          
        do ij = im+2, this%jmax
          imj = im*(this%jmax+1)-im*(im+1)/2+ij+1
          
          this%ab(1,imj) = sqrt(         (2*ij+one)*(ij-im-1)*(ij+im-1)/((2*ij-3)*(ij-im)*(ij+im))) / this%ab(2,imj-1)
          this%ab(2,imj) = sqrt((2*ij-1)*(2*ij+one)                    /(         (ij-im)*(ij+im)))
        end do
      end if
    end do
    
  end procedure compute_coeffs_sub
  
end submodule coeffs