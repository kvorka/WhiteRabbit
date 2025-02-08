submodule (lege_poly) coeffs
  implicit none; contains
  
  module procedure compute_coeffs_sub
    integer :: im, ij, imj, ima
    
    allocate( this%emj((this%jmax+3)*(this%jmax+2)/2) )
      
    do im = 0, this%jmax+1
      do ij = im, this%jmax+1
        this%emj(im*(this%jmax+2)-im*(im+1)/2+ij+1) = sqrt((ij**2-im**2)/(4*ij**2-qone))
      end do
    end do
    
    allocate( this%fmj(3,this%nrma) ) ; ima = 0
    
    do im = 0, this%jmax
      !j = m
        imj = im*(this%jmax+2)-(im-2)*(im+1)/2
        ima = ima+1
        
        if ( im == 0 ) then
          this%fmj(1,ima) = qone / sqrt(4*pi)
        else
          this%fmj(1,ima) = -sqrt( (2*im+qone) / (2*im) )
        end if
      
      if ( im < this%jmax ) then
        do ij = 1, (this%jmax-im)/2
          imj = imj+2
          ima = ima+1
          
          if ( ij == 1 ) then
            this%fmj(1,ima) =               1 / ( this%emj(imj) * this%emj(imj-1) )
            this%fmj(2,ima) = this%emj(imj-1) / ( this%emj(imj)                   )
            this%fmj(3,ima) = zero
          else
            this%fmj(1,ima) =                                                             1 / ( this%emj(imj) * this%emj(imj-1) )
            this%fmj(2,ima) = ( this%emj(imj-1)**2 + this%emj(imj-2)**2                   ) / ( this%emj(imj) * this%emj(imj-1) )
            this%fmj(3,ima) = (                      this%emj(imj-2)    * this%emj(imj-3) ) / ( this%emj(imj) * this%emj(imj-1) )
          end if
        end do
        
        if ( mod((this%jmax-im),2) /= 0 ) then
          imj = imj+2
          ima = ima+1
          
          if ( im == this%jmax-1 ) then
            this%fmj(1,ima) =               1 / ( this%emj(imj) * this%emj(imj-1) )
            this%fmj(2,ima) = this%emj(imj-1) / ( this%emj(imj)                   )
            this%fmj(3,ima) = zero
          else
            this%fmj(1,ima) =                                                             1 / ( this%emj(imj) * this%emj(imj-1) )
            this%fmj(2,ima) = ( this%emj(imj-1)**2 + this%emj(imj-2)**2                   ) / ( this%emj(imj) * this%emj(imj-1) )
            this%fmj(3,ima) = (                      this%emj(imj-2)    * this%emj(imj-3) ) / ( this%emj(imj) * this%emj(imj-1) )
          end if
        end if
      end if
    end do
    
  end procedure compute_coeffs_sub
  
end submodule coeffs
