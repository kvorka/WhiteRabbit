submodule (lege_poly) coeffs
  implicit none; contains
  
  module procedure compute_coeffs_sub
    integer                     :: im, ij, imj, ima
    real(kind=qbl), allocatable :: qamj(:), qemj(:)
    
    allocate( this%emj((this%jmax+3)*(this%jmax+2)/2), qemj((this%jmax+3)*(this%jmax+2)/2) )
      
    do im = 0, this%jmax+1
      do ij = im, this%jmax+1
        qemj(im*(this%jmax+2)-im*(im+1)/2+ij+1)     = sqrt((ij**2-im**2)/(4*ij**2-qone))
        this%emj(im*(this%jmax+2)-im*(im+1)/2+ij+1) = sqrt((ij**2-im**2)/(4*ij**2-qone))
      end do
    end do
    
    allocate( this%amj(this%nrma), qamj(this%nrma) ); ima = 0
    
    do im = 0, this%jmax
      !j = m
        imj = im*(this%jmax+2)-(im-2)*(im+1)/2
        ima = ima+1
        
        qamj(ima)     = qone
        this%amj(ima) = real( qamj(ima), kind=dbl )
      
      do ij = 1, (this%jmax-im)/2
        imj = imj+2
        ima = ima+1
        
        qamj(ima)     = qone / ( qemj(imj) * qemj(imj-1) ) / qamj(ima-1)
        this%amj(ima) = real( qamj(ima), kind=dbl )
      end do
      
      if ( mod((this%jmax-im),2) /= 0 ) then
        imj = imj+2
        ima = ima+1
        
        qamj(ima)     = qone / ( qemj(imj) * qemj(imj-1) ) / qamj(ima-1)
        this%amj(ima) = real( qamj(ima), kind=dbl )
      end if
    end do
    
    allocate( this%fmj(2,this%nrma) ) ; ima = 0
    
    do im = 0, this%jmax
      !j = m
        imj = im*(this%jmax+2)-(im-2)*(im+1)/2
        ima = ima+1
        
        if ( im == 0) then
          this%fmj(1,ima) = qone
          this%fmj(2,ima) = qone / sqrt(4*qpi)
        else
          this%fmj(1,ima) = qone
          this%fmj(2,ima) = -sqrt( (2*im+qone) / (2*im) )
        end if
      
      do ij = 1, (this%jmax-im)/2
        imj = imj+2
        ima = ima+1
        
        this%fmj(1,ima) = real( qamj(ima-1)**2, kind=dbl )
        this%fmj(2,ima) = real( ( qemj(imj-1)**2 + qemj(imj-2)**2 ) * qamj(ima-1)**2, kind=dbl )
      end do
      
      if ( mod((this%jmax-im),2) /= 0 ) then
        imj = imj+2
        ima = ima+1
        
        this%fmj(1,ima) = real( qamj(ima-1)**2, kind=dbl )
        this%fmj(2,ima) = real( ( qemj(imj-1)**2 + qemj(imj-2)**2 ) * qamj(ima-1)**2, kind=dbl )
      end if
    end do
    
    deallocate( qemj, qamj )
    
  end procedure compute_coeffs_sub
  
end submodule coeffs