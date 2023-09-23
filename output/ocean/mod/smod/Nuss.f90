submodule(OutputOceanMod) Nuss
  implicit none
  
  contains
  
  subroutine nuss_curve_sub()
    integer         :: n, error
    real(kind=dbl)  :: t, dt, Nuss, Re, sumNuss, sumRe
    
    n = 0
      sumNuss = 0._dbl
      sumRe   = 0._dbl
    
    open(unit=1, file=path_nuss, status='old', action='read')
    
    do
      read(1,*,iostat=error) t, dt, Nuss, Re
      
      if ( error /= 0 ) then
        exit
      else if ( t > tNuss ) then
        n = n + 1
          sumNuss = sumNuss + Nuss
          sumRe   = sumRe + Re
      end if
    end do
    
    close(1)
    
    open(unit=8, file='nuss', status='new', action='write')
      write(8,'(1f4.1)') sumNuss / n , sumRe / n
    close(8)
    
  end subroutine nuss_curve_sub

end submodule Nuss