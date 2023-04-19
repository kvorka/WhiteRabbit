submodule(OutputOceanMod) Nuss
  implicit none

  contains

  subroutine nuss_curve_sub()
    integer         :: n, error
    real(kind=dbl)  :: t, dt, Nuss, Re, zonRe, sumNuss, sumRe, sumzonRe

    n = 0
      sumNuss = 0.0d0
      sumRe = 0.0d0
      sumzonRe = 0.0d0

    open(unit=1, file=path_nuss, status='old', action='read')
      do
        read(1,*,iostat=error) t, dt, Nuss, Re, zonRe; if (error /= 0) exit

        n = n + 1
          sumNuss = sumNuss + Nuss
          sumRe = sumRe + Re
          sumzonRe = sumzonRe + zonRe
      end do
    close(1)

    open(unit=8, file='nuss', status='new', action='write')
      write(8,'(1f4.1)') sumNuss / n
    close(8)
    
  end subroutine nuss_curve_sub

end submodule Nuss