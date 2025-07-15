submodule (oceantides) init_ubnd
  implicit none; contains
  
  module procedure init_ubnd_oceanTides_sub
    integer                        :: i
    complex(kind=dbl), allocatable :: u201(:), u203(:), u221(:), u223(:)
    
    allocate( this%v201(this%n_iter), &
            & this%v203(this%n_iter), &
            & this%v221(this%n_iter), &
            & this%v223(this%n_iter)  )
    
    allocate( u201(this%n_iter), u203(this%n_iter), u221(this%n_iter), u223(this%n_iter) )
      open(unit=1, file='code/ocean/OceanTides/upper_bound_disp.txt', status='old', action='read')
        do
          read(1,*) i, u201(i), u203(i), u221(i), u223(i)
          if (i == this%n_iter) exit
        end do
      close(1)
      
      u201(:) = u201(:) / D_ud_ocean; u203(:) = u203(:) / D_ud_ocean
      u221(:) = u221(:) / D_ud_ocean; u223(:) = u223(:) / D_ud_ocean
      
      this%v201(1) = ( u201(1) - u201(this%n_iter) ) / this%dt; this%v203(1) = ( u203(1) - u203(this%n_iter) ) / this%dt
      this%v221(1) = ( u221(1) - u221(this%n_iter) ) / this%dt; this%v223(1) = ( u223(1) - u223(this%n_iter) ) / this%dt
      
      do i = 2, this%n_iter
        this%v201(i) = ( u201(i) - u201(i-1) ) / this%dt; this%v203(i) = ( u203(i) - u203(i-1) ) / this%dt
        this%v221(i) = ( u221(i) - u221(i-1) ) / this%dt; this%v223(i) = ( u223(i) - u223(i-1) ) / this%dt
      end do
      
    deallocate( u201, u203, u221, u223 )
    
  end procedure init_ubnd_oceanTides_sub
  
end submodule init_ubnd