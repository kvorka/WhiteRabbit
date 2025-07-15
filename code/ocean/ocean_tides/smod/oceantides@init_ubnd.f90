submodule (oceantides) init_ubnd
  implicit none; contains
  
  module procedure init_ubnd_oceanTides_sub
    integer                        :: i
    complex(kind=dbl), allocatable :: u_t_up(:,:)
    
    allocate( this%v_t_up(4,this%n_iter), u_t_up(4,this%n_iter) )
      
      !! Read the time evolution of the upper boundary vector (spherical part) on degree/order: 2/0, 2/2
      open(unit=1, file='code/ocean/OceanTides/upper_bound_disp.txt', status='old', action='read')
        do
          read(1,*) i, u_t_up(1:4,i)
          if (i == this%n_iter) exit
        end do
      close(1)
      
      !! Periodic time derivative and non-dimensionalization
      this%v_t_up(:,1) = ( u_t_up(:,1) - u_t_up(:,this%n_iter) ) / this%dt / D_ud_ocean
      
      !! Time derivatives inside the period and non-dimensionalization
      do i = 2, this%n_iter
        this%v_t_up(:,i) = ( u_t_up(:,i) - u_t_up(:,i-1) ) / this%dt / D_ud_ocean
      end do
      
    deallocate( u_t_up )
    
  end procedure init_ubnd_oceanTides_sub
  
end submodule init_ubnd