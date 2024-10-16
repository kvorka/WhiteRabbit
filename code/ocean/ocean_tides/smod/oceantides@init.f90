submodule (oceantides) init
  implicit none; contains
  
  module subroutine init_oceanTides_sub(this)
    class(T_oceanTides), intent(inout) :: this
    integer                            :: i
    complex(kind=dbl), allocatable     :: u201(:), u203(:), u221(:), u223(:)
    
    call this%init_ocean_sub()
    
    this%dt = 2 * pi / this%n_iter ; this%number_of_periods = 0
    
    call this%init_eq_torr_sub()
      allocate( this%ntorr(this%jms,2:this%nd) )
      this%ntorr = czero
      
      call this%prepare_mat_torr_sub( ijstart=1 , ijend=this%jmax )
    
    call this%init_eq_mech_sub()
      allocate( this%nsph1(this%jms,2:this%nd), this%nsph2(this%jms,2:this%nd) )
      this%nsph1 = czero
      this%nsph2 = czero
      
      call this%prepare_mat_mech_sub( ijstart=1 , ijend=this%jmax )
    
    allocate( this%v201(this%n_iter), this%v203(this%n_iter), this%v221(this%n_iter), this%v223(this%n_iter) )
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
    
  end subroutine init_oceanTides_sub
  
end submodule init