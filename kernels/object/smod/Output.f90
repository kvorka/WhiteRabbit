submodule (PhysicalObject) Output
  implicit none ; contains
  
  module subroutine vypis_sub(this, filenum, path, quantity)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: filenum
    character(len=*),        intent(in) :: path, quantity
    integer                             :: i, ijm
    
    select case (quantity)
      case ('temperature')
        open(unit=filenum, file=path//'/Temp-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do i = 1, this%nd+1
            write(filenum,*) this%rad_grid%rr(i), this%sol%temp_jm_fn(i)
          end do
        close(filenum)
      
      case ('velocity')
        open(unit=filenum, file=path//'/Velc-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do i = 1, this%nd+1
            write(filenum,*) this%rad_grid%rr(i), this%sol%velocity_jml_fn(i)
          end do
        close(filenum)
      
      case ('flux')
        open(unit=filenum, file=path//'/Flux-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%flux_up(ijm)
          end do
        close(filenum)
      
      case ('topo')
        open(unit=filenum, file=path//'/Topo_dn-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%sol%t_dn(ijm) * this%D_ud
          end do
        close(filenum)
        
        open(unit=filenum, file=path//'/Topo_up-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%sol%t_up(ijm) * this%D_ud
          end do
        close(filenum)
      
      case('shape')
        open(unit=filenum, file=path//'/Shape_dn-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%sol%u_dn(ijm) * this%D_ud
          end do
        close(filenum)

        open(unit=filenum, file=path//'/Shape_up-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%sol%u_up(ijm) * this%D_ud
          end do
        close(filenum)
    end select
    
  end subroutine vypis_sub
  
end submodule Output