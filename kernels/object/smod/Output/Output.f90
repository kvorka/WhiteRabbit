submodule (PhysicalObject) Output
  implicit none ; contains
  
  module subroutine vypis_sub(this, filenum, path, quantity)
    class(T_physicalObject), intent(in) :: this
    integer,                 intent(in) :: filenum
    character(len=*),        intent(in) :: path, quantity
    integer                             :: ir, ijm
    complex(kind=dbl),      allocatable :: field(:)
    
    select case (quantity)
      case ('temperature')
        open(unit=filenum, file=path//'/Temp-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          allocate( field(this%jms) )
            do ir = 1, this%nd+1
              call this%temp_rr_ijm_sub(ir, field); write(filenum,*) this%rad_grid%rr(ir), field
            end do
          deallocate( field )
        close(filenum)
      
      case ('velocity')
        open(unit=filenum, file=path//'/Velc-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          allocate( field(this%jmv) )
            do ir = 1, this%nd+1
              call this%v_rr_ijml_sub(ir, field); write(filenum,*) this%rad_grid%rr(ir), field
            end do
          deallocate( field )
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
            write(filenum,*) ijm, this%bnd%t_dn(ijm) * this%D_ud
          end do
        close(filenum)
        
        open(unit=filenum, file=path//'/Topo_up-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%bnd%t_up(ijm) * this%D_ud
          end do
        close(filenum)
      
      case('shape')
        open(unit=filenum, file=path//'/Shape_dn-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%bnd%u_dn(ijm) * this%D_ud
          end do
        close(filenum)

        open(unit=filenum, file=path//'/Shape_up-'//trim(adjustl(int2str_fn(this%poc)))//'.dat', status='new', action='write')
          do ijm = 1, this%jms
            write(filenum,*) ijm, this%bnd%u_up(ijm) * this%D_ud
          end do
        close(filenum)
    end select
    
  end subroutine vypis_sub
  
end submodule Output