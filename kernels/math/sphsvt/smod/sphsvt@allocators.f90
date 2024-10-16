submodule (sphsvt) allocators
  implicit none; contains
  
  module pure subroutine allocate_scalars_sub(this, ns, cscal)
    class(T_sphsvt),                intent(in)  :: this
    integer,                        intent(in)  :: ns
    complex(kind=dbl), allocatable, intent(out) :: cscal(:)
    
    allocate( cscal(ns*this%jms2) )
      call zero_carray_sub( ns*this%jms2, cscal(1) )
    
  end subroutine allocate_scalars_sub
  
  module pure subroutine allocate_vectors_sub(this, nv, cvec)
    class(T_sphsvt),                intent(in)  :: this
    integer,                        intent(in)  :: nv
    complex(kind=dbl), allocatable, intent(out) :: cvec(:)
    
    allocate( cvec(nv*this%jmv1) )
      call zero_carray_sub( nv*this%jmv1, cvec(1) )
    
  end subroutine allocate_vectors_sub
  
end submodule allocators