submodule (sphsvt) vec_to_vec
  implicit none; contains
  
  module pure subroutine vec2vec_jml_to_jml_sub(this, cjml, cab, ncab, cabpadding)
    class(T_sphsvt),   intent(in)    :: this
    integer,           intent(in)    :: ncab, cabpadding
    complex(kind=dbl), intent(in)    :: cjml(*)
    complex(kind=dbl), intent(inout) :: cab(ncab,*)
    integer                          :: ijm
    
    do concurrent( ijm = 1:this%jmv )
      cab(cabpadding,ijm) = cjml(ijm)
    end do
    
  end subroutine vec2vec_jml_to_jml_sub
  
end submodule vec_to_vec