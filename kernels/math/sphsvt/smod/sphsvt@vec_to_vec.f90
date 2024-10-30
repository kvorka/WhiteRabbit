submodule (sphsvt) vec_to_vec
  implicit none; contains
  
  module procedure vec2vec_jml_to_jml_sub
    integer :: ijm
    
    do concurrent( ijm = 1:this%jmv )
      cab(cabpadding,ijm) = cjml(ijm)
    end do
    
  end procedure vec2vec_jml_to_jml_sub
  
end submodule vec_to_vec