submodule (math) copy
  implicit none; contains
  
  module procedure copy_rarray_sub
    integer :: i
    
    do concurrent  ( i = 1:length )
      arrto(i) = arrfrom(i)
    end do
    
  end procedure copy_rarray_sub
  
end submodule copy