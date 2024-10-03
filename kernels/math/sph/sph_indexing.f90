module sph_indexing
  use math
  implicit none; public; contains
  
  pure integer function jm(ij, im)
    integer, intent(in) :: ij, im
    
    jm = ij*(ij+1)/2+im+1
    
  end function jm
  
  pure integer function jml(ij, im, il)
    integer, intent(in) :: ij, im, il
    
    jml = 3*(ij*(ij+1)/2+im)+il
    
  end function jml
  
  pure integer function jml2(ij, im, il)
    integer, intent(in) :: ij, im, il
    
    jml2 = 5*(ij*(ij+1)/2+im)+il-1
    
  end function jml2
  
end module sph_indexing