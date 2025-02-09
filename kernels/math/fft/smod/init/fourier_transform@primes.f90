submodule (fourier_transform) primes
  implicit none; contains
  
  module procedure prime_adjustement_sub
    integer :: n2, n3, n5, nn0, ntemp, ntemp1
    
    if ( mod(n,2) /= 0 ) then
      nn = n + 1
    else
      nn = n
    end if
    
    nn0 = nn
    
    n2 = 0
    n3 = 0
    n5 = 0
    
    do
      do
        if ( mod(nn,2) == 0 ) then
          nn = nn / 2
          n2 = n2 + 1
        else
          exit
        end if
      end do
      
      do
        if ( mod(nn,3) == 0 ) then
          nn = nn / 3
          n3 = n3 + 1
        else
          exit
        end if
      end do
      
      do
        if ( mod(nn,5) == 0 ) then
          nn = nn / 5
          n5 = n5 + 1
        else
          exit
        end if
      end do
      
      if ( nn == 1 ) then
        exit
      else
        nn = nn + 1
      end if
    end do
    
    nn = 2**n2 * 3**n3 * 5**n5
    ntemp1 = nn
    
    do
      if ( n2 >= 4 ) then
        ntemp = 15 * nn / 16

        if ( ntemp >= nn0 ) then
          nn = ntemp
          n2 = n2 - 4
          n3 = n3 + 1
          n5 = n5 + 1
        end if
      end if
      
      if ( n3 >= 1 ) then
        ntemp = 2 * nn / 3

        if ( ntemp >= nn0 ) then
          nn = ntemp
          n2 = n2 + 1
          n5 = n3 - 1
        end if
      end if
      
      if ( n2 >= 3 ) then
        ntemp = 5 * nn / 8

        if ( ntemp >= nn0 ) then
          nn = ntemp
          n2 = n2 - 3
          n5 = n5 + 1
        end if
      end if
      
      if ( n2 >= 2 ) then
        ntemp = 3 * nn / 4

        if ( ntemp >= nn0 ) then
          nn = ntemp
          n2 = n2 - 2
          n3 = n3 + 1
        end if
      end if
      
      if ( ( n2 >= 2 ) .and. ( n3 >= 1 ) ) then
        ntemp = 5 * nn / 6
        
        if ( ntemp >= nn0 ) then
          nn = ntemp
          n2 = n2 - 1
          n3 = n3 - 1
          n5 = n5 + 1
        end if
      end if
      
      if ( ( n5 >= 1 ) .and. ( n2 >= 2 ) ) then
        ntemp = 9 * nn / 10
        
        if ( ntemp >= nn0 ) then
          nn = ntemp
          n2 = n2 - 1
          n3 = n3 + 1
          n5 = n5 - 1
        end if
      end if
      
      if ( nn == ntemp1 ) then
        exit
      else
        ntemp1 = nn
      end if
    end do
    
  end procedure prime_adjustement_sub
  
end submodule primes