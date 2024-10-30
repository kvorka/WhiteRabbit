submodule (solution) temperature
  implicit none; contains
  
  module procedure temp_fn
    integer :: is
    
    temp_fn = czero
    
    if ( this%inittemp ) then
      is = 3*(ir-1)+1 ; temp_fn = this%temp(is,ijm)
    end if
    
  end procedure temp_fn
  
  module procedure temp_rr_many1_sub
    integer :: ir
    
    do concurrent ( ir = 1:this%nd+1 )
      temp1(ir) = this%temp(3*(ir-1)+1,ijm)
    end do
    
  end procedure temp_rr_many1_sub
  
  module procedure temp_jm_many1_sub
    integer :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is,ijm)
    end do
    
  end procedure temp_jm_many1_sub
  
  module procedure temp_jm_many2_sub
    integer :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is  ,ijm)
      temp2(ijm) = this%temp(is+3,ijm)
    end do
    
  end procedure temp_jm_many2_sub
  
  module procedure temp_jm_many3_sub
    integer :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is  ,ijm)
      temp2(ijm) = this%temp(is+3,ijm)
      temp3(ijm) = this%temp(is+6,ijm)
    end do
    
  end procedure temp_jm_many3_sub
  
  module procedure temp_jm_many4_sub
    integer :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is  ,ijm)
      temp2(ijm) = this%temp(is+3,ijm)
      temp3(ijm) = this%temp(is+6,ijm)
      temp4(ijm) = this%temp(is+9,ijm)
    end do
    
  end procedure temp_jm_many4_sub
  
end submodule temperature