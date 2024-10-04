submodule (solution) temperature
  implicit none; contains
  
  module pure complex(kind=dbl) function temp_fn(this, ir, ijm)
    class(T_solution), intent(in) :: this
    integer,           intent(in) :: ir, ijm
    integer                       :: is
    
    temp_fn = czero
    
    if ( this%inittemp ) then
      is = 3*(ir-1)+1 ; temp_fn = this%temp(is,ijm)
    end if
    
  end function temp_fn
  
  module pure subroutine temp_rr_many1_sub(this, ijm, temp1)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ijm
    complex(kind=dbl), intent(out) :: temp1(:)
    integer                        :: ir
    
    do concurrent ( ir = 1:this%nd+1 )
      temp1(ir) = this%temp(3*(ir-1)+1,ijm)
    end do
    
  end subroutine temp_rr_many1_sub
  
  module pure subroutine temp_jm_many1_sub(this, ir, temp1)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: temp1(:)
    integer                        :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is,ijm)
    end do
    
  end subroutine temp_jm_many1_sub
  
  module pure subroutine temp_jm_many2_sub(this, ir, temp1, temp2)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: temp1(:), temp2(:)
    integer                        :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is  ,ijm)
      temp2(ijm) = this%temp(is+3,ijm)
    end do
    
  end subroutine temp_jm_many2_sub
  
  module pure subroutine temp_jm_many3_sub(this, ir, temp1, temp2, temp3)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: temp1(:), temp2(:), temp3(:)
    integer                        :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is  ,ijm)
      temp2(ijm) = this%temp(is+3,ijm)
      temp3(ijm) = this%temp(is+6,ijm)
    end do
    
  end subroutine temp_jm_many3_sub
  
  module pure subroutine temp_jm_many4_sub(this, ir, temp1, temp2, temp3, temp4)
    class(T_solution), intent(in)  :: this
    integer,           intent(in)  :: ir
    complex(kind=dbl), intent(out) :: temp1(:), temp2(:), temp3(:), temp4(:)
    integer                        :: ijm, is
    
    is = 3*(ir-1)+1
    
    do concurrent ( ijm = 1:this%jms )
      temp1(ijm) = this%temp(is  ,ijm)
      temp2(ijm) = this%temp(is+3,ijm)
      temp3(ijm) = this%temp(is+6,ijm)
      temp4(ijm) = this%temp(is+9,ijm)
    end do
    
  end subroutine temp_jm_many4_sub
  
end submodule temperature