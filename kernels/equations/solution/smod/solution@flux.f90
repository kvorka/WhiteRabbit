submodule (solution) flux
  implicit none; contains

  module procedure flux_fn
    integer :: is
    
    flux_fn = czero ; is = 3*(ir-1)+1
    
    select case (il)
      case (-1)
        flux_fn = this%temp(is+1,ijm)
      case ( 0)
        flux_fn = czero
      case (+1)
        flux_fn = this%temp(is+2,ijm)
    end select
    
  end procedure flux_fn
  
  module procedure flux_r_many1_sub
    integer :: ir, is
    
    do concurrent ( ir = 1:this%nd )
      is = 3*(ir-1)+2
      
      flux1(1,ir) = this%temp(is,ijm)
      flux1(2,ir) = czero
      flux1(3,ir) = this%temp(is+1,ijm)
    end do
    
  end procedure flux_r_many1_sub
  
  module procedure flux_jml_many1_sub
    integer :: is, ijm, ijml
    
    is = 3*(ir-1)+2
    
    ijm = 1
      ijml = 1
        flux1(1) = this%temp(is+1,1)
    
    do concurrent ( ijm = 2:this%jms )
      ijml = 3*(ijm-1)-1
       
      flux1(ijml  ) = this%temp(is,ijm)
      flux1(ijml+1) = czero
      flux1(ijml+2) = this%temp(is+1,ijm)
    end do
    
  end procedure flux_jml_many1_sub
  
  module procedure flux_jml_many2_sub
    integer :: is, ijm, ijml
    
    is = 3*(ir-1)+2
    
    ijm = 1
      ijml = 1
        flux1(1) = this%temp(is+1,1)
        flux2(1) = this%temp(is+4,1)
    
    do concurrent ( ijm = 2:this%jms )
      ijml = 3*(ijm-1)-1
       
      flux1(ijml  ) = this%temp(is,ijm)
      flux1(ijml+1) = czero
      flux1(ijml+2) = this%temp(is+1,ijm)
      
      flux2(ijml  ) = this%temp(is+3,ijm)
      flux2(ijml+1) = czero
      flux2(ijml+2) = this%temp(is+4,ijm)
    end do
    
  end procedure flux_jml_many2_sub
  
  module procedure flux_jml_many3_sub
    integer :: is, ijm, ijml
    
    is = 3*(ir-1)+2
    
    ijm = 1
      ijml = 1
        flux1(1) = this%temp(is+1,1)
        flux2(1) = this%temp(is+4,1)
        flux3(1) = this%temp(is+7,1)
    
    do concurrent ( ijm = 2:this%jms )
      ijml = 3*(ijm-1)-1
       
      flux1(ijml  ) = this%temp(is,ijm)
      flux1(ijml+1) = czero
      flux1(ijml+2) = this%temp(is+1,ijm)
      
      flux2(ijml  ) = this%temp(is+3,ijm)
      flux2(ijml+1) = czero
      flux2(ijml+2) = this%temp(is+4,ijm)
      
      flux3(ijml  ) = this%temp(is+6,ijm)
      flux3(ijml+1) = czero
      flux3(ijml+2) = this%temp(is+7,ijm)
    end do
    
  end procedure flux_jml_many3_sub
  
  module procedure flux_jml_many4_sub
    integer :: is, ijm, ijml
    
    is = 3*(ir-1)+2
    
    ijm = 1
      ijml = 1
        flux1(1) = this%temp(is+ 1,1)
        flux2(1) = this%temp(is+ 4,1)
        flux3(1) = this%temp(is+ 7,1)
        flux4(1) = this%temp(is+10,1)
    
    do concurrent ( ijm = 2:this%jms )
      ijml = 3*(ijm-1)-1
       
      flux1(ijml  ) = this%temp(is,ijm)
      flux1(ijml+1) = czero
      flux1(ijml+2) = this%temp(is+1,ijm)
      
      flux2(ijml  ) = this%temp(is+3,ijm)
      flux2(ijml+1) = czero
      flux2(ijml+2) = this%temp(is+4,ijm)
      
      flux3(ijml  ) = this%temp(is+6,ijm)
      flux3(ijml+1) = czero
      flux3(ijml+2) = this%temp(is+7,ijm)
      
      flux4(ijml  ) = this%temp(is+9,ijm)
      flux4(ijml+1) = czero
      flux4(ijml+2) = this%temp(is+10,ijm)
    end do
    
  end procedure flux_jml_many4_sub
  
end submodule flux