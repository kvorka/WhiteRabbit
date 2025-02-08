submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure c2r_mj_to_mj_sub
    integer :: i, m, j, mj, ma
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = cab(i,mj+1) * this%emj(mj+1)
          rcab(i,2,ma) = cab(i,mj)
        end do
      
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = this%emj(mj) * cab(i,mj-1) + this%emj(mj+1) * cab(i,mj+1)
          rcab(i,2,ma) =                cab(i,mj)
        end do
      end do
      
      !j == this%jmax
      if ( mod((this%jmax),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = this%emj(mj) * cab(i,mj-1)
          rcab(i,2,ma) =                cab(i,mj)
        end do
      
      else
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = this%emj(mj+1) * cab(i,mj)
        end do
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = cab(i,mj+1) * this%emj(mj+m+1)
          rcab(i,2,ma) = cab(i,mj)
        end do
      
      do j = 1, (this%jmax-1-m)/2
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = this%emj(mj+m) * cab(i,mj-1) + this%emj(mj+m+1) * cab(i,mj+1)
          rcab(i,2,ma) =                  cab(i,mj)
        end do
      end do
      
      !j == this%jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = this%emj(mj+m) * cab(i,mj-1)
          rcab(i,2,ma) =                  cab(i,mj)
        end do
      
      else
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(i,1,ma) = this%emj(mj+m+1) * cab(i,mj)
        end do
      end if
    end do
    
    m = this%jmax
      !j == m
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(i,2,ma) = cab(i,mj)
        end do
      
  end procedure c2r_mj_to_mj_sub
  
end submodule c2r