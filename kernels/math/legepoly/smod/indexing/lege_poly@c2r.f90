submodule (lege_poly) c2r
  implicit none; contains
  
  module procedure c2r_mj_to_mj_sub
    integer :: i, m, j, mj, ma
    
    m = 0
      !j == m
        ma = 1
        mj = 1
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = cab(i,mj+1)%re * this%emj(mj+1)
          rcab(2,i,1,ma) = cab(i,mj+1)%im * this%emj(mj+1)
          rcab(1,i,2,ma) = cab(i,mj)%re
          rcab(2,i,2,ma) = cab(i,mj)%im
        end do
      
      do j = 1, (this%jmax-1)/2
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = this%emj(mj) * cab(i,mj-1)%re + this%emj(mj+1) * cab(i,mj+1)%re
          rcab(2,i,1,ma) = this%emj(mj) * cab(i,mj-1)%im + this%emj(mj+1) * cab(i,mj+1)%im
          rcab(1,i,2,ma) =                cab(i,mj)%re
          rcab(2,i,2,ma) =                cab(i,mj)%im
        end do
      end do
      
      !j == this%jmax
      if ( mod((this%jmax),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = this%emj(mj) * cab(i,mj-1)%re
          rcab(2,i,1,ma) = this%emj(mj) * cab(i,mj-1)%im
          rcab(1,i,2,ma) =                cab(i,mj)%re
          rcab(2,i,2,ma) =                cab(i,mj)%im
        end do
      
      else
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = this%emj(mj+1) * cab(i,mj)%re
          rcab(2,i,1,ma) = this%emj(mj+1) * cab(i,mj)%im
        end do
      end if
    
    do m = 1, this%jmax-1
      !j == m
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = cab(i,mj+1)%re * this%emj(mj+m+1)
          rcab(2,i,1,ma) = cab(i,mj+1)%im * this%emj(mj+m+1)
          rcab(1,i,2,ma) = cab(i,mj)%re
          rcab(2,i,2,ma) = cab(i,mj)%im
        end do
      
      do j = 1, (this%jmax-1-m)/2
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = this%emj(mj+m) * cab(i,mj-1)%re + this%emj(mj+m+1) * cab(i,mj+1)%re
          rcab(2,i,1,ma) = this%emj(mj+m) * cab(i,mj-1)%im + this%emj(mj+m+1) * cab(i,mj+1)%im
          rcab(1,i,2,ma) =                  cab(i,mj)%re
          rcab(2,i,2,ma) =                  cab(i,mj)%im
        end do
      end do
      
      !j == this%jmax
      if ( mod((this%jmax-m),2) == 0 ) then
        ma = ma+1
        mj = mj+2
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = this%emj(mj+m) * cab(i,mj-1)%re
          rcab(2,i,1,ma) = this%emj(mj+m) * cab(i,mj-1)%im
          rcab(1,i,2,ma) =                  cab(i,mj)%re
          rcab(2,i,2,ma) =                  cab(i,mj)%im
        end do
      
      else
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,1,ma) = this%emj(mj+m+1) * cab(i,mj)%re
          rcab(2,i,1,ma) = this%emj(mj+m+1) * cab(i,mj)%im
        end do
      end if
    end do
    
    m = this%jmax
      !j == m
        ma = ma+1
        mj = mj+1
        
        do concurrent ( i = 1:ncab )
          rcab(1,i,2,ma) = cab(i,mj)%re
          rcab(2,i,2,ma) = cab(i,mj)%im
        end do
      
  end procedure c2r_mj_to_mj_sub
  
end submodule c2r