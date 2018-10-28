      subroutine outerR(R,M,k,n)
      implicit none
      integer n,i,j,k,h
      double precision R(k,n), M(n,n)

      do 10 i=1,n
        do 20 j=1,n
          do 30 h=1,k
            M(i,j) = M(i,j) + R(h,i)*R(h,j)
 30       continue
 20     continue
 10   continue
      return 
      end
