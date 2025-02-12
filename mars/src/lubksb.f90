     SUBROUTINE LUBKSB(a, indx, b)
     USE PRECISION
     IMPLICIT none

     REAL(RP), INTENT(IN), DIMENSION(:,:) :: a
     INTEGER, INTENT(IN), DIMENSION(:) :: indx
     REAL(RP), INTENT(INOUT), DIMENSION(:) :: b

     INTEGER :: i,ii,ll,j,n
     REAL(RP) :: sum

      n = SIZE(b)
      ii = 0
      do 12 i = 1, n
      ll = indx(i)
      sum = b(ll)
      b(ll) = b(i)
      if (ii .ne. 0) then
      do 11 j = ii, i - 1
      sum = sum - (a(i,j) * b(j))
   11 continue
      else if (sum .ne. 0._RP) then
      ii = i
      end if
      b(i) = sum
   12 continue
      do 14 i = n, 1, -1
      sum = b(i)
      if (i .lt. n) then
      do 13 j = i + 1, n
      sum = sum - (a(i,j) * b(j))
   13 continue
      end if
      b(i) = sum / a(i,i)
   14 continue
      return 
      END SUBROUTINE LUBKSB
