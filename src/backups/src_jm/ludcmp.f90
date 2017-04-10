       SUBROUTINE LUDCMP(a, indx, d)
       USE PRECISION
       IMPLICIT none
       REAL(RP), INTENT(INOUT), DIMENSION(:,:) :: a
       INTEGER, INTENT(OUT) :: indx(:)
       REAL(RP), INTENT(OUT) :: d

       REAL(RP), PARAMETER :: tiny = 1.E-60_RP
       INTEGER i,j,k,imax,n
       REAL(RP) vv(SIZE(a,1)),aamax,sum,dum

       INTRINSIC ABS
    
       n = SIZE(a,1)
       d = 1._RP
       do 12 i = 1, n
       aamax = 0._RP
       do 11 j = 1, n
          if (abs(a(i,j)) > aamax) aamax = ABS(a(i,j))
   11  continue
       if (aamax == 0._RP) stop 'Singular matrix.'
       vv(i) = 1._RP / aamax
   12  continue
       do 19 j = 1, n
         if (j .gt. 1) then
            do 14 i = 1, j - 1
               sum = a(i,j)
               if (i > 1) then
                  do 13 k = 1, i - 1
                     sum = sum - (a(i,k) * a(k,j))
13                continue
                  a(i,j) = sum
               end if
14          continue
         end if
         aamax = 0._RP
         do 16 i = j, n
            sum = a(i,j)
            if (j .gt. 1) then
               do 15 k = 1, j - 1
                  sum = sum - (a(i,k) * a(k,j))
15             continue
               a(i,j) = sum
            end if
            dum = vv(i) * abs(sum)
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            end if
16       continue
         if (j .ne. imax) then
            do 17 k = 1, n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
17          continue
            d = - d
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (j .ne. n) then
            if (a(j,j) == 0._RP) a(j,j) = tiny
            dum = 1._RP / a(j,j)
            do 18 i = j + 1, n
               a(i,j) = a(i,j) * dum
18          continue
          end if
19     continue
       if (a(n,n) == 0._RP) a(n,n) = tiny
       return
       END SUBROUTINE LUDCMP
                         
