   SUBROUTINE TRIDAG (a, b, c, r, u)
   USE PRECISION
   IMPLICIT none

   REAL(RP), INTENT(IN), DIMENSION(:) ::  a, b, c, r
   REAL(RP), INTENT(OUT), DIMENSION(:) ::  u
   INTEGER :: n, j
   REAL(RP), DIMENSION(SIZE(a)) :: gam
   REAL(RP) :: bet
   REAL(RP), PARAMETER :: zero=0.

   n = SIZE(a)
   u = zero

   IF (b(1) == zero) STOP 'tridag: error'
   bet=b(1)
   u(1)=r(1)/bet
   DO j = 2, n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      IF (bet == zero) STOP 'tridag failed'
      u(j)=(r(j)-a(j)*u(j-1))/bet
   END DO
   DO j = n-1, 1, -1
      u(j)=u(j)-gam(j+1)*u(j+1)
   END DO

   RETURN
   END SUBROUTINE TRIDAG
