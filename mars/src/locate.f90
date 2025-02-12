! This function locates the value x in a sorted list xtab. The earlier index will be picked.

   FUNCTION LOCATE(xtab, x)
   USE PRECISION
   IMPLICIT NONE
   REAL(RP), DIMENSION(:), INTENT(IN) :: xtab
   REAL(RP), INTENT(IN) :: x
   INTEGER :: LOCATE
   INTEGER :: n,jl,jm,ju
   LOGICAL :: ascnd ! is list ascending?
   n=size(xtab)
   ascnd = (xtab(n) >= xtab(1))
   jl=0
   ju=n+1
   do
      if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd .eqv. (x >= xtab(jm))) then
         jl=jm
      else
         ju=jm
      end if
   end do
   if (x == xtab(1)) then
      LOCATE=1
   else if (x == xtab(n)) then
      LOCATE=n-1 !?
   else
      LOCATE=jl
   end if
   END FUNCTION LOCATE
