SUBROUTINE READ_SOL_BND(file_sol,wav_sol,flx_sol)

  USE PRECISION
  USE CONSTANTS
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: file_sol
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: wav_sol
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: flx_sol
  INTEGER :: nl, nw, nwav_sol
  REAL(RP) :: f1, f2, f3, f4, f5, eng
  CHARACTER(len=128) :: header

  OPEN(unit=66,file=file_sol,status='old',action='read')
     DO nl = 1, 2
        READ(66,"(A)") header
     END DO
     nw = 0
     DO
        nw = nw + 1
        READ(66,*,END=10) f1, f2
    END DO
10  CONTINUE
    nwav_sol = nw-1
    ALLOCATE(wav_sol(nwav_sol),flx_sol(nwav_sol))
    REWIND(66)
    DO nl = 1, 2
        READ(66,"(A)") header
     END DO
     DO nw = 1, nwav_sol
        READ(66,*) f1, f2
        wav_sol(nw) = f1                    
        flx_sol(nw) = f2             
     END DO
  CLOSE(unit=66)

  RETURN
END SUBROUTINE READ_SOL_BND
