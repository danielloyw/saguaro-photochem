SUBROUTINE NORM_SOL(wav_sol_lres,flx_sol_lres,wav_sol_hres,flx_sol_hres)

  USE PRECISION
  USE CONSTANTS
  USE SUBS, ONLY : LOCATE
  IMPLICIT NONE

  REAL(RP), INTENT(IN), DIMENSION(:) :: wav_sol_lres
  REAL(RP), INTENT(IN), DIMENSION(:) :: flx_sol_lres
  REAL(RP), INTENT(IN), DIMENSION(:) :: wav_sol_hres
  REAL(RP), INTENT(INOUT), DIMENSION(:) :: flx_sol_hres
   
  INTEGER :: nwav_sol_lres, nwav_sol_hres, nw, i1, i2, nf
  REAL(RP) :: sm, rat,delwl, delwh

  nwav_sol_lres = SIZE(wav_sol_lres)
  nwav_sol_hres = SIZE(wav_sol_hres)
  delwl = wav_sol_lres(2)-wav_sol_lres(1)
  delwH = wav_sol_hres(2)-wav_sol_hres(1)
  DO nw = 1, nwav_sol_lres
     IF((wav_sol_lres(nw)-half*delwl >= wav_sol_hres(1)).and.(wav_sol_lres(nw)+half*delwl <= wav_sol_hres(nwav_sol_hres))) THEN
        i1 = LOCATE(wav_sol_hres,wav_sol_lres(nw)-half*delwl)+1
        i2 = LOCATE(wav_sol_hres,wav_sol_lres(nw)+half*delwl)
!
!  Add in endpoints someday
!
        sm = zero
        DO nf = i1, i2
           sm = sm + flx_sol_hres(nf)
        END DO
        rat = (delwl/delwh)*flx_sol_lres(nw)/sm
        DO nf = i1, i2
           flx_sol_hres(nf) = rat*flx_sol_hres(nf)
        END DO
     END IF
  END DO

  RETURN
END SUBROUTINE NORM_SOL
