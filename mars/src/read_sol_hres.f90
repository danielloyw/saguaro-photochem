SUBROUTINE READ_SOL_HRES(file_sol_hres,wav_sol_hres,flx_sol_hres)

  USE PRECISION
  USE CONSTANTS
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: file_sol_hres
  REAL(RP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wav_sol_hres
  REAL(RP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: flx_sol_hres

!  INTEGER, PARAMETER :: nwav_hres =  500000
  INTEGER, PARAMETER :: nwav_hres =  7934

  !
  !  .. Read High Resolution Solar Spectrum
  !

!  OPEN(unit=67,file='../data/nrleuv_ref_qs_hr_v2.txt',status='old',action='read')
!     DO nl = 1, 12
!        READ(67,"(A)") header
!     END DO
!     DO nf = 1, nwav_sol_hres
!        READ(67,*) f1, f2
!        wav_sol_hres(nf) = f1
!        flx_sol_hres(nf) = f2/0.05_RP          !  divide by wavelength interval in angstroms to get flux in ph/cm2/s/angstrom
!     END DO
!  CLOSE(unit=67)

  ALLOCATE(wav_sol_hres(nwav_hres),flx_sol_hres(nwav_hres))
  OPEN(unit=63,file=file_sol_hres,status='old',action='read')
     READ(63,*) wav_sol_hres
     READ(63,*) flx_sol_hres
  CLOSE(unit=63)

  RETURN
END SUBROUTINE READ_SOL_HRES
