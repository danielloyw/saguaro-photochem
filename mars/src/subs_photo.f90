MODULE SUBS_PHOTO


  INTERFACE
     SUBROUTINE READ_PHOTO(name,file_photo,nbrnch,loab,lopr,ionize,enrgI,charge_state,phrct,wcrs,xcrs,brat)
       USE PRECISION
       USE CONSTANTS
       USE SUBS, ONLY : FIND_NAME
       IMPLICIT NONE
       CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: name
       CHARACTER(len=*), INTENT(IN) :: file_photo
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: nbrnch
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: loab
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:) :: lopr
       LOGICAL, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: ionize
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: enrgI
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: charge_state
       CHARACTER(len=*), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: phrct
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: wcrs
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: xcrs
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:) :: brat
     END SUBROUTINE READ_PHOTO
  END INTERFACE

  INTERFACE
     SUBROUTINE READ_PHOTOB(name,nbrnchB,loabB,loprB,ionizeB,enrgIB,charge_stateB,phrctB,wcrsB,wcrsB_low,delwB,xcrsB,bratB)
       USE PRECISION
       USE CONSTANTS
       USE SUBS, ONLY : FIND_NAME, INTRP
       IMPLICIT NONE
       CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: name
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:):: nbrnchB
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:):: loabB
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:):: loprB
       LOGICAL, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: ionizeB
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: enrgIB
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: charge_stateB
       CHARACTER(len=*), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: phrctB
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: wcrsB, wcrsB_low
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: delwB
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: xcrsB
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:) :: bratB
     END SUBROUTINE READ_PHOTOB
  END INTERFACE

  INTERFACE
     SUBROUTINE READ_JVALS(name,file_jvals,nbrnchJ,loabJ,loprJ,phrctJ,srateJ)
       USE PRECISION
       USE CONSTANTS
       USE SUBS, ONLY : FIND_NAME
       IMPLICIT NONE
       CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: name
       CHARACTER(len=*), INTENT(IN) :: file_jvals
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: nbrnchJ
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: loabJ
       INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:) :: loprJ
       CHARACTER(len=*), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: phrctJ
       REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: srateJ
     END SUBROUTINE READ_JVALS
  END INTERFACE

  INTERFACE
     SUBROUTINE READ_SOL_HRES(file_sol_hres,wav_sol_hres,flx_sol_hres)
       USE PRECISION
       USE CONSTANTS
       IMPLICIT NONE
       CHARACTER(len=*), INTENT(IN) :: file_sol_hres
       REAL(RP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wav_sol_hres
       REAL(RP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: flx_sol_hres
     END SUBROUTINE READ_SOL_HRES
  END INTERFACE

  INTERFACE
     SUBROUTINE READ_SOL(file_sol,wav_sol,flx_sol)
       USE PRECISION
       USE CONSTANTS
       IMPLICIT NONE
       CHARACTER(len=*), INTENT(IN) :: file_sol
       REAL(RP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wav_sol
       REAL(RP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: flx_sol
     END SUBROUTINE READ_SOL
  END INTERFACE

  INTERFACE
     SUBROUTINE NORM_SOL(wav_sol_lres,flx_sol_lres,wav_sol_hres,flx_sol_hres)
       USE PRECISION
       USE CONSTANTS
       USE SUBS, ONLY : LOCATE
       IMPLICIT NONE
       REAL(RP), INTENT(IN), DIMENSION(:) :: wav_sol_lres
       REAL(RP), INTENT(IN), DIMENSION(:) :: flx_sol_lres
       REAL(RP), INTENT(IN), DIMENSION(:) :: wav_sol_hres
       REAL(RP), INTENT(INOUT), DIMENSION(:) :: flx_sol_hres
     END SUBROUTINE NORM_SOL
  END INTERFACE

END MODULE SUBS_PHOTO
