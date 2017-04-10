SUBROUTINE READ_JVALS(name,file_jvals,nbrnchJ,loabJ,loprJ,phrctJ,srateJ)

  USE PRECISION
  USE CONSTANTS
  USE SUBS, ONLY : FIND_NAME

  IMPLICIT NONE

  !  .. External Variables

  CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: name
  CHARACTER(len=*), INTENT(IN) :: file_jvals
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: nbrnchJ
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: loabJ
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:) :: loprJ
  CHARACTER(len=*), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: phrctJ
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: srateJ

  !  .. Internal Variables

  INTEGER, PARAMETER :: nprmax = 4
  INTEGER :: nabs, nbrmax, ndum, nsp
  CHARACTER(len=256) :: cline
  CHARACTER(len=12) :: xname
  CHARACTER(len=12), DIMENSION(6) :: fm
  REAL(RP) :: rdum
  INTEGER :: na, nb, j, np, nm

  !  .. Main Section

  nsp = SIZE(name)

  OPEN(unit=65, file=file_jvals, status='old', action='read')
     READ(65,*) nabs, nbrmax
     ALLOCATE(nbrnchJ(nabs),loabJ(nabs),loprJ(nprmax,nbrmax,nabs),phrctJ(nbrmax,nabs),srateJ(nbrmax,nabs))
     loabJ(:) = 0
     loprJ(:,:,:) = 0     
     DO na = 1, nabs
        READ(65,'(A12, I4)') xname, ndum
        loabJ(na) = FIND_NAME(xname,name)
        nbrnchJ(na) = ndum
        DO nb = 1, nbrnchJ(na)
           READ (65,'(A)') cline
           phrctJ(nb,na) = cline
           READ (cline,"(5(A12,3X),A12,ES11.3)") (fm(j),j=1,6), rdum
           srateJ(nb,na) = rdum
           np = 0
           DO j = 3, 6
              nm = FIND_NAME(fm(j),name)
              IF ((nm > 0) .and. (nm <= nsp)) THEN
                 np = np + 1
                 loprJ(np,nb,na) = nm
              END IF
           END DO
        END DO
     END DO
     
  CLOSE(unit=65)

  RETURN
END SUBROUTINE READ_JVALS
