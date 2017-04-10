SUBROUTINE READ_PHOTO(name,file_photo,nbrnch,loab,lopr,ionize,enrgI,charge_state,phrct,wcrs,xcrs,brat)

  USE PRECISION
  USE CONSTANTS
  USE SUBS, ONLY : FIND_NAME

  IMPLICIT NONE

  !  .. External Variables

  CHARACTER(len=*), INTENT(IN), DIMENSION(0:) :: name
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

  !  .. Internal Variables

  INTEGER, PARAMETER :: nprmax = 4
  INTEGER :: ncrs, nabs, nbrmax, ndum, nsp
  INTEGER :: na, nb, nw, j, nm, np
  REAL(RP) :: rdum
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: xdum
  CHARACTER(len=256) :: header, cline
  CHARACTER(len=12) :: xname
  CHARACTER(len=12), DIMENSION(6) :: fm
!  CHARACTER(len=1) :: iret

  nsp = SIZE(name)-1

  !  .. Main Section

  OPEN(unit=65, file=file_photo, status='old', action='read')

     READ(65,*) ncrs, nabs, nbrmax
     ALLOCATE(wcrs(ncrs),xcrs(ncrs,nabs),nbrnch(nabs),brat(ncrs,nbrmax,nabs),ionize(nbrmax,nabs),   &
          loab(nabs),lopr(nprmax,nbrmax,nabs),phrct(nbrmax,nabs),enrgI(nbrmax,nabs),                &
          charge_state(nbrmax,nabs))

     ALLOCATE(xdum(ncrs,nabs))
     ionize = .false.
     loab(:) = 0
     lopr(:,:,:) = 0
     enrgI(:,:) = zero
     
     READ(65,'(A)') header                               
     READ(65,'(10ES11.3)') wcrs
     DO na = 1, nabs
        READ(65,'(A12, I4)') xname, ndum
!        WRITE(*,"(A12)") xname
        loab(na) = FIND_NAME(xname,name)
        nbrnch(na) = ndum
        READ(65,'(A)') header
        READ(65,'(10ES11.3)') (xcrs(nw,na),nw=1,ncrs)
        READ(65,'(A)') header
        READ(65,'(10ES11.3)') (xdum(nw,na),nw=1,ncrs)
        DO nb = 1, nbrnch(na)
           READ (65,'(A)') cline
           phrct(nb,na) = cline
           READ (cline,"(5(A12,3X),A12,F6.1)") (fm(j),j=1,6), rdum
           enrgI(nb,na) = rdum
           charge_state(nb,na) = zero
           lopr(1:4,nb,na) = 0 
           np = 0
           DO j = 3, 6
              nm = FIND_NAME(fm(j),name)
              IF ((nm > 0) .and. (nm <= nsp)) THEN
                 np = np + 1
                 lopr(np,nb,na) = nm
              END IF
              IF (nm == nsp) THEN 
                 charge_state(nb,na) = charge_state(nb,na) + 1._RP
                 ionize(nb,na) = .true.
              END IF
           END DO
           WRITE(*,"(A87,F8.3,2X,4A12,F4.0)") phrct(nb,na), enrgI(nb,na), (name(lopr(np,nb,na)),np = 1, 4), charge_state(nb,na)
           READ(65,'(10ES11.3)') (brat(nw,nb,na),nw=1,ncrs)
        END DO
     END DO
     
  CLOSE(unit=65)
!  READ(*,"(A)") iret
  DEALLOCATE(xdum)

  RETURN
END SUBROUTINE READ_PHOTO
