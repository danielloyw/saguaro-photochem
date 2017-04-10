SUBROUTINE READ_MOLECULES

  !  .. Modules

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_NAME

  !  .. Local Variables Declarations

  IMPLICIT NONE
  INTEGER :: nx1, nx2, nm, nx, ndum, n, ichk, neq
  CHARACTER(len=128) :: header

  !  .. Read neutral and ionize molecules

  OPEN(unit=61,file='nmolecules.dat',status='old',action='read')
  OPEN(unit=62,file='imolecules.dat',status='old',action='read')

  READ(61,*) neutrmax
  READ(62,*) nionsmax

  nsp = neutrmax + nionsmax + 1

  ALLOCATE(name(0:nsp),istat(nsp),mmw(nsp),nhyd(nsp),ncar(nsp),n14n(nsp),                           &
       n15n(nsp),noxy(nsp),dtype(nsp),ad(nsp),sd(nsp),phi(nsp),sd_2(nsp),                           &
       sd_3(nsp),ibnd(nsp,2),bval(nsp,2),ichrg(nsp))

  READ(61,"(A)") header
  DO nm = 1, neutrmax
     READ(61,961) ndum, name(nm), istat(nm), ichrg(nm), mmw(nm), nhyd(nm),                          &
          ncar(nm), n14n(nm), n15n(nm), noxy(nm), dtype(nm), ad(nm), sd(nm),                        &
          phi(nm), sd_2(nm),sd_3(nm), ibnd(nm,1), bval(nm,1), ibnd(nm,2), bval(nm,2)
  END DO

961 FORMAT(I4,1X,A12,1X,I1,1X,I2,1X,F7.1,5I4,2X,I4,1X,ES9.2,1X,F6.3,1X,ES9.2,1X,                    &
         F5.1,1X,F6.1,6X,I1,2X,ES10.3,3X,I1,2X,ES10.3)

  READ(62,"(A)") header
  ichk = 0
  DO n = 1, nionsmax
     nm = neutrmax + n
     READ(62,962) ndum, name(nm), istat(nm), ichrg(nm), mmw(nm), nhyd(nm),                          &
          ncar(nm), n14n(nm),  n15n(nm), noxy(nm), ibnd(nm,1), bval(nm,1),                          &
          ibnd(nm,2), bval(nm,2)
    IF(istat(nm) == 1) ichk = ichk + 1 
  END DO
962 FORMAT(I4,1X,A12,1X,I1,1X,I2,1X,F7.1,5I4,2(I3,ES11.3))

  CLOSE(unit=61)
  CLOSE(unit=62)
     
  name(0)   = '            '
  name(nsp) = 'E           '
  mmw(nsp) = zero

  nx1 = 0; nx2 = 0
  DO nm = 1, nsp-1
     IF(istat(nm) == 1) nx1=nx1+1
     IF(istat(nm) == 2) nx2=nx2+1
  END DO
  nchem = nx1
  ndiff = nx2

  IF(ichk > 0) THEN
     lions = .true.
     neq = nchem + 1
  ELSE
     lions = .false.
     neq = nchem
  END IF

  ALLOCATE(locp(neq),lopc(0:nsp),ldcp(ndiff),ldpc(0:nsp))
  locp(:) = 0; lopc(:) = 0; ldcp(:) = 0; ldpc(:) = 0

  nx = 0
  DO nm = 1, nsp-1
     IF(istat(nm) == 1) THEN
        nx=nx+1
        locp(nx) = nm
        lopc(nm) = nx
     END IF
  END DO
  
  IF(lions) THEN
     locp(neq) = nsp
     lopc(nsp) = neq
  END IF

  WRITE(*,"(' CHEMICAL SPECIES')")
  IF(SIZE(locp) > 0) THEN
     WRITE(*,"(10(2X,A12,1X))") (name(locp(nx)),nx=1,nchem)
  ELSE
     WRITE(*,"('  NONE')")
  END IF

  nx = 0
  DO nm = 1, nsp-1
     IF(istat(nm) == 2) THEN
        nx=nx+1
        ldcp(nx) = nm
        ldpc(nm) = nx
     END IF
  END DO

  WRITE(*,"(' DIFFUSING SPECIES')")
  IF(SIZE(ldcp) > 0) THEN
     WRITE(*,"(10(2X,A12,1X))") (name(ldcp(nx)),nx=1,ndiff)
  ELSE
     WRITE(*,"('  NONE')")
  END IF

  iN2 =  FIND_NAME('N2          ',name)
  iCH4 = FIND_NAME('CH4         ',name)
  iELE = nsp

RETURN
END SUBROUTINE READ_MOLECULES
