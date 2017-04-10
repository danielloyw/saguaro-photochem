SUBROUTINE READ_ATMOS

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_NAME, LOCATE

  IMPLICIT NONE

  !
  !  .. Internal Variables
  !

  REAL(RP), PARAMETER :: DMIN = 1.0E-50_RP
  INTEGER :: nx, nm, nz, nl, nmol
  CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: namex
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indx_den
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: found
  CHARACTER(len=128) :: header
  CHARACTER(len=12) :: xname
  REAL(RP) :: edz

  !
  !  .. Main Section
  !

  ALLOCATE(found(nsp))
  found(1:nsp) = .false.
  OPEN (unit=63,file='atm1D.in',status='old',action='read')
  
     READ (63, *) nlev, nmol
     READ (63,'(A)') header

     !  .. Allocate Arrays

     ALLOCATE(namex(nmol),indx_den(nmol),z(nlev),rz(nlev),grv(nlev),ek(nlev),tn(nlev),ti(nlev),     &
          te(nlev),prs(nlev),mass(nlev),rho(nlev),den(nlev,0:nsp),den_old(nlev,0:nsp),xmol(nlev,nsp))
     
     !  .. Read Molecules names, associate with reactions.dat

     READ(63,"(10(2X,A12,1X))") (namex(nx),nx=1,nmol)
     DO nx = 1, nmol
        nm = FIND_NAME(namex(nx),name)
        indx_den(nx) = nm
        IF(nm /= 0) THEN
           found(nm) = .true.
        END IF
     END DO

     !  .. Read Atmosphere Arrays

     READ (63,'(A)') header
     READ (63,*) (z(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (rz(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (grv(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (tn(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (te(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (prs(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (rho(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (mass(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (ek(nz), nz=1, nlev)
     READ (63,'(A)') header
     READ (63,*) (den(nz,0), nz=1, nlev) ! den(*,0) is for overall density
     DO nx = 1, nmol
        nm = indx_den(nx)
        READ (63,'(A)') xname
        READ (63,*) (den(nz,nm), nz=1, nlev)
     END DO
        
  CLOSE(unit=63)

  ! ####################################################
  !   scale CO density to 50 ppm
  ! ####################################################
!!!
!  nm = FIND_NAME('H2           ',name)
!  den(1:nlev,nm) = DMIN

!  nm = FIND_NAME('H            ',name)
!  den(1:nlev,nm) = DMIN
 

!  DO nm = 1, neutrmax
!     IF(.not. found(nm)) THEN
!        den(1:nlev,nm) = DMIN
!        WRITE(*,"(' Assign density of ',ES11.3,' to ',A12)") DMIN, name(nm)
!     END IF
!  END DO
!!!

  DO nx = 1, ndiff
     nm = ldcp(nx)
     IF(ibnd(nm,1) == 3) THEN
        den(1,nm) = bval(nm,1)
     END IF
  END DO

  !  .. Set ion density = min val if species not found
      
  DO nm = neutrmax+1, nsp-1
     IF(.not. found(nm)) THEN
        den(1:nlev,nm) = DMIN
        WRITE(*,"(' Assign density of ',ES11.3,' to ',A12)") DMIN, name(nm)
     END IF
  END DO

  !  .. Set electron density = sum of ion density

  DO nz = 1, nlev
     den(nz,nsp) = SUM(den(nz,neutrmax+1:nsp-1))
  END DO

 !  ..  Set neutral density = min val if below minimum

  WHERE (den < DMIN) den=DMIN

  DO nz = 1, nlev
     den(nz,0) = SUM(den(nz,1:nsp))
  END DO

  DO nm = 1, nsp-1
     IF(istat(nm) == 3) THEN
        WRITE(*,"(' Assigning ',A12,' density')") name(nm)
        den(1:nlev,nm) = bval(nm,1)*den(1:nlev,0)
     END IF
  END DO

  !  .. Calculate mole fraction

  DO nz = 1, nlev
     DO nm = 1, nsp
        xmol(nz,nm) = den(nz,nm)/den(nz,0)
     END DO
  END DO


  !  .. Reset pressure to be consistent with N & T

  DO nz = 1, nlev
     prs(nz) = rkb*tn(nz)*den(nz,0)
  END DO

  !  .. Calculate mass density, mean mass, and gas constant
   
  DO nz = 1, nlev
     IF (nionsmax > 0) THEN
        rho(nz) = amu*DOT_PRODUCT(mmw(1:nsp),den(nz,1:nsp))
        mass(nz) = rho(nz)/(den(nz,0)+two*den(nz,nsp))/amu ! ions have 2x the scale height
     ELSE
        rho(nz) = amu*DOT_PRODUCT(mmw(1:nsp),den(nz,1:nsp))
        mass(nz) = rho(nz)/den(nz,0)/amu
     END IF
  END DO  

  CALL HYDROST

  nibot = LOCATE(z,zbot)

  ! ############################################################################
  ! #                                                                          #
  ! #          AUXILIARY ALTITUDE ARRAYS                                       #
  ! #                                                                          #
  ! ############################################################################


  ALLOCATE(rmid(nlev),drp(nlev),dr(nlev),rp2(nlev),rm2(nlev),ht(nlev,0:nsp),   &
       col(nlev))

  !  .. Define some altitude-related arrays

  DO nl = 1, nlev-1
     rmid(nl) = half*(rz(nl+1)+rz(nl))
     drp(nl) = rz(nl+1)-rz(nl)
     rp2(nl) = (rmid(nl)/rz(nl))**2
  END DO

  DO nl = 2, nlev
     dr(nl) = rmid(nl)-rmid(nl-1)
     rm2(nl) = (rmid(nl-1)/rz(nl))**2
  END DO

  !  .. Calculate scale heights

  DO nm = 1, nsp-1                                    ! neutrals
     DO nl = 1, nlev
        ht(nl,nm) = rkb*tn(nl)/grv(nl)/mmw(nm)/amu
     END DO
  END DO

  DO nl = 1, nlev                                    ! mean scale height
     ht(nl,0) = rkb*tn(nl)/grv(nl)/mass(nl)/amu
  END  DO

  !  .. Calculate column density (vertically integrated)

  col(nlev) = den(nlev,0)*ht(nlev,0)
  DO nl = nlev-1, 1, -1
        col(nl)=col(nl+1)+half*(den(nl,0)+den(nl+1,0))*drp(nl)
  END DO

  ti(1:nlev) = 230._RP

  !  .. Calculate eddy profile

!  DO nz = 1, nlev
!     edz = ed1*(pk0/prs(nz))**gek
!     ek(nz) = ed0*edz/(ed0+edz)
!     IF(ek(nz) < ed1) ek(nz) = ed1
!     IF(ek(nz) > ed0) ek(nz) = ed0
!  END DO

   nibot = LOCATE(z,zbot)
   nzion = nlev-nibot+1

  DEALLOCATE(namex,found,indx_den)
  RETURN
END SUBROUTINE READ_ATMOS
