SUBROUTINE READ_PHOTOB(name,nbrnchB,loabB,loprB,ionizeB,enrgIB,charge_stateB,phrctB,wcrsB,delwB,xcrsB,bratB)
     
  USE PRECISION
  USE CONSTANTS
  USE SUBS, ONLY : FIND_NAME, INTRP

  IMPLICIT NONE

  !
  !  .. External Variables
  !

  CHARACTER(len=*), INTENT(IN), DIMENSION(0:) :: name
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:):: nbrnchB
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:):: loabB
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:):: loprB
  LOGICAL, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: ionizeB
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: enrgIB
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: charge_stateB
  CHARACTER(len=*), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: phrctB
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: wcrsB
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: delwB
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: xcrsB
  REAL(RP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:) :: bratB

  !
  !  .. Internal Variables
  !

  INTEGER, PARAMETER :: ncrsB = 540000 ! wavelengths
  INTEGER, PARAMETER :: nabsB = 3
  INTEGER, PARAMETER :: nbrmaxB = 12
  INTEGER, PARAMETER :: nprmaxB =  4

  CHARACTER(len=128) :: header, cline
  CHARACTER(len=12), DIMENSION(6) :: fm
  INTEGER :: nwav_co2, nbrnch_co2, nsp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wav_co2, crs_co2_tot, crs_co2_ion
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: brat_co2
  REAL(RP) :: rdum
  INTEGER :: na, nf, nw, nb, nm, np, j
!  CHARACTER(len=1) :: iret


  !
  !  .. Allocate External Arrays
  !

  ALLOCATE(nbrnchB(nabsB),loabB(nabsB),loprB(nprmaxB,nbrmaxB,nabsB),phrctB(nbrmaxB,nabsB),          &
       ionizeB(nbrmaxB,nabsB),charge_stateB(nbrmaxB,nabsB),wcrsB(ncrsB),enrgIB(nbrmaxB,nabsB),      &
       delwB(ncrsB),xcrsB(ncrsB,nabsB),bratB(ncrsB,nbrmaxB,nabsB))

  loprB(:,:,:) = 0
  loabB(:) = 0
  enrgIB(:,:) = 0
  charge_stateB(:,:) = 0

  !
  !  .. N2 High Resolution Lewis et al. Cross Section for 800-1000 Angstroms
  !

  na = 1
  OPEN(Unit=64,file='../data/photons/A.14.e18h6.getlines.output.150K;type=i',status='old',action='read')
     READ(64,"(A)") header
     READ(64,"(A)") header
     DO nf = 1, ncrsB
        nw = ncrsB-nf+1
        READ(64,*) wcrsB(nw), xcrsB(nw,na) ! wavelengths, cross sections
        wcrsB(nw) = 1.E8_RP/wcrsB(nw) ! wavelengths
     END DO
  CLOSE(unit=64)

  loabB(na) = FIND_NAME('N2          ',name)
  nbrnchB(na) = 1
  nb = 1
  phrctB(nb,na) = 'N2           + hv           = N2D          + N            +              +             '
  ionizeB(nb,na) = .false.
  charge_stateB(nb,na) = zero
  enrgIB(nb,na) = zero
  bratB(1:ncrsB,nb,na) = one
  loprB(1,nb,na) = FIND_NAME('N2D         ',name)
  loprB(2,nb,na) = FIND_NAME('N           ',name)
  loprB(3:4,nb,na) = 0

  !
  !  .. 29N2 High Resolution Cross Section in Band Region
  !

  na = 2
  OPEN(Unit=65,file='../data/photons/A.1415.e18h6.getlines.output.150K;type=i',status='old',action='read')
     READ(65,"(A)") header
     READ(65,"(A)") header
     DO nf = 1, ncrsB
        nw = ncrsB-nf+1
        READ(65,*) wcrsB(nw), xcrsB(nw,na)
        wcrsB(nw) = 1.E8_RP/wcrsB(nw)
     END DO
  CLOSE(unit=65)


  loabB(na) = FIND_NAME('N2I         ',name)
  nbrnchB(na) = 2

  nb = 1
  phrctB(nb,na) = 'N2I          + hv           = N2DI         + N            +              +             '
  ionizeB(nb,na) = .false.
  charge_stateB(nb,na) = zero
  enrgIB(nb,na) = zero
  bratB(1:ncrsB,nb,na) = half
  loprB(1,nb,na) = FIND_NAME('N2DI        ',name)
  loprB(2,nb,na) = FIND_NAME('N           ',name)
  loprB(3:4,nb,na) = 0

  nb = 2
  phrctB(nb,na) = 'N2I          + hv           = N2D          + NI           +              +             '
  ionizeB(nb,na) = .false.
  charge_stateB(nb,na) = zero
  enrgIB(nb,na) = zero
  bratB(1:ncrsB,2,na) = half
  loprB(1,nb,na) = FIND_NAME('N2D         ',name)
  loprB(2,nb,na) = FIND_NAME('NI          ',name)
  loprB(3:4,nb,na) = 0

  !
  !  .. Wavelength intervals for cross section grid
  !

  DO nf = 2, ncrsB-1                               ! Need this to integrate over wavelength
     delwB(nf) = half*(wcrsB(nf+1)-wcrsB(nf-1))
  END DO
  delwB(1) = delwB(2)
  delwB(ncrsB) = delwB(ncrsB-1)
  
  !
  !  .. Read CO2 cross section and interpolate to high res grid
  !

  nsp = SIZE(name)-1

  na = 3
  loabB(na) = FIND_NAME('CO2         ',name)
  OPEN(unit=66,file='../data/photons/PHOTO_CO2.DAT',status='unknown',action='read')
     READ(66,*) nwav_co2, nbrnch_co2 ! wavelengths, reactions
     ALLOCATE(wav_co2(nwav_co2),crs_co2_tot(nwav_co2),crs_co2_ion(nwav_co2),brat_co2(nwav_co2,nbrnch_co2))
     READ(66,'(A)') header
     READ(66,*) wav_co2 ! wavelength scale
     READ(66,'(A)') header
     READ(66,*) crs_co2_tot ! total absorption cross section
     READ(66,'(A)') header
     READ(66,*) crs_co2_ion ! total ionization cross section
     DO nb = 1, nbrnch_co2
        READ(66,'(A)') cline
        phrctB(nb,na) = cline
        READ (cline,"(5(A12,3X),A12,F6.1)") (fm(j),j=1,6), rdum
        enrgIB(nb,na) = rdum
        np = 0
        charge_stateB(nb,na) = zero
        loprB(1:4,nb,na) = 0 
        DO j = 3, 6
           nm = FIND_NAME(fm(j),name)
           IF ((nm > 0) .and. (nm <= nsp)) THEN
              np = np + 1
              loprB(np,nb,na) = nm
           END IF
           IF (nm == nsp) THEN 
              charge_stateB(nb,na) = charge_stateB(nb,na) + one
              ionizeB(nb,na) = .true.
           END IF
        END DO
        READ(66,'(10ES11.3)')  brat_co2(1:nwav_co2,nb)
     END DO
  CLOSE(unit=66)
  nbrnchB(na) = nbrnch_co2


  CALL INTRP(wav_co2,crs_co2_tot(1:nwav_co2),wcrsB(1:ncrsB),xcrsB(1:ncrsB,na))
  DO nb = 1, nbrnchB(na)
     CALL INTRP(wav_co2,brat_co2(1:nwav_co2,nb),wcrsB(1:ncrsB),bratB(1:ncrsB,nb,na))
  END DO

  DEALLOCATE(wav_co2,crs_co2_tot,crs_co2_ion,brat_co2)

  RETURN
END SUBROUTINE READ_PHOTOB
