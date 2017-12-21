SUBROUTINE READ_PHOTOB(name,nbrnchB,loabB,loprB,ionizeB,enrgIB,charge_stateB,phrctB,wcrsB,delwB,xcrsB,bratB)
     
  USE PRECISION
  USE CONSTANTS
  USE SUBS, ONLY : FIND_NAME, INTRP, LOCATE

  IMPLICIT NONE

  !
  !  .. External Variables
  !

  CHARACTER(len=*), INTENT(IN), DIMENSION(0:) :: name
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:):: nbrnchB
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:):: loabB
  INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:,:):: loprB ! product index (product #, branch #, species #)
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

  INTEGER, PARAMETER :: ncrsB = 1500001 ! wavelengths
  INTEGER, PARAMETER :: nwav_n2 = 540000
  INTEGER, PARAMETER :: nwav_coA = 16000
  INTEGER, PARAMETER :: nwav_coB = 840000

  CHARACTER(len=128) :: header, cline
  CHARACTER(len=12), DIMENSION(6) :: fm
  INTEGER :: nwav_co2, nbrnch_co2, nsp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wav_co, crs_co_diss, crs_co_tot, wav_n2, crs_n2
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: brat_co
  REAL(RP) :: rdum
  INTEGER :: na, nf, nw, nb, nm, np, j

  ! Read main cross-section file
  
  OPEN(unit=65, file='../data/photons/photoB.dat',status='old',action='read')

    READ(65,*) ncrs, nabs, nbrmaxB ! file wavelengths, species, maximum branches
    nabsB = nabs + 4
    ALLOCATE(nbrnchB(nabsB),loabB(nabsB),loprB(nprmaxB,nbrmaxB,nabsB),phrctB(nbrmaxB,nabsB),          &
       ionizeB(nbrmaxB,nabsB),charge_stateB(nbrmaxB,nabsB),wcrsB(ncrsB),enrgIB(nbrmaxB,nabsB),      &
       delwB(ncrsB),xcrsB(ncrsB,nabsB),bratB(ncrsB,nbrmaxB,nabsB))

    loprB(:,:,:) = 0
    loabB(:) = 0 ! index of reactant
    enrgIB(:,:) = zero
    charge_stateB(:,:) = 0
    ionizeB = .false.
    
    ! Create wavelength scale
  
    delwB(:) = 2.E-4_RP
    DO nw = 0, ncrsB-1
       wcrsB(nw) = 790._RP + nw * delwB(nw)
    END DO

    ALLOCATE(wcrs(ncrs),xcrs(ncrs),xdum(ncrs),brat(ncrs))
    
    READ(65,'(A)') header                               
    READ(65,'(10ES11.3)') wcrs ! wavelength scale
    DO na = 1, nabs
       READ(65,'(A12, I4)') xname, ndum
       loabB(na) = FIND_NAME(xname,name)
       nbrnchB(na) = ndum ! number of reactions/branches
       READ(65,'(A)') header
       READ(65,*) (xcrs(nw),nw=1,ncrs) !'(10ES11.3)' ! total absorption cross section
       READ(65,'(A)') header
       READ(65,*) (xdum(nw),nw=1,ncrs) !'(10ES11.3)' ! total ionization cross section
       DO nb = 1, nbrnchB(na)
          READ (65,'(A)') cline
          phrctB(nb,na) = cline
          READ (cline,"(5(A12,3X),A12,F7.1)") (fm(j),j=1,6), rdum ! reactants/products
          enrgIB(nb,na) = rdum
          charge_stateB(nb,na) = zero
          loprB(1:4,nb,na) = 0
          np = 0
          DO j = 3, 6 ! products
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
          READ(65,'(10ES11.3)') (brat(nw),nw=1,ncrs) ! branching ratios
          CALL INTRP(wcrs,brat,wcrsB,bratB(1:ncrsB,nb,na))
       END DO
       CALL INTRP(wcrs,xcrs,wcrsB,xcrsB(1:ncrsB,na))
    END DO
     
  CLOSE(unit=65)
  DEALLOCATE(wcrs,xcrs,xdum,brat)
  
  !
  !  .. N2 High Resolution Lewis et al. Cross Section for 800-1000 Angstroms
  !
  ALLOCATE(wav_n2(nwav_n2),crs_n2(nwav_n2))
  na = nabs + 1
  OPEN(Unit=64,file='../data/photons/photoB-28N2.dat',status='old',action='read')
     READ(64,"(A)") header
     READ(64,"(A)") header
     DO nf = 1, nwav_n2
        nw = nwav_n2-nf+1
        READ(64,*) wav_n2(nw), crs_n2(nw) ! wavelengths, cross sections
        wav_n2(nw) = 1.E8_RP/wav_n2(nw) ! convert cm-1 to angstrom
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

  CALL INTRP(wav_n2,crs_n2,wcrsB,xcrsB(1:ncrsB,na))

  !
  !  .. 29N2 High Resolution Cross Section in Band Region
  !

  na = nabs + 2
  OPEN(Unit=65,file='../data/photons/photoB-29N2.dat',status='old',action='read')
     READ(65,"(A)") header
     READ(65,"(A)") header
     DO nf = 1, nwav_n2
        nw = nwav_n2-nf+1
        READ(65,*) wav_n2(nw), crs_n2(nw)
        wav_n2(nw) = 1.E8_RP/wav_n2(nw)
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

  CALL INTRP(wav_n2,crs_n2,wcrsB,xcrsB(1:ncrsB,na))
  DEALLOCATE(wav_n2, crs_n2)

  !
  !  .. Read CO cross section and interpolate to high res grid
  !

  na = nabs + 3
  nbrnchB(na) = 2
  ALLOCATE(wav_co(nwav_coA+nwav_coB),crs_co_diss(nwav_coA+nwav_coB),crs_co_tot(nwav_coA+nwav_coB),brat_co(nwav_coA+nwav_coB,nbrnchB(na)))
  OPEN(Unit=67,file='../data/photons/photoB-12C16O_300K_89.6-91.2.dat',status='old',action='read')
     READ(67,"(A)") header
     READ(67,"(A)") header
     DO nf = 1, nwav_coA
        READ(67,*) wav_co(nf), crs_co_diss(nf) ! wavelengths, cross sections
        wav_co(nf) = wav_co(nf)*10._RP
        crs_co_tot(nf) = crs_co_diss(nf)
     END DO
  CLOSE(unit=67)
  OPEN(Unit=68,file='../data/photons/photoB-12C16O_300K_91.2-108.dat',status='old',action='read')
     READ(68,"(A)") header
     READ(68,"(A)") header
     READ(68,"(A)") header
     READ(68,"(A)") header
     DO nf = nwav_coA+1, nwav_coA+nwav_coB
        READ(68,*) wav_co(nf), crs_co_tot(nf), crs_co_diss(nf) ! wavelengths, total cross section, dissociative cross section
        wav_co(nf) = wav_co(nf)*10._RP
     END DO
  CLOSE(unit=68)

  loabB(na) = FIND_NAME('CO          ',name)

  nb = 1
  phrctB(nb,na) = 'CO           + hv           = C            + O            +              +             '
  ionizeB(nb,na) = .false.
  charge_stateB(nb,na) = zero
  enrgIB(nb,na) = zero
  DO nf = 1, nwav_coA+nwav_coB
     brat_co(nf,nb) = crs_co_diss(nf)/crs_co_tot(nf)
  END DO
  loprB(1,nb,na) = FIND_NAME('C           ',name)
  loprB(2,nb,na) = FIND_NAME('O           ',name)
  loprB(3:4,nb,na) = 0

  nb = 2
  phrctB(nb,na) = 'CO           + hv           = COP          + E            +              +             '
  ionizeB(nb,na) = .true.
  charge_stateB(nb,na) = one
  enrgIB(nb,na) = 885.6_RP
  DO nf = 1, nwav_coA+nwav_coB
     brat_co(nf,nb) = one - brat_co(nf,1)
  END DO
  loprB(1,nb,na) = FIND_NAME('COP         ',name)
  loprB(2,nb,na) = FIND_NAME('E           ',name)
  loprB(3:4,nb,na) = 0

  CALL INTRP(wav_co,crs_co_tot,wcrsB,xcrsB(1:ncrsB,na))
  DO nb = 1, nbrnchB(na)
     CALL INTRP(wav_co,brat_co(:,nb),wcrsB,bratB(1:ncrsB,nb,na))
  END DO

  
  DEALLOCATE(wav_co, crs_co_diss, crs_co_tot, brat_co)

  RETURN
END SUBROUTINE READ_PHOTOB
