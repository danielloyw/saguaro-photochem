SUBROUTINE PHOTO 

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS_PHOTO
  USE SUBS, ONLY : INTRP

  IMPLICIT  none

  !
  !  .. Internal Variables
  !

  CHARACTER(len=128) :: file_photoA, file_photoC, file_jvals, file_sol_ionz, file_sol_uv,           &
       file_sol_bnd, file_sol_hres
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: delwB
  CHARACTER(len=87), ALLOCATABLE, DIMENSION(:,:) :: phrctA, phrctB, phrctC, phrctJ
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: loprA, loprB, loprC, loprJ
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wav_sol_ionz, flx_sol_ionz,wav_sol_uv,flx_sol_uv,          &
       wav_sol_hres, flx_sol_hres, wav_sol_bnd,flx_sol_bnd
  INTEGER :: na, nb, nw, nph, np

  !
  !  .. Initialization
  !
    
  file_photoA = '../data/photons/photoA.dat'
  file_photoC = '../data/photons/photoC.dat'
  file_jvals = '../data/photons/jvals.dat'
  file_sol_ionz = '../data/solar/sorce_ionz_T40.dat'
  file_sol_bnd = '../data/solar/sorce_bnd_T40.dat'
  file_sol_uv = '../data/solar/sorce_uv_T40.dat'
  file_sol_hres = '../data/solar/sumer.new'

  !
  !  .. Read photolysis cross sections and branching ratios for 3 wavelength regions
  !

  !
  !  .. Region A: 0-800 Angstroms
  !
 
  CALL READ_PHOTO(name,file_photoA,nbrnchA,loabA,loprA,ionizeA,enrgIA,charge_stateA,phrctA,         &   !!!arguments
    wcrsA,xcrsA,bratA)

  ncrsA = SIZE(wcrsA)
  nabsA = SIZE(xcrsA,2)
  nbrmaxA = SIZE(bratA,2)

  !
  !  .. Region B: 800-1000 Angstroms 
  !
 
  CALL READ_PHOTOB(name,nbrnchB,loabB,loprB,ionizeB,enrgIB,charge_stateB,phrctB,wcrsB,delwB,xcrsB,  &
       bratB)

  ncrsB = SIZE(wcrsB)
  nabsB = SIZE(xcrsB,2)
  nbrmaxB = SIZE(bratB,2)

  ! 
  !  .. Region C: 1000-3000 Angstroms
  !
 
  CALL READ_PHOTO(name,file_photoC,nbrnchC,loabC,loprC,ionizeC,enrgIC,charge_stateC,phrctC,wcrsC,   & 
    xcrsC,bratC)

  ncrsC = size(wcrsC)
  nabsC = size(xcrsC,2)
  nbrmaxC = SIZE(bratC,2)
  !
  !  .. Optically thin absorbers
  !

  CALL READ_JVALS(name,file_jvals,nbrnchJ,loabJ,loprJ,phrctJ,srateJ)                               

  nabsJ = SIZE(srateJ,2)
  nbrmaxJ = SIZE(phrctJ,1)

  ALLOCATE(fsolA(ncrsA),fsolB(ncrsB),fsolC(ncrsC))

  !
  !  .. Read solar spectra wav < 800 A
  !

  CALL READ_SOL(file_sol_ionz,wav_sol_ionz,flx_sol_ionz)
  
  DO nw = 1, ncrsA
     fsolA(nw) = flx_sol_ionz(nw)
  END DO

  ! .. Band region 800 A < wav < 1000 A

  CALL READ_SOL(file_sol_bnd,wav_sol_bnd,flx_sol_bnd)
  flx_sol_bnd = 0.1_RP*flx_sol_bnd   ! convert from per bin to per A

  CALL READ_SOL_HRES(file_sol_hres,wav_sol_hres,flx_sol_hres)

  !
  !  .. Normalize High Res Spectrum to Low Res Spectrum
  !
  
  CALL NORM_SOL(wav_sol_bnd,flx_sol_bnd,wav_sol_hres,flx_sol_hres)

  !
  !  .. Interpolate High Res Spectrum to Low Res Grid
  !
  
  CALL INTRP(wav_sol_hres,flx_sol_hres,wcrsB,fsolB)

  !
  !  .. Convert High Res Spectrum from ph/cm2/s/A to ph/cm2/s/bin
  !

  DO nw = 1, ncrsB
     fsolB(nw) = delwB(nw)*fsolB(nw)
  END DO

  !  .. wav > 1000

  CALL READ_SOL(file_sol_uv,wav_sol_uv,flx_sol_uv)

  DO nw = 1, ncrsC
     fsolC(nw) = flx_sol_uv(nw)
  END DO
  
  !
  !  .. Scale to orbital distance
  !

  fsolA = fsolA/DH**2                  
  fsolB = fsolB/DH**2                  
  fsolC = fsolC/DH**2                  
  srateJ = srateJ/DH**2

  !
  !  .. Check solar fluxes
  !

!  WRITE(*,"(' Solar Check: ',3ES11.3)") (DH**2)*SUM(fsolA(1:ncrsA)),(DH**2)*SUM(fsolB(1:ncrsB)),   &
!       (DH**2)*SUM(fsolC(1:ncrsC))


  ! ################################################################################################
  ! #                                                                                              #
  ! #                                    Put Regions Together                                      #
  ! #                                                                                              #
  !#################################################################################################

  nph = 0
  DO na = 1, nabsA
     DO nb = 1, nbrnchA(na)
        nph = nph + 1
     END DO
  END DO
  DO na = 1, nabsB
     DO nb = 1, nbrnchB(na)
        nph = nph + 1
     END DO
  END DO
  DO na = 1, nabsC
     DO nb = 1, nbrnchC(na)
        nph = nph + 1
     END DO
  END DO
  DO na = 1, nabsJ
     DO nb = 1, nbrnchJ(na)
        nph = nph + 1
     END DO
  END DO
  nphrt = nph

  ALLOCATE(irpt(5,nphrt),ptitle(nphrt)) 
  nph = 0
  DO na = 1, nabsA
     DO nb = 1, nbrnchA(na)
        nph = nph + 1
        irpt(1,nph) = loabA(na)
        irpt(2:5,nph) = loprA(1:4,nb,na)
        ptitle(nph) = phrctA(nb,na)
     END DO
  END DO

  DO na = 1, nabsB
     DO nb = 1, nbrnchB(na)
        nph = nph + 1
        irpt(1,nph) = loabB(na)
        irpt(2:5,nph) = loprB(1:4,nb,na)
        ptitle(nph) = phrctB(nb,na)
     END DO
  END DO

  DO na = 1, nabsC
     DO nb = 1, nbrnchC(na)
        nph = nph + 1
        irpt(1,nph) = loabC(na)
        irpt(2:5,nph) = loprC(1:4,nb,na)
        ptitle(nph) = phrctC(nb,na)
     END DO
  END DO

  DO na = 1, nabsJ
     DO nb = 1, nbrnchJ(na)
        nph = nph + 1
        irpt(1,nph) = loabJ(na)
        irpt(2:5,nph) = loprJ(1:4,nb,na)
        ptitle(nph) = phrctJ(nb,na)
     END DO
  END DO

  ALLOCATE(trnA(ncrsA,nlev),trnB(ncrsB,nlev),trnC(ncrsC,nlev), & 
     prtA(ncrsA,nbrmaxA,nabsA),prtB(ncrsB,nbrmaxB,nabsB),prtC(ncrsC,nbrmaxC,nabsC))
  
  ! #################################################################################################

  DEALLOCATE(loprA, loprB, loprC, loprJ, phrctA, phrctB, phrctC, phrctJ, delwB, wav_sol_ionz,       &
       flx_sol_ionz, wav_sol_bnd, flx_sol_bnd, wav_sol_uv, flx_sol_uv, wav_sol_hres, flx_sol_hres)

  ALLOCATE(rph(nphrt,nlev),rpt(nphrt,nlev))


  RETURN

END SUBROUTINE PHOTO
