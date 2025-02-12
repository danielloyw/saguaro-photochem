PROGRAM PCHEM1D

  !  .. Modules

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : LOCATE, INTRP, FIND_NAME

  !  .. Local Variable Declarations

  IMPLICIT none

  REAL(RP) :: tol_diff, tol_chem, test_chem, test_diff, test_time
  INTEGER :: ntim_diff, ntim_chem, iter_chem, iter_diff ! max diffusion time steps (overall number of loops), max reaction time steps per diffusion step, iteration count for reactions, iteration count for diffusion
  INTEGER :: nm, nl, nx, np, nln, nmn, nm1, nm2, nm3, nm4, nm5
  INTEGER :: isol, iprnt
  REAL(RP) :: dnm, sm, tst1

  ! #################################################################################################
  ! #                                                                                               #
  ! #                                               INPUT                                           #
  ! #                                                                                               #
  ! #################################################################################################

  lcompo_sol = .false.
  lchemeq_sol = .false.

  tol_diff = 1.E-10_RP
  tol_chem = 1.E-10_RP

  lprnt = .true.
  zbot = 0.E5_RP
  RSHADOW = RPLANET + zbot

!  WRITE(*,"(' Enter run ID ..........................  > ')",ADVANCE='NO')
!  READ(*,"(A)") runID

  OPEN(unit=20,file='runid.out',status='old',action='read')
     READ(20,"(A)") runid
  CLOSE(unit=20)

!  WRITE(*,"(' Enter time step and # of steps ........  > ')", ADVANCE='NO')
!  READ(*,*) tstep_diff, ntim_diff

  OPEN(unit=21,file='tctl.out',status='old',action='read')
     READ(21,*) tstep_diff, ntim_diff, tstep_chem, ntim_chem
  CLOSE(unit=21)


  !  .. Read in parameters for photolysis averaging

  OPEN(unit=60,file='pchem1D.sol',status='old',action='read')
     READ(60,*) cos_sza
     READ(60,*) diurnal_average
  CLOSE(unit=60)

  !   .. Read in parameters for run

  OPEN(unit=60,file='pchem1D.ctl',status='old',action='read')
     READ(60,*) ed1,ed0,pk0,gek !3e7, 3e7, 1e5, 0
     READ(60,*) ascl ! 1
     READ(60,*) ihetero ! 0
     READ(60,*) next ! 0
     IF(next > 0) THEN
     ALLOCATE(name_ext(next),fext(next),zext(next),hext(next))
     DO nm = 1, next
        READ(60,'(A12,3ES11.3)') name_ext(nm), fext(nm), zext(nm), hext(nm)
     END DO
     END IF
     READ(60,*) iaer ! 0
     READ(60,*) iprnt, isol ! frequency to write results to disc, frequency to recalculate photolysis and photoelectrons
  CLOSE(unit=60)


  ! #################################################################################################
  ! #                                                                                               #
  ! #                                     PREPARE CALCULATIONS                                      #
  ! #                                                                                               #
  ! #################################################################################################

  CALL READ_MOLECULES           !  Read in molecule data

  CALL READ_REACTIONS           !  Read in reaction rate data

  CALL READ_ATMOS               !  Read model atmosphere

  CALL DIFCO                    !  Calculate diffusion coefficients

  CALL PHOTO                    !  Set up photolyis calculation

  CALL ELCTRN                   !  Set up suprathermal electron calculation

  CALL READ_AEROSOL             !  Read Aerosol Model

  CALL READ_RAYLEIGH            !  Calculate Rayleigh optical depth

  CALL PATHS1D                  !  Calculate path lengths for solar absorption

  CALL RATECO                   !  Calculate reaction rate coeff

  CALL EXT_PROD                 !  Specified production rates

  ALLOCATE(pr_ph(nlev,0:nsp),ls_ph(nlev,0:nsp),     & !  Production and Loss from photon processes
           pr_pe(nlev,0:nsp),ls_pe(nlev,0:nsp),     & !  Production and Loss from electron processes
           pr_chem(nlev,0:nsp),ls_chem(nlev,0:nsp), & !  Production and Loss from chemical processes
           pr(nlev,0:nsp),ls(nlev,0:nsp),           & !  Net production and loss
           div_flx(nlev,0:nsp),dNdt(nlev,0:nsp),    & !  Flux divergence and net balance
           rcdn(nlev,0:nsp),flx(nlev,0:nsp))          !  Condensation loss and flux


  ! #################################################################################################
  ! #                                                                                               #
  ! #                                        TIME STEP LOOP                                         #
  ! #                                                                                               #
  ! #################################################################################################


  WRITE(*,"(' START TIME LOOP')")
  iter_diff = 0; test_diff = 10._RP*tol_diff
  TIME: DO

     iter_diff = iter_diff + 1
     IF(iter_diff > ntim_diff) EXIT

     den_old(1:nlev,0:nsp) = den(1:nlev,0:nsp)

     IF (MOD(iter_diff-1,isol) == 0) THEN

        !  .. Photolysis Reactions

     !  Occasionally calculate photolysis rate using high resolution in region B (lcrsB = .true.)
        lcrsA = .true.; lcrsB = .true.; lcrsC = .true.; lcrsJ = .true.

        CALL SOLDEP1

        !  .. Suprathermal Electron Reactions

        CALL ELDEP1

     END IF

      ! ---------------------------------------------------------------------------------------------
     !
     !                                 CHEMICAL EQUILIBRIUM SPECIES
     !
     ! ---------------------------------------------------------------------------------------------

     lcrsA = .true.; lcrsB = .false.; lcrsC = .true.; lcrsJ = .true.

     iter_chem = 0; test_chem = 10._RP*tol_chem
     CHEM: DO

        iter_chem = iter_chem + 1
        IF((iter_chem > ntim_chem) .or. (test_chem < tol_chem)) EXIT

        den_old = den

        IF(lchemeq_sol) THEN
           CALL SOLDEP1
!           CALL ELDEP1
        END IF

        CALL CHEMEQ

        !  .. Check Convergence

        test_chem = zero; nmn = 1; nln = 1
        DO nx = 1, nchem
           nm = locp(nx)
           DO nl = 1, nlev
              IF(den(nl,nm) > 1.0E+00_RP) THEN
                 dnm = MAX(pr(nl,nm),ls(nl,nm))
                 IF(dnm > 1.E-6_RP) THEN
                    tst1 = ABS(pr(nl,nm)-ls(nl,nm)-dNdt(nl,nm))/dnm
                 ELSE
                    tst1 = zero
                 END IF
                 IF(tst1 > test_chem) THEN
                    test_chem = tst1
                    nln = nl
                    nmn = nm
                 END IF
              END IF
           END DO
        END DO

        IF (lprnt .and. (nchem > 0)) THEN
           WRITE(*,"(10X,'  CHEMEQ: ITER = ',I4,2X,A12,' Z = ',F8.1,' N = ',ES11.3,' P = ',ES11.3,' L = ',ES11.3,' dNdt = ',ES11.3,' Err = ',ES11.3,' Relerr = ',ES11.3)") &
                iter_chem,name(nmn),1.E-5_RP*z(nln),den(nln,nmn),pr(nln,nmn), ls(nln,nmn), dNdt(nln,nmn),pr(nln,nmn)-ls(nln,nmn)-dNdt(nln,nmn),test_chem
        END IF

     END DO CHEM

     ! ---------------------------------------------------------------------------------------------
     !
     !                                   DIFFUSING SPECIES
     !
     ! ---------------------------------------------------------------------------------------------

     lcrsA = .true.; lcrsB = .false.; lcrsC = .true.; lcrsJ = .true.
     CALL COMPO

!     CALL HYDROST
     ! .. Update stuff

     DO nl = 1, nlev
        den(nl,0) = sum(den(nl,1:nsp))
        DO nm = 1, nsp
           xmol(nl,nm) = den(nl,nm)/den(nl,0)
        END DO
     END DO

     DO nl = 1, nlev
        IF (nionsmax > 0) THEN
           rho(nl) = amu*DOT_PRODUCT(mmw(1:nsp),den(nl,1:nsp))
           mass(nl) = rho(nl)/(den(nl,0)+two*den(nl,nsp))/amu
        ELSE
           rho(nl) = amu*DOT_PRODUCT(mmw(1:nsp),den(nl,1:nsp))
           mass(nl) = rho(nl)/den(nl,0)/amu
        END IF
     END DO


     ! ---------------------------------------------------------------------------------------------
     !
     !                                   CHECK CONVERGENCE
     !
     ! ---------------------------------------------------------------------------------------------

     test_time = zero; nln = 1; nmn = 1

     IF(ndiff > 0) THEN
        DO nx = 1, ndiff
           nm = ldcp(nx)
           DO nl = 1, nlev
              IF(den(nl,nm) > 1.E0_RP) THEN
                 tst1 = ABS(den(nl,nm)-den_old(nl,nm))/den(nl,nm)
              ELSE
                 tst1 = zero
              END IF
              IF(tst1 > test_time) THEN
                 test_time = tst1
                 nln = nl
                 nmn = nm
              END IF
           END DO
        END DO
     END IF

     IF(nchem > 0) THEN
        DO nx = 1, nchem
           nm = locp(nx)
           DO nl = 1, nlev
              IF(den(nl,nm) > 1.E0_RP) THEN
                 tst1 = ABS(den(nl,nm)-den_old(nl,nm))/den(nl,nm)
              ELSE
                 tst1 = zero
              END IF
              IF(tst1 > test_time) THEN
                 test_time = tst1
                 nln = nl
                 nmn = nm
              END IF
           END DO
        END DO
     END IF

     WRITE(*,910) iter_diff,name(nmn),1.E-5_RP*z(nln),den(nln,nmn),den(nln,nmn)-den_old(nln,nmn),test_time

     !  .. Write to Disc

     IF (MOD(iter_diff-1,iprnt) == 0) THEN
        CALL COMPOUT
     END IF

  END DO TIME


  CALL COMPOUT                    !  .. Final output


910 FORMAT(10X,' STEP=',I6,' SPECIES = ',A12,' Z = ',ES12.4,' N = ',ES11.3,' DelN =',ES11.3,' Rerr = ',ES11.3)
911 FORMAT(5X,' STEP=',I6,2X,A12,' DNDT = ',ES11.3,' Pr = ',ES11.3,' Ls = ',   &
      ES11.3,' div_flx = ',ES11.3,' Z = ',ES12.4,' Test = ',ES11.3)
912 FORMAT(10X,A12,' DNDT = ',ES11.3,' Pr = ',ES11.3,' Ls = ',   &
      ES11.3,' div_flx = ',ES11.3,' Z = ',ES12.4,' N = ',ES11.3)
913 FORMAT(10X,'  CHEM : ITER=',I3,2X,' Species = ',A12,' Alt = ',F8.1, &
         ' Den = ',ES11.3' Pr = ',ES11.3,' Ls = ',ES11.3,' dNdt = ',ES11.3)
  STOP

END PROGRAM