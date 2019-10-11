MODULE GLOBAL_VARIABLES

  USE PRECISION

  !  .. Planet Parameters

  REAL(RP), PARAMETER :: GM = 4.283E19            !  Grav Const X Planet Mass      
  REAL(RP), PARAMETER :: RPLANET = 3389.5E5        !  Radius of Planet              
  REAL(RP), PARAMETER :: DH = 1.5200_RP           !  Planet Heliocentric Distance (AU)

  !  .. Photolysis averaging

  REAL(RP) :: cos_sza 
  REAL(RP) :: diurnal_average

  !  .. Run control parameters

  CHARACTER(len=3) :: runID
  LOGICAL :: lprnt
  REAL(RP) :: ascl
  INTEGER :: iaer
  REAL(RP) :: tstep_chem, tstep_diff ! timesteps for reactions and diffusion

  ! .. Molecules

  INTEGER :: neutrmax   ! number of neutrals
  INTEGER :: nionsmax   ! number of ions
  INTEGER :: nsp        ! total species
  INTEGER :: nchem      ! total species in chem eqm
  INTEGER :: ndiff      ! total species in diffusion eqm
  CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: name  ! name of species
  INTEGER, ALLOCATABLE, DIMENSION(:) :: istat       ! how is species treated in model
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: mmw        ! molecular weight (in amu)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhyd        ! number of H
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ncar        ! number of C
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n14n        ! number of N-14
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n15n        ! number of N-15
  INTEGER, ALLOCATABLE, DIMENSION(:) :: noxy        ! number of O
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dtype       ! how is diffusion treated
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ad         ! diffusion parameter A
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sd         ! diffusion parameter s_1
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: phi        ! diffusion parameter phi
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sd_2       ! diffusion parameter s_2
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sd_3       ! diffusion parameter s_3
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ibnd      ! boundary type (species #, bottom/top)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: bval     ! boundary density value (species #, bottom/top)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ichrg       ! charge
  INTEGER, ALLOCATABLE, DIMENSION(:) :: locp        ! index conversion from list of chemical eqm species onto list of total species
  INTEGER, ALLOCATABLE, DIMENSION(:) :: lopc        ! index conversion from list of total species onto list of chemical eqm species
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ldcp        ! index conversion from list of diffusing species onto list of total species
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ldpc        ! index conversion from list of total species onto list of diffusing species
  INTEGER :: iN2, iCO2, iELE                        ! index of N2, CO2 and e
  LOGICAL :: lions      ! are there ions in the model?

  !  .. Reactions 

  CHARACTER(len=73), ALLOCATABLE, DIMENSION(:) :: ctitle        ! chemical equation
  INTEGER :: nrct   ! number of reactions
  INTEGER, ALLOCATABLE, DIMENSION(:) :: itype           ! reaction type
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: irct          ! reactant index (reactant #, reaction #)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rck          ! reaction rate constants (constant 1-10, reaction #)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ntab            ! number of tabulated reactions
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ntmp_rct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: tmp_rct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nprs_rct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: plog_rct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: rct_tab

  !  .. Atmosphere

  INTEGER :: nlev, nibot            ! number of elevation bins, index of zbot
  REAL(RP) :: zbot                  ! lower altitude bound for model
  REAL(RP) :: ed1,ed0,pk0,gek       ! coefficients for calculating eddy coefficient
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: z              ! altitude
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rz             ! radius
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: grv            ! gravity
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: tn             ! neutral temp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ti             ! ion temp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: te             ! electron temp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: prs            ! pressure
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: mass           ! mean molecular weight
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rho            ! mass density
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ek             ! eddy coeff
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: col            ! column density
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: den, den_old ! number density
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: xmol         ! mole density
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: ht           ! scale height
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rmid           ! middle of rz bins
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: dr             ! rmid(i)-rmid(i-1)
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: drp            ! rz(i)-rz(i-1)
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rp2            ! (rmid(nl)/rz(nl))**2
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rm2            ! (rmid(nl-1)/rz(nl))**2

  !  .. DIFCO

  REAL(RP), PARAMETER :: polN2 = 17.6E-25               ! Polarizability of N2
  REAL(RP), PARAMETER :: polCO2 = 26.3E-26              ! Polarizability of CH4 = 26.0E-26, need to know CO2 = 26.3
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: df           ! diffusion coefficient
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: alpha
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: beta
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: a            ! divergence coefficient for bin below
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: b            ! divergence coefficient for bin
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: c            ! divergence coefficient for bin above
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: alphax
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: betax

  !  .. RATECO

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rt, rct      ! rate (reaction #, altitude bin #),

  !  .. PATHS

  REAL(RP) :: RSHADOW   ! planet radius + zbot
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: ds           ! path through atmosphere for RT (altitude, tangent altitude)
  INTEGER :: illum, nbot                                ! day(1)/twilight(0)/night(-1), index of lowest illuminated altitude
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ntan            ! altitude index for tangent altitude

  !  .. Aerosols

  INTEGER :: naer, nawv, ihetero
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: kaer, tau_aer
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: surfarea


  !  .. Rayleigh Scattering
 
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: crs_ray ! scattering cross-section in spectral range C
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: tau_ray ! tau (altitude, wavelength)

  !  .. External Production

  INTEGER :: next
  CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: name_ext
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: zext, hext, fext
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: prext ! external production
   
  !  .. PHOTO

  LOGICAL :: lcrsA, lcrsB, lcrsC, lcrsJ                 ! calculate photolysis?
  INTEGER :: nphrt                                      ! total number of photolysis reactions
  INTEGER :: ncrsA, ncrsB, ncrsB_low, ncrsC                        ! number of wavelength bins for spectral ranges A, B and C
  INTEGER :: nabsA, nabsB, nabsC, nabsJ                 ! number of species for spectral ranges A, B and C
  INTEGER :: nbrmaxA, nbrmaxB, nbrmaxC, nbrmaxJ         ! maximum number of branches for each species
  INTEGER, ALLOCATABLE, DIMENSION(:) :: loabA, loabB, loabC, loabJ                      ! photolyzed species index
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnchA, nbrnchB, nbrnchC, nbrnchJ              ! number of photo reactions/branches
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wcrsA, wcrsB, wcrsB_low, wcrsC                            ! wavelength scale
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: fsolA, fsolB, fsolC                            ! solar flux for spectral ranges A, B and C
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: trnA, trnB, trnC                            ! transmission of solar flux for spectral ranges A, B and C
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: prtA, prtB, prtC                            ! non-diurnally averaged photo rate for spectral ranges A, B and C (wavelength, branch #, species #)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: xcrsA, xcrsB, xcrsC                          ! total absorption cross sections (wavelength, species #)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: bratA, bratB, bratC                        ! branch ratio (wavelength, branch #, species #)
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ionizeA, ionizeB, ionizeC                     ! are there any ions?
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: enrgIA, enrgIB, enrgIC                       ! threshold energy (branch #, species #) in angstroms
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: charge_stateA, charge_stateB, charge_stateC  ! number of ions (branch #, species #)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: srateJ                                       ! specified optically thin photolysis rates
  CHARACTER(len=87), ALLOCATABLE, DIMENSION(:) :: ptitle                                ! formula for each photolysis reaction
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: irpt                                          ! indices for reactants/products

  !  .. Variables used in ELCTRN, ELDEP1

  INTEGER :: nabs_el_thk, nabs_el_thn, nelb, nert ! number of species (thick), number of species (thin), number of bins, number of electron impact reactions
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: elctreV ! mean bin energy
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: elctDeV ! bin width
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_tot_elast ! elastic cross section (energy bin, species)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_tot_inelast ! inelastic cross section (energy bin, species)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: eCS ! cross-sections (energy bin, species, state)
  CHARACTER(len=10), ALLOCATABLE, DIMENSION(:,:) :: state ! name of excited state (species, state)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: enrgE ! energy level of excited state (species, state)
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipath ! number of states (species, excited/dissociation/ionization)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: Sel ! electron production/source (altitude, bin)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_exc ! total cross-section for excitation and dissociation
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_ion ! total cross-section for ionization
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: pS ! 0?
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: brat_el ! cross section for thin?
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnch_el ! number of branches for dissociation + ionization (species)
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: crs_tot_inel ! total cross section for thin
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: sum_cs_rees_ion
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: sum_cs_rees_sec ! total cross-section for secondary electron production
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: esrc ! electron production rate (reaction, altitude, energy bin)
  CHARACTER(len=87), ALLOCATABLE, DIMENSION(:) :: etitle ! equation (reaction)
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iert ! reactant/product index (reactant/product #, reaction)
  
  ! .. SOLDEP1

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rph, rpt ! absorption rate per molecule (reaction, altitude), total absorption rate (reaction, altitude)
 
  !  .. ELDEP1

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eflux, eph, rpe

  !  .. CHEMEQ & COMPO

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr_chem, ls_chem !  Production and Loss from chemical processes
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rcdn             !  Condensation loss
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr_ph, ls_ph     !  Production and Loss from photon processes
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr_pe, ls_pe     !  Production and Loss from electron processes
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr, ls           !  Net production and loss
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: flx              !  Flux
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: div_flx          !  Flux divergence
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: dNdt             !  Net balance

  !  .. COMPOUT

END MODULE GLOBAL_VARIABLES
