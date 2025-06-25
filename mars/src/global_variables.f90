! Defines all global variables that are used in program
module global_variables
  use types, only: wp => dp
  use constants

  !----------------------------------------------------------------------------
  !  Planet Parameters
  !----------------------------------------------------------------------------

  ! Planet = Mars
  real(wp), parameter :: m_planet = 6.417E26    ! planet mass
  real(wp), parameter :: GM = Gcnst * m_planet  ! grav const * planet mass
  real(wp), parameter :: rPlanet = 3389.5E5     ! planet radius
  real(wp), parameter :: DH = 1.520_wp          ! orbit semimajor axis in AU

  !----------------------------------------------------------------------------
  !  Run Parameters
  !----------------------------------------------------------------------------

  !  .. Photolysis averaging

  real(wp) :: cos_sza
  real(wp) :: diurnal_average

  !  .. Run control parameters

  character(len=3) :: runID
  logical :: lprnt
  logical :: lcompo_sol, lchemeq_sol
  real(wp) :: ascl ! not used
  integer :: iaer
  real(wp) :: tstep_chem, tstep_diff ! timesteps for reactions and diffusion

  !----------------------------------------------------------------------------
  !  Species
  !----------------------------------------------------------------------------

  ! number of neutrals
  integer :: n_neu
  ! number of ions
  integer :: n_ion
  ! total number of species (neutrals + ions + electrons)
  integer :: n_sp
  ! total species in chemical equilibrium
  integer :: n_chem
  ! total species in diffusive equilibrium
  integer :: n_diff
  ! species name (0th index = total)
  character(len=12), allocatable, dimension(:) :: sp_list
  ! how is species treated in model: dim=(species #)
  integer, allocatable, dimension(:) :: istat
  ! molecular weight (in amu): dim=(species #)
  real(wp), allocatable, dimension(:) :: mmw
  ! number of H, C, N-14, N-15, O in species: dim=(species #)
  integer, allocatable, dimension(:) :: nhyd, ncar, n14n, n15n, noxy

  ! how is diffusion treated: dim=(species #)
  integer, allocatable, dimension(:) :: dtype
  ! diffusion parameters A, s_1, phi, s_2, s_3: dim=(species #)
  real(wp), allocatable, dimension(:) :: Ad, sd, phi, sd_2, sd_3

  ! boundary type: dim=(species #, bottom/top)
  integer, allocatable, dimension(:,:) :: ibnd
  ! boundary velocity value (+ve upward): dim=(species #, bottom/top) 
  real(wp), allocatable, dimension(:,:) :: bval
  ! charge: dim=(species #)
  integer, allocatable, dimension(:) :: chrg
  
  ! are there ions in the model?
  logical :: have_ions
  
  ! Index mapping between lists
  ! list of chemical species -> list of all species: dim=(chemical species #)
  integer, allocatable, dimension(:) :: im_chem_all
  ! list of all species -> list of chemical species: dim=(species #)
  integer, allocatable, dimension(:) :: im_all_chem
  ! list of diffusing species -> list of all species: dim=(diffusive species #)
  integer, allocatable, dimension(:) :: im_diff_all
  ! list of all species -> list of diffusing species: dim=(species #)
  integer, allocatable, dimension(:) :: im_all_diff
  
  ! index of N2, CO2 and e
  integer :: iN2, iCO2, iELE

  !----------------------------------------------------------------------------
  !  Reactions
  !----------------------------------------------------------------------------

  ! number of reactions
  integer :: n_rct
  ! chemical equation: dim=(reaction #)
  character(len=73), allocatable, dimension(:) :: ctitle
  ! reaction type: dim=(reaction #)
  integer, allocatable, dimension(:) :: chem_type
  ! reactant index: dim=(reactant #, reaction #)
  integer, allocatable, dimension(:,:) :: irct
  ! reaction rate constants: dim=(constant 1-10, reaction #)
  real(wp), allocatable, dimension(:,:) :: rk
  
  ! Tabulated reaction rates variables (deactivated)
  ! integer, allocatable, dimension(:) :: ntab
  ! integer, allocatable, dimension(:) :: nprs_rct
  ! integer, allocatable, dimension(:) :: ntmp_rct
  ! real(wp), allocatable, dimension(:,:) :: plog_rct
  ! real(wp), allocatable, dimension(:,:) :: tmp_rct
  ! real(wp), allocatable, dimension(:,:,:) :: rct_tab

  !----------------------------------------------------------------------------
  !  Atmosphere
  !----------------------------------------------------------------------------

  ! number of elevation bins
  integer :: n_z
  ! altitude to opaque atmosphere layer (e.g., aerosols)
  real(wp) :: z_bot
  ! index of z_bot in n_z
  integer :: iz_bot
  ! coefficients for calculating eddy coefficient
  real(wp) :: ed1,ed0,pk0,gek
  ! altitude: unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: z
  ! radius (altitude + rPlanet): unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: rz
  ! gravity: unit=cm s-2; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: grv
  ! neutral temperature: unit=K; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: Tn
  ! ion temperature: unit=K; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: Ti
  ! electron temperature: unit=K; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: Te
  ! pressure: unit=dyne cm-2 | 0.1 Pa | 1E-6 bar; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: prs
  ! mass density: unit=g cm-3; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: rho
  ! mean molecular weight: unit=amu; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: mass
  ! eddy diffusion coefficient: unit=cm2 s-1; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: eK
  ! column density above specified altitude: unit=cm-2; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: col
  ! number density: unit=cm-3; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: den, den_old
  ! volume mixing ratio: unit=unitless; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: vmr
  ! scale height: unit=cm; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: Ht
  ! middle of rz bins: unit=cm; dim=(# of altitude levels)
  real(wp), allocatable, dimension(:) :: rmid
  ! rmid(i)-rmid(i-1): unit=cm; dim=(# of altitude levels)
  real(wp), allocatable, dimension(:) :: dr
  ! rz(i)-rz(i-1): unit=cm; dim=(# of altitude levels)
  real(wp), allocatable, dimension(:) :: drp
  ! (rmid(nl)/rz(nl))**2: unit=cm; dim=(# of altitude levels)
  real(wp), allocatable, dimension(:) :: rp2
  ! (rmid(nl-1)/rz(nl))**2: unit=cm; dim=(# of altitude levels)
  real(wp), allocatable, dimension(:) :: rm2

  !----------------------------------------------------------------------------
  !  DIFCO
  !----------------------------------------------------------------------------

  real(wp), parameter :: polN2 = 17.6E-25               ! Polarizability of N2
  real(wp), parameter :: polCO2 = 2.63E-24              ! Polarizability of CO2
  real(wp), allocatable, dimension(:,:) :: df           ! diffusion coefficient
  !
  real(wp), allocatable, dimension(:,:) :: alpha
  !
  real(wp), allocatable, dimension(:,:) :: beta
  ! divergence coefficient for bin below
  real(wp), allocatable, dimension(:,:) :: a
  ! divergence coefficient for bin
  real(wp), allocatable, dimension(:,:) :: b
  ! divergence coefficient for bin above
  real(wp), allocatable, dimension(:,:) :: c
  !
  real(wp), allocatable, dimension(:,:) :: alphax
  !
  real(wp), allocatable, dimension(:,:) :: betax
  !
  real(wp), allocatable, dimension(:) :: ekp
  !
  real(wp), allocatable, dimension(:,:) :: dfp

  !----------------------------------------------------------------------------
  !  RATECO
  !----------------------------------------------------------------------------

  ! rate (reaction #, altitude bin #),
  real(wp), allocatable, dimension(:,:) :: rt, rct

  !----------------------------------------------------------------------------
  !  PATHS
  !----------------------------------------------------------------------------

  ! planet radius + z_bot
  real(wp) :: RSHADOW
  ! path through atmosphere for RT (altitude, tangent altitude)
  real(wp), allocatable, dimension(:,:) :: ds
  ! illumination condition: 1=day, 0=twilight, -1=night
  integer :: illum
  ! index of lowest illuminated altitude
  integer :: nbot
  ! altitude index for tangent altitude
  integer, allocatable, dimension(:) :: ntan

  !----------------------------------------------------------------------------
  !  Aerosols
  !----------------------------------------------------------------------------

  integer :: naer, nawv, ihetero !ihetero not used
  real(wp), allocatable, dimension(:,:) :: kaer, tau_aer
  real(wp), allocatable, dimension(:) :: surfarea

  !----------------------------------------------------------------------------
  !  Rayleigh Scattering
  !----------------------------------------------------------------------------

  ! scattering cross-section in spectral range C
  real(wp), allocatable, dimension(:) :: crs_ray
  ! tau (altitude, wavelength)
  real(wp), allocatable, dimension(:,:) :: tau_ray

  !----------------------------------------------------------------------------
  !  External Production
  !----------------------------------------------------------------------------

  integer :: next
  character(len=12), allocatable, dimension(:) :: name_ext
  real(wp), allocatable, dimension(:) :: zext, hext, fext
  real(wp), allocatable, dimension(:,:) :: prext ! external production

  !----------------------------------------------------------------------------
  !  PHOTO
  !----------------------------------------------------------------------------

  ! calculate photolysis?
  logical :: lcrsA, lcrsB, lcrsC, lcrsJ
  ! total number of photolysis reactions
  integer :: nphrt
  ! number of wavelength bins for spectral ranges A, B and C
  integer :: ncrsA, ncrsB, ncrsB_low, ncrsC
  ! number of species for spectral ranges A, B and C
  integer :: nabsA, nabsB, nabsC, nabsJ
  ! maximum number of branches for each species
  integer :: nbrmaxA, nbrmaxB, nbrmaxC, nbrmaxJ
  ! photolyzed species index
  integer, allocatable, dimension(:) :: loabA, loabB, loabC, loabJ
  ! number of photo reactions/branches
  integer, allocatable, dimension(:) :: nbrnchA, nbrnchB, nbrnchC, nbrnchJ
  ! wavelength scale
  real(wp), allocatable, dimension(:) :: wcrsA, wcrsB, wcrsB_low, wcrsC
  ! solar flux for spectral ranges A, B and C
  real(wp), allocatable, dimension(:) :: fsolA, fsolB, fsolC
  ! transmission of solar flux for spectral ranges A, B and C
  real(wp), allocatable, dimension(:,:) :: trnA, trnB, trnC
  ! non-diurnally averaged photo rate for spectral ranges A, B and C
  ! (wavelength, branch #, species #)
  real(wp), allocatable, dimension(:,:,:) :: prtA, prtB, prtC
  ! total absorption cross sections (wavelength, species #)
  real(wp), allocatable, dimension(:,:) :: xcrsA, xcrsB, xcrsC
  ! branch ratio (wavelength, branch #, species #)
  real(wp), allocatable, dimension(:,:,:) :: bratA, bratB, bratC
  ! are there any ions?
  logical, allocatable, dimension(:,:) :: ionizeA, ionizeB, ionizeC
  ! threshold energy (branch #, species #) in angstroms
  real(wp), allocatable, dimension(:,:) :: enrgIA, enrgIB, enrgIC
  ! number of ions (branch #, species #)
  real(wp), allocatable, dimension(:,:) :: charge_stateA
  real(wp), allocatable, dimension(:,:) :: charge_stateB
  real(wp), allocatable, dimension(:,:) :: charge_stateC
  ! specified optically thin photolysis rates
  real(wp), allocatable, dimension(:,:) :: srateJ
  ! formula for each photolysis reaction
  character(len=87), allocatable, dimension(:) :: ptitle
  ! indices for reactants/products
  integer, allocatable, dimension(:,:) :: irpt

  !----------------------------------------------------------------------------
  !  ELCTRN, ELDEP1
  !----------------------------------------------------------------------------

  ! number of species (thick), number of species (thin)
  integer :: nabs_el_thk, nabs_el_thn
  ! number of bins, number of electron impact reactions
  integer :: nelb, nert
  ! mean bin energy
  real(wp), allocatable, dimension(:) :: elctreV
  ! bin width
  real(wp), allocatable, dimension(:) :: elctDeV
  ! elastic cross section (energy bin, species)
  real(wp), allocatable, dimension(:,:) :: eCS_tot_elast
  ! inelastic cross section (energy bin, species)
  real(wp), allocatable, dimension(:,:) :: eCS_tot_inelast
  ! cross-sections (energy bin, species, state)
  real(wp), allocatable, dimension(:,:,:) :: eCS
  ! name of excited state (species, state)
  character(len=10), allocatable, dimension(:,:) :: state
  ! energy level of excited state (species, state)
  real(wp), allocatable, dimension(:,:) :: enrgE
  ! number of states (species, excited/dissociation/ionization)
  integer, allocatable, dimension(:,:) :: ipath
  ! electron production/source (altitude, bin)
  real(wp), allocatable, dimension(:,:) :: Sel
  ! total cross-section for excitation and dissociation
  real(wp), allocatable, dimension(:,:) :: eCS_exc
  ! total cross-section for ionization
  real(wp), allocatable, dimension(:,:) :: eCS_ion
  ! 0?
  real(wp), allocatable, dimension(:) :: pS
  ! cross section for thin?
  real(wp), allocatable, dimension(:,:,:) :: brat_el
  ! number of branches for dissociation + ionization (species)
  integer, allocatable, dimension(:) :: nbrnch_el
  ! total cross section for thin
  real(wp), allocatable, dimension(:,:) :: crs_tot_inel
  !
  real(wp), allocatable, dimension(:,:,:) :: sum_cs_rees_ion
  ! total cross-section for secondary electron production
  real(wp), allocatable, dimension(:,:,:) :: sum_cs_rees_sec
  ! electron production rate (reaction, altitude, energy bin)
  real(wp), allocatable, dimension(:,:,:) :: esrc
  ! equation (reaction)
  character(len=87), allocatable, dimension(:) :: etitle
  ! reactant/product index (reactant/product #, reaction)
  integer, allocatable, dimension(:,:) :: iert

  !----------------------------------------------------------------------------
  !  SOLDEP1
  !----------------------------------------------------------------------------

  ! absorption rate per molecule (reaction, altitude)
  real(wp), allocatable, dimension(:,:) :: rph
  ! total absorption rate (reaction, altitude)
  real(wp), allocatable, dimension(:,:) :: rpt

  !----------------------------------------------------------------------------
  !  ELDEP1
  !----------------------------------------------------------------------------

  real(wp), allocatable, dimension(:,:) :: eflux, eph, rpe

  !----------------------------------------------------------------------------
  !  CHEMEQ & COMPO
  !----------------------------------------------------------------------------

  !  Production and Loss from chemical processes
  real(wp), allocatable, dimension(:,:) :: pr_chem, ls_chem
  !  Condensation loss
  real(wp), allocatable, dimension(:,:) :: rcdn
  !  Production and Loss from photon processes
  real(wp), allocatable, dimension(:,:) :: pr_ph, ls_ph
  !  Production and Loss from electron processes
  real(wp), allocatable, dimension(:,:) :: pr_pe, ls_pe
  !  Net production and loss
  real(wp), allocatable, dimension(:,:) :: pr, ls
  !  Flux
  real(wp), allocatable, dimension(:,:) :: flx
  !  Flux divergence
  real(wp), allocatable, dimension(:,:) :: div_flx
  !  Net balance
  real(wp), allocatable, dimension(:,:) :: dNdt

end module global_variables
