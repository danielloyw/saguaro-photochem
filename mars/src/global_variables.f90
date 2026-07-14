! Defines all global variables that are used in program
module global_variables
  use types, only: wp => dp
  use constants

  !----------------------------------------------------------------------------
  !  Planet Parameters
  !----------------------------------------------------------------------------

  ! Planet = Mars
  ! planet mass: unit=g
  real(wp), parameter :: m_planet = 6.417E26
  ! grav const * planet mass
  real(wp), parameter :: GM = Gcnst * m_planet
  ! planet radius: unit=cm
  real(wp), parameter :: rPlanet = 3389.5E5
  ! orbit semimajor axis: unit=AU
  real(wp), parameter :: DH = 1.520_wp

  !----------------------------------------------------------------------------
  !  Run Parameters
  !----------------------------------------------------------------------------

  ! name of run
  character(len=3) :: runID
  ! time step size for diffusion loop iterations: unit=s
  real(wp) :: tstep_diff
  ! time step size for chemistry loop iterations: unit=s
  real(wp) :: tstep_chem
  ! recalculate photo rates in diffusion loop?
  logical :: recal_photo_diff
  ! recalculate photo rates in chemistry loop?
  logical :: recal_photo_chem
  ! print intermediate chemistry loops?
  logical :: print_chem
  
  !----------------------------------------------------------------------------
  !  Species
  !----------------------------------------------------------------------------

  ! number of neutrals
  integer :: n_neu
  ! number of ions
  integer :: n_ion
  ! total number of species (neutrals + ions + electrons)
  integer :: n_sp
  ! total species for diffusion loop
  integer :: n_diff
  ! total species for chemistry loop
  integer :: n_chem
  ! species name (0th index = total)
  character(len=12), allocatable, dimension(:) :: sp_list
  ! how is species treated in model: dim=(species #)
  integer, allocatable, dimension(:) :: istat
  ! molecular weight: unit=amu; dim=(species #)
  real(wp), allocatable, dimension(:) :: mmw
  ! number of H, C, N-14, N-15, O in species: dim=(species #)
  integer, allocatable, dimension(:) :: nhyd, ncar, n14N, n15N, noxy

  ! how is diffusion treated: dim=(species #)
  integer, allocatable, dimension(:) :: dtype
  ! parameters for calculating diffusion constant: dim=(species #)
  real(wp), allocatable, dimension(:) :: Ad, sd1, phi, sd2, sd3

  ! boundary type: dim=(species #, bottom/top)
  integer, allocatable, dimension(:,:) :: ibnd
  ! boundary velocity value (+ve upward): unit=cm s-1;
  ! dim=(species #, bottom/top) 
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
  !  Atmosphere
  !----------------------------------------------------------------------------

  ! number of elevation bins
  integer :: n_z
  ! altitude to opaque atmosphere layer (e.g., aerosols): unit=cm
  real(wp) :: z_bot
  ! index of z_bot in n_z
  integer :: iz_bot
  ! altitude: unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: z
  ! radius (altitude + rPlanet): unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: rz
  ! gravity: unit=cm s-2; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: grv
  ! neutral temperature: unit=K; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: Tn
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
  real(wp), allocatable, dimension(:,:) :: den
  ! number density from previous TIME loop iteration: 
  ! unit=cm-3; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: den_old
  ! volume mixing ratio: unit=unitless; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: vmr
  ! scale height: unit=cm; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: Ht
  ! middle of rz bins, (rz(i+1) + rz(i))/2: unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: rz_mid
  ! rz_mid(i)-rz_mid(i-1): unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: dr
  ! rz(i)-rz(i-1): unit=cm; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: dr_mid

  !----------------------------------------------------------------------------
  !  Diffusion
  !----------------------------------------------------------------------------

  ! diffusion coefficient: unit=cm2 s-1; dim(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: Df
  ! flux equation constants: unit=cm s-1; 
  ! dim(altitude level, diffusive species #)
  real(wp), allocatable, dimension(:,:) :: alpha, beta
  ! divergence coefficient for bin below: unit=s-1; 
  ! dim(altitude level, diffusive species #)
  real(wp), allocatable, dimension(:,:) :: a
  ! divergence coefficient for bin: unit=s-1; 
  ! dim(altitude level, diffusive species #)
  real(wp), allocatable, dimension(:,:) :: b
  ! divergence coefficient for bin above: unit=s-1; 
  ! dim(altitude level, diffusive species #)
  real(wp), allocatable, dimension(:,:) :: c

  !----------------------------------------------------------------------------
  !  Chemical Reactions
  !----------------------------------------------------------------------------

  ! number of reactions
  integer :: n_rct
  ! chemical equation: dim=(reaction #)
  character(len=73), allocatable, dimension(:) :: ctitle
  ! reaction type: dim=(reaction #)
  integer, allocatable, dimension(:) :: chem_type
  ! reactant index: dim=(reactant #, reaction #)
  integer, allocatable, dimension(:,:) :: ireactant
  
  ! Tabulated reaction rates variables (deactivated)
  ! integer, allocatable, dimension(:) :: ntab
  ! integer, allocatable, dimension(:) :: nprs_rct
  ! integer, allocatable, dimension(:) :: ntmp_rct
  ! real(wp), allocatable, dimension(:,:) :: plog_rct
  ! real(wp), allocatable, dimension(:,:) :: tmp_rct
  ! real(wp), allocatable, dimension(:,:,:) :: rct_tab

  ! rate coefficients: unit variable; dim=(reaction #, altitude level)
  real(wp), allocatable, dimension(:,:) :: rk
  ! reaction rate: unit=cm-3 s-1; dim=(reaction #, altitude level)
  real(wp), allocatable, dimension(:,:) :: rate_rct
  
  !----------------------------------------------------------------------------
  !  Photo
  !----------------------------------------------------------------------------

  ! cosine of solar zenith angle
  real(wp) :: cos_sza
  ! sine of solar zenith angle
  real(wp) :: sin_sza
  ! factor for diurnal averaging (usually 0.5 for day-night)
  real(wp) :: diurnal_average
  
  ! calculate photolysis?
  logical :: cal_photoA, cal_photoB, cal_photoC, cal_photoJ
  ! total number of photolysis reactions
  integer :: n_prct
  ! wavelength scale for spectral regions A, B, C (lres for low-resolution): 
  ! unit=angstrom; dim=(wavelength bin)
  real(wp), allocatable, dimension(:) :: waveA, waveB, waveB_lres, waveC
  ! number of wavelength bins for spectral regions A, B, C (lres for
  ! low-resolution)
  integer :: n_waveA, n_waveB, n_waveB_lres, n_waveC
  ! number of photo species for spectral regions A, B, C, and specified
  ! reactions
  integer :: n_sp_photoA, n_sp_photoB, n_sp_photoC, n_sp_photoJ
  ! index mapping from list of photo species -> list of all species: 
  ! dim=(photolyzed species #)
  integer, allocatable, dimension(:) :: im_photoA_all, im_photoB_all
  integer, allocatable, dimension(:) :: im_photoC_all, im_photoJ_all
  ! number of photo reactions/branches: dim=(photolyzed species #)
  integer, allocatable, dimension(:) :: n_branchA, n_branchB
  integer, allocatable, dimension(:) :: n_branchC, n_branchJ
  ! equation for photo reaction: dim=(photo reaction #)
  character(len=87), allocatable, dimension(:) :: ptitle
  ! index mapping from list of species in photo reactions -> list of all
  ! species: dim=(species # in reaction, photo reaction #)
  integer, allocatable, dimension(:,:) :: im_photo_all
  ! threshold energy: unit=angstrom; dim=(branch #, photolyzed species #)
  real(wp), allocatable, dimension(:,:) :: enrgIA, enrgIB, enrgIC
  ! number of ions produced in reaction (branch #, photolyzed species #)
  real(wp), allocatable, dimension(:,:) :: charge_stateA
  real(wp), allocatable, dimension(:,:) :: charge_stateB
  real(wp), allocatable, dimension(:,:) :: charge_stateC
  ! total absorption cross sections: unit=cm2
  ! dim=(wavelength bin, photolyzed species #)
  real(wp), allocatable, dimension(:,:) :: csA, csB, csC
  ! branch ratio: dim=(wavelength bin, branch #, photolyzed species #)
  real(wp), allocatable, dimension(:,:,:) :: branch_ratioA
  real(wp), allocatable, dimension(:,:,:) :: branch_ratioB
  real(wp), allocatable, dimension(:,:,:) :: branch_ratioC
  ! specified optically thin photolysis rates: unit=
  ! dim=(branch #, photolyzed species #)
  real(wp), allocatable, dimension(:,:) :: srateJ
  
  ! solar flux for spectral ranges A, B, C: 
  ! unit=photons cm-2 s-1 (bin width in nm)-1; dim=(wavelength bin)
  real(wp), allocatable, dimension(:) :: sol_fluxA, sol_fluxB, sol_fluxC
  ! transmission of solar flux for spectral regions A, B, C: unit=unitless;
  ! dim=(wavelength bin, altitude level)
  real(wp), allocatable, dimension(:,:) :: trnA, trnB, trnC
  ! non-diurnally averaged photo rate for spectral regions A, B, C:
  ! unit=s-1 (bin width in nm)-1; dim=(wavelength bin, branch #, species #)
  real(wp), allocatable, dimension(:,:,:) :: prtA, prtB, prtC
  
  ! Rayleigh scattering cross-section in spectral range C: unit=cm2
  real(wp), allocatable, dimension(:) :: cs_ray
  ! tau: unit=unitless; dim=(altitude level, wavelength bin)  
  real(wp), allocatable, dimension(:,:) :: tau_ray
  
  ! illumination condition: 1=day, 0=twilight, -1=night
  integer :: illum
  ! index of lowest illuminated altitude level
  integer :: ibot
  ! path through atmosphere for RT: unit=cm;
  ! dim=(altitude level, tangent altitude level)
  real(wp), allocatable, dimension(:,:) :: ds
  ! altitude index for tangent altitude
  integer, allocatable, dimension(:) :: itan
  
  ! absorption rate per molecule (reaction, altitude): unit=s-1
  ! dim=(photo reaction #, altitude level)
  real(wp), allocatable, dimension(:,:) :: freq_prct
  ! photo reaction rate: unit=cm-3 s-1; dim=(photo reaction #, altitude level)
  real(wp), allocatable, dimension(:,:) :: rate_prct

  !----------------------------------------------------------------------------
  !  Electrons
  !----------------------------------------------------------------------------

  ! number of electron impact reactions
  integer :: n_erct
  ! energy scale for electrons: unit=eV; dim=(energy bin)
  real(wp), allocatable, dimension(:) :: e_enrg
  ! number of bins in scale
  integer :: n_e_enrg
  ! energy bin width: unit=eV; dim=(energy bin)
  real(wp), allocatable, dimension(:) :: e_denrg
  ! number of species (thick), number of species (thin)
  integer :: n_sp_e_thk, n_sp_e_thn
  ! number of states: dim=(excited/dissociation/ionization, e species #)
  integer, allocatable, dimension(:,:) :: n_states
  ! number of dissociation + ionization reactions: dim=(e species #)
  integer, allocatable, dimension(:) :: n_branch_e
  ! equation for electron reaction: dim=(electron reaction #)
  character(len=87), allocatable, dimension(:) :: etitle
  ! energy level of excited state: unit=eV; dim=(state, e species #)
  real(wp), allocatable, dimension(:,:) :: enrg_state
  ! index mapping from list of electron impact species -> list of all species
  integer, allocatable, dimension(:) :: im_e_reactant_all
  ! index mapping from list of species in electron reactions -> list of
  ! all species: dim=(species # in reaction, electron reaction #)
  integer, allocatable, dimension(:,:) :: im_e_all
  ! electron cross sections: unit=cm2; dim=(energy bin, state, e species #)
  real(wp), allocatable, dimension(:,:,:) :: ecs
  ! total cross-section for ionization: unit=cm2;
  ! dim=(energy bin, state, e species #)
  real(wp), allocatable, dimension(:,:,:) :: cs_ion
  ! total cross-section for secondary electron production: unit=cm2
  ! dim=(energy bin, state, e species #)
  real(wp), allocatable, dimension(:,:,:) :: cs_sec
  ! total cross-section for excitation and dissociation: unit=cm2;
  ! dim=(energy bin, e species #)
  ! real(wp), allocatable, dimension(:,:) :: ecs_exc
  ! total cross-section for ionization: unit=cm2; dim=(energy bin, e species #)
  ! real(wp), allocatable, dimension(:,:) :: ecs_ion

  ! photoelectron production rate: unit=cm-3 s-1 eV-1;
  ! dim=(photo reaction, altitude level, energy bin)
  real(wp), allocatable, dimension(:,:,:) :: esrc
  ! photoelectron production rate: unit=cm-3 s-1 eV-1;
  ! dim=(altitude, energy bin)
  real(wp), allocatable, dimension(:,:) :: Sel
  ! omnidirectional electron flux: unit=cm-2 s-1 eV-1;
  ! dim=(altitude level, energy bin)
  real(wp), allocatable, dimension(:,:) :: eflux
  ! electron reaction frequency: unit=s-1;
  ! dim=(electron reaction #, altitude level)
  real(wp), allocatable, dimension(:,:) :: freq_erct
  ! electron reaction rate: unit=cm-3 s-1;
  ! dim=(electron reaction #, altitude level)
  real(wp), allocatable, dimension(:,:) :: rate_erct
  ! real(wp), allocatable, dimension(:) :: pS
  
  !----------------------------------------------------------------------------
  !  Balance
  !----------------------------------------------------------------------------

  ! production and loss from chemical processes: unit=cm-3 s-1;
  ! dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: pr_c, ls_c
  ! production and loss from photo processes: unit=cm-3 s-1;
  ! dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: pr_p, ls_p
  ! production and loss from electron processes: unit=cm-3 s-1;
  ! dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: pr_e, ls_e
  ! loss from condensation: unit=cm-3 s-1;
  ! dim=(altitude level, species #)
  ! real(wp), allocatable, dimension(:,:) :: ls_cdn
  ! net production and loss: unit=cm-3 s-1;
  ! dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: pr, ls
  ! diffusive flux: unit=cm-2 s-1; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: flux
  ! flux divergence: unit=cm-3 s-1; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: div_flux
  ! net balance: unit=cm-3 s-1; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: dNdt

end module
