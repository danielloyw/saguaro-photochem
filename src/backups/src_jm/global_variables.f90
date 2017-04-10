MODULE GLOBAL_VARIABLES

  USE PRECISION

  !  .. Planet Parameters

  REAL(RP), PARAMETER :: GM = 8.978E18            !  Grav Const X Planet Mass      
  REAL(RP), PARAMETER :: RPLANET = 2575.E5        !  Radius of Planet              
  REAL(RP), PARAMETER :: DH = 9.0576_RP           !  Planet Heliocentric Distance

  !  .. Photolysis averaging

  REAL(RP) :: cos_sza 
  REAL(RP) :: diurnal_average

  !  .. Run control parameters

  CHARACTER(len=3) :: runID
  LOGICAL :: lprnt
  REAL(RP) :: ascl
  INTEGER :: iaer
  REAL(RP) :: tstep_chem, tstep_diff

  ! .. Molecules

  INTEGER :: neutrmax
  INTEGER :: nionsmax
  INTEGER :: nsp
  INTEGER :: nchem
  INTEGER :: ndiff
  CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: name
  INTEGER, ALLOCATABLE, DIMENSION(:) :: istat
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: mmw
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhyd
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ncar
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n14n
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n15n
  INTEGER, ALLOCATABLE, DIMENSION(:) :: noxy
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dtype
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ad
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sd
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: phi
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sd_2
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sd_3
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ibnd
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: bval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ichrg
  INTEGER, ALLOCATABLE, DIMENSION(:) :: locp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: lopc
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ldcp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ldpc
  INTEGER :: iN2, iCH4, iELE
  LOGICAL :: lions

  !  .. Reactions 

  CHARACTER(len=256) :: nreact_file, ireact_file, treact_file
  CHARACTER(len=73), ALLOCATABLE, DIMENSION(:) :: ctitle
  INTEGER :: nrct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: itype
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: irct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rck
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ntab
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ntmp_rct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: tmp_rct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nprs_rct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: plog_rct
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: rct_tab

  !
  !  .. Atmosphere
  !

  INTEGER :: nlev, nibot, nzion
  REAL(RP) :: zbot
  REAL(RP) :: ed1,ed0,pk0,gek
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: z
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rz
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: grv
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: tn
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ti
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: te
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: prs
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: mass
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rho
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ek
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: col
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: den, den_old
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: ht
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rmid
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: dr
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: drp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rp2
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: rm2

  !  .. DIFCO

  REAL(RP), PARAMETER :: polN2 = 17.6E-25                              !  Polarizability of N2 
  REAL(RP), PARAMETER :: polCH4 = 26.0E-26                             !  Polarizability of CH4
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: df
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: alpha
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: beta
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: a
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: b
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: c
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: alphax
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: betax

  !  .. RATECO

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rt, rct

  !  .. PATHS

  REAL(RP) :: RSHADOW
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: ds
  INTEGER :: illum, nbot
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ntan

  !  .. Aerosols

  INTEGER :: naer, nawv, ihetero
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: kaer, tau_aer
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: surfarea


  !  .. Rayleigh Scattering
 
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: crs_ray
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: tau_ray

  !  .. External Production

  INTEGER :: next
  CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: name_ext
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: zext, hext, fext
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: prext
   
  !  .. PHOTO

  LOGICAL :: lcrsA, lcrsB, lcrsC, lcrsJ
  INTEGER :: nphrt, ncrsA, ncrsB, ncrsC, nabsA, nabsB, nabsC, nabsJ, nbrmaxA, nbrmaxB, nbrmaxC, nbrmaxJ
  INTEGER, ALLOCATABLE, DIMENSION(:) :: loabA
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnchA
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wcrsA
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: fsolA
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: xcrsA
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: bratA
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ionizeA
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: enrgIA
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: charge_stateA
  INTEGER, ALLOCATABLE, DIMENSION(:) :: loabB
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnchB
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wcrsB
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: fsolB
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: xcrsB
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: bratB
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ionizeB
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: enrgIB
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: charge_stateB
  INTEGER, ALLOCATABLE, DIMENSION(:) :: loabC
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnchC
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wcrsC
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: fsolC
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: xcrsC
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: bratC
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ionizeC
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: enrgIC
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: charge_stateC
  INTEGER, ALLOCATABLE, DIMENSION(:) :: loabJ
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnchJ
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: srateJ
  CHARACTER(len=87), ALLOCATABLE, DIMENSION(:) :: ptitle
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: irpt

  !  .. Variables used in ELCTRN, ELDEP1

  INTEGER :: nabs_el_thk, nabs_el_thn, nelb, nert
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: elctreV
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: elctDeV
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_tot_elast
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_tot_inelast
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: eCS
  CHARACTER(len=10), ALLOCATABLE, DIMENSION(:,:) :: state
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: enrgE
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipath
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: Sel
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_exc
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eCS_ion
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: pS
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: brat_el
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nbrnch_el
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: crs_tot_inel
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: sum_cs_rees_ion
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: sum_cs_rees_sec
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: esrc
  CHARACTER(len=87), ALLOCATABLE, DIMENSION(:) :: etitle
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iert
  
  ! .. SOLDEP1

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rph, rpt   
 
  !  .. ELDEP1

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eflux, eph, rpe

  !  .. CHEMEQ & COMPO

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr_chem, ls_chem
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rcdn
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr_ph, ls_ph
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr_pe, ls_pe
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pr, ls
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: flx
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: div_flx
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: dNdt

  !  .. COMPOUT

END MODULE GLOBAL_VARIABLES
