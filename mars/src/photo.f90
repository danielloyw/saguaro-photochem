module photo_mod

contains

subroutine photo
! This subroutine sets up the calculation for photolysis reactions, including
! the reading in of the photo reactions (and their cross sections) and the
! solar flux. Reaction rates are calculated using cross sections from 0-2000 A,
! divided into 3 wavelength regions:
! Region A: 0-750 A
! Region B: 750-1080 A
! Region C: 1080-2000 A
! Region B has significantly higher spectral resolution to account for sharp
! peaks in the CO and N2 cross sections from pre-dissociating states. 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : intrp, locate

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  implicit none
  ! input file paths
  character(len=128) :: file_photoA, file_photoC, file_jvals
  character(len=128) ::  file_solA, file_solB, file_solC, file_sol_hres
  
  ! high resolution wavelength bin size for region B
  real(wp), allocatable, dimension(:) :: dwaveB
  ! equation for photo reaction
  character(len=87), allocatable, dimension(:,:) :: ptitleA, ptitleB
  character(len=87), allocatable, dimension(:,:) :: ptitleC, ptitleJ
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, allocatable, dimension(:,:,:) :: im_photoA_product_all
  integer, allocatable, dimension(:,:,:) :: im_photoB_product_all
  integer, allocatable, dimension(:,:,:) :: im_photoC_product_all
  integer, allocatable, dimension(:,:,:) :: im_photoJ_product_all
  ! wavelength scale for solar fluxes over regions A, B and C
  real(wp), allocatable, dimension(:) :: sol_waveA, sol_waveB_lres, sol_waveC
  ! low-resolution solar fluxes over region B
  real(wp), allocatable, dimension(:) :: sol_fluxB_lres
  ! native high-resolution wavelength scale for solar fluxes
  real(wp), allocatable, dimension(:) :: sol_wave_hres
  ! high-resolution solar fluxes with native wavelength scale
  real(wp), allocatable, dimension(:) :: sol_flux_hres

  ! loop variables
  integer :: i_sp, i_branch, i_wave
  ! temporary / dummy variables
  integer :: ti1, nr
  
  !----------------------------------------------------------------------------
  !  Paths to data files
  !----------------------------------------------------------------------------
    
  file_photoA   = '../data/photons/photoA.dat'
  file_photoC   = '../data/photons/photoC.dat'
  file_jvals    = '../data/photons/jvals.dat'
  file_solA     = '../data/solar/solar_fluxA.dat'
  file_solB     = '../data/solar/solar_fluxB.dat'
  file_solC     = '../data/solar/solar_fluxC.dat'
  file_sol_hres = '../data/solar/sumer.dat'
  
  !----------------------------------------------------------------------------
  !  Create high-resolution wavelength scale for region B
  !----------------------------------------------------------------------------
  
  n_waveB = 1650000
  allocate(waveB(n_waveB), dwaveB(n_waveB))
  dwaveB = 2.E-4_wp
  do concurrent (i_wave = 1:n_waveB)
    waveB(i_wave) = 750._wp + (i_wave-half) * dwaveB(i_wave)
  end do
  
  !----------------------------------------------------------------------------
  !  Read photolysis cross sections, branching ratios, and rates
  !----------------------------------------------------------------------------
  
  ! Region A (0-750 A)
  call read_photo(file_photoA, sp_list, & ! input variables
    im_photoA_all, n_branchA, ptitleA, enrgIA, is_ionizationA, charge_stateA, & 
    im_photoA_product_all, waveA, crsA, branch_ratioA)
  
  n_waveA = size(waveA)
  n_sp_photoA = size(crsA, 2)
  n_branch_maxA = size(branch_ratioA, 2)
  
  ! Region B (750-1080 A)
  call read_photoB(sp_list, waveB, & ! input variables
    im_photoB_all, n_branchB, ptitleB, enrgIB, is_ionizationB, charge_stateB, &
    im_photoB_product_all, waveB_lres, crsB, branch_ratioB)
  
  n_waveB_lres = size(waveB_lres)
  n_sp_photoB = size(crsB, 2)
  n_branch_maxB = size(branch_ratioB, 2)
  
  ! Region C (1080-2000 A)
  call read_photo(file_photoC, sp_list, & ! input variables
    im_photoC_all, n_branchC, ptitleC, enrgIC, is_ionizationC, charge_stateC, &
    im_photoC_product_all, waveC, crsC, branch_ratioC)
  
  n_waveC = size(waveC)
  n_sp_photoC = size(crsC, 2)
  n_branch_maxC = size(branch_ratioC, 2)
  
  ! Specified rates of species with no optical effects
  call read_Jvals(file_jvals, sp_list, & ! input variables
    im_photoJ_all, n_branchJ, ptitleJ, im_photoJ_product_all, srateJ)
  
  n_sp_photoJ = size(srateJ,2)
  n_branch_maxJ = size(srateJ,1)

  !----------------------------------------------------------------------------
  !  Read solar spectra, and interpolate onto standard high-resolution 
  !  wavelength scale for region B
  !----------------------------------------------------------------------------
  
  ! --- Region A (0-750 A) ----------------------------------------------------
  call read_sol(file_solA, sol_waveA, sol_fluxA)

  ! --- Region B (750-1080 A) -------------------------------------------------
  allocate(sol_fluxB(n_waveB))
  call read_sol(file_solB, sol_waveB_lres, sol_fluxB_lres)
  
  ! get high-resolution solar fluxes and normalize to low-resolution spectrum
  call read_sol(file_sol_hres, sol_wave_hres, sol_flux_hres)
  call norm_sol(sol_waveB_lres, sol_fluxB_lres, sol_wave_hres, sol_flux_hres)

  ! high-resolution wavelength scale from solar fluxes spans 800-1150 A, so
  ! long-wavelength end goes beyond range of region B
  
  ! find index in standard wavelength scale that corresponds to start of
  ! high-resolution solar spectra
  ti1 = locate(sol_wave_hres(1), waveB)
  ! ensure starting index is in range for interpolation
  if (waveB(ti1) < sol_wave_hres(1)) ti1=ti1+1
  
  ! interpolate onto standard wavelength scale
  call intrp(sol_waveB_lres, sol_fluxB_lres, &
    waveB(1:ti1-1), sol_fluxB(1:ti1-1))
  call intrp(sol_wave_hres, sol_flux_hres, &
    waveB(ti1:n_waveB), sol_fluxB(ti1:n_waveB))
  
  ! unit conversion of high-resolution spectrum to new bin size
  ! 10 A -> 2E-4 A
  do concurrent (i_wave=1:n_waveB)
    sol_fluxB(i_wave) = sol_fluxB(i_wave) * dwaveB(i_wave) / ten
  end do
  
  ! --- Region C (1080-2000 A) ------------------------------------------------
  call read_sol(file_solC, sol_waveC, sol_fluxC)
  
  ! --- Scale to Mars orbital distance (fluxes are 1-AU values) ---------------
  sol_fluxA = sol_fluxA/DH**2
  sol_fluxB = sol_fluxB/DH**2
  sol_fluxC = sol_fluxC/DH**2
  srateJ = srateJ/DH**2

  !----------------------------------------------------------------------------
  !  Combine regions
  !----------------------------------------------------------------------------
  
  ! count number of photo reactions
  nr = 0
  do i_sp = 1, n_sp_photoA
    do i_branch = 1, n_branchA(i_sp)
      nr = nr + 1
    end do
  end do
  do i_sp = 1, n_sp_photoB
    do i_branch = 1, n_branchB(i_sp)
      nr = nr + 1
    end do
  end do
  do i_sp = 1, n_sp_photoC
    do i_branch = 1, n_branchC(i_sp)
      nr = nr + 1
    end do
  end do
  do i_sp = 1, n_sp_photoJ
    do i_branch = 1, n_branchJ(i_sp)
      nr = nr + 1
    end do
  end do
  n_prct = nr
  
  allocate(im_photo_all(5,n_prct), ptitle(n_prct)) 
  nr = 0
  do i_sp = 1, n_sp_photoA
    do i_branch = 1, n_branchA(i_sp)
      nr = nr + 1
      ptitle(nr) = ptitleA(i_branch,i_sp)
      im_photo_all(1,nr) = im_photoA_all(i_sp)
      im_photo_all(2:5,nr) = im_photoA_product_all(1:4,i_branch,i_sp)
    end do
  end do

  do i_sp = 1, n_sp_photoB
    do i_branch = 1, n_branchB(i_sp)
      nr = nr + 1
      ptitle(nr) = ptitleB(i_branch,i_sp)
      im_photo_all(1,nr) = im_photoB_all(i_sp)
      im_photo_all(2:5,nr) = im_photoB_product_all(1:4,i_branch,i_sp)
    end do
  end do

  do i_sp = 1, n_sp_photoC
    do i_branch = 1, n_branchC(i_sp)
      nr = nr + 1
      ptitle(nr) = ptitleC(i_branch,i_sp)
      im_photo_all(1,nr) = im_photoC_all(i_sp)
      im_photo_all(2:5,nr) = im_photoC_product_all(1:4,i_branch,i_sp)
    end do
  end do

  do i_sp = 1, n_sp_photoJ
    do i_branch = 1, n_branchJ(i_sp)
      nr = nr + 1
      ptitle(nr) = ptitleJ(i_branch,i_sp)
      im_photo_all(1,nr) = im_photoJ_all(i_sp)
      im_photo_all(2:5,nr) = im_photoJ_product_all(1:4,i_branch,i_sp)
    end do
  end do
  allocate(trnA(n_waveA,n_z), prtA(n_waveA,n_branch_maxA,n_sp_photoA))
  allocate(trnB(n_waveB,n_z), prtB(n_waveB,n_branch_maxB,n_sp_photoB))
  allocate(trnC(n_waveC,n_z), prtC(n_waveC,n_branch_maxC,n_sp_photoC))
  allocate(rph(n_prct,n_z),rpt(n_prct,n_z))  

  deallocate(dwaveB)
  deallocate(ptitleA, ptitleB, ptitleC, ptitleJ)
  deallocate(im_photoA_product_all, im_photoB_product_all)
  deallocate(im_photoC_product_all, im_photoJ_product_all)
  deallocate(sol_waveA, sol_waveB_lres, sol_waveC, sol_wave_hres)
  deallocate(sol_fluxB_lres, sol_flux_hres)

end subroutine photo


subroutine read_photo(file_photo, sp_list, & ! input variables
  im_photo_all, n_branch, ptitle, enrgI, is_ionization, charge_state, &
  im_photo_product_all, wave, crs, branch_ratio)
! This subroutine reads the cross sections and branching ratios in spectral
! regions A and C for photo reactions. 

  use types, only: wp => dp
  use constants
  use utils, only : find_name

  implicit none

  !----------------------------------------------------------------------------
  !  External variables
  !----------------------------------------------------------------------------
  
  ! input file name
  character(len=*), intent(in) :: file_photo
  ! list of all species
  character(len=*), intent(in), dimension(:) :: sp_list
  
  ! index mapping from list of photo species -> list of all species
  integer, intent(out), allocatable, dimension(:) :: im_photo_all
  ! number of branches
  integer, intent(out), allocatable, dimension(:) :: n_branch
  ! equation for photo reaction
  character(len=*), intent(out), allocatable, dimension(:,:) :: ptitle
  ! threshold energy in angstroms
  real(wp), intent(out), allocatable, dimension(:,:) :: enrgI
  ! is reaction an ionization reaction?
  logical, intent(out), allocatable, dimension(:,:) :: is_ionization
  ! number of ions produced in reaction
  real(wp), intent(out), allocatable, dimension(:,:) :: charge_state
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, intent(out), allocatable, dimension(:,:,:) :: im_photo_product_all
  ! wavelength scale for cross sections
  real(wp), intent(out), allocatable, dimension(:) :: wave
  ! cross sections
  real(wp), intent(out), allocatable, dimension(:,:) :: crs
  ! branch ratio
  real(wp), intent(out), allocatable, dimension(:,:,:) :: branch_ratio

  !----------------------------------------------------------------------------
  !  Internal variables
  !----------------------------------------------------------------------------
  
  ! file unit number for input file
  integer :: fid
  ! total number of species in model
  integer :: n_sp
  ! number of species with listed cross sections
  integer :: n_sp_photo
  ! number of wavelength bins
  integer :: n_wave
  ! maximum branches for each species with cross sections
  integer :: n_branch_max
  ! maximum number of products in reaction formula
  integer, parameter :: n_product_max = 4

  ! loop variables
  integer :: i_sp, i_branch, i_wave, i, np
  ! temporary / dummy variables
  character(len=256) :: ts_line
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  real(wp), allocatable, dimension(:) :: t_crs
  integer :: tn_isp
  
  !----------------------------------------------------------------------------
  !  Read reactions, cross sections and branching ratios
  !----------------------------------------------------------------------------
  
  n_sp = size(sp_list)-1 ! since 0th index = total

  open(newunit=fid, file=file_photo, status='old', action='read')

  read(fid,*) n_wave, n_sp_photo, n_branch_max

  allocate(wave(n_wave), crs(n_wave,n_sp_photo), n_branch(n_sp_photo))
  allocate(ptitle(n_branch_max,n_sp_photo))
  allocate(enrgI(n_branch_max,n_sp_photo))
  allocate(is_ionization(n_branch_max,n_sp_photo))
  allocate(charge_state(n_branch_max,n_sp_photo))
  allocate(im_photo_all(n_sp_photo))
  allocate(im_photo_product_all(n_product_max,n_branch_max,n_sp_photo))
  allocate(branch_ratio(n_wave,n_branch_max,n_sp_photo))

  allocate(t_crs(n_wave))
  
  ! initialize variables
  enrgI = zero
  is_ionization = .false.
  charge_state = zero
  im_photo_all = 0
  im_photo_product_all = 0
  
  ! read wavelength scale
  read(fid,'(A)') ts_line
  read(fid,'(10ES11.3)') wave
  
  do i_sp = 1, n_sp_photo
    read(fid,'(A12, I4)') ts_species1, n_branch(i_sp) ! header for species
    im_photo_all(i_sp) = find_name(ts_species1, sp_list)
    
    read(fid,'(A)') ts_line ! "Total Absorption Cross Section (cm^2)"
    read(fid,*) (crs(i_wave,i_sp), i_wave=1,n_wave)
    
    read(fid,'(A)') ts_line ! "Total Ionization Cross Section (cm^2)"
    read(fid,*) (t_crs(i_wave), i_wave=1,n_wave)
    
    do i_branch = 1, n_branch(i_sp)
      ! read reaction formula
      read (fid,'(A)') ts_line
      ptitle(i_branch,i_sp) = ts_line
      read (ts_line,"(5(A12,3X), A12, F7.1)") &
        (ts_species2(i), i=1,6), enrgI(i_branch,i_sp)
      
      ! link products to species list
      np = 0
      do i = 3, 6 ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          np = np + 1
          im_photo_product_all(np,i_branch,i_sp) = tn_isp
        end if
        ! determine if ionization reaction (is electron produced?)
        if (tn_isp == n_sp) then
          is_ionization(i_branch,i_sp) = .true.
          charge_state(i_branch,i_sp) = charge_state(i_branch,i_sp) + one
        end if
      end do
      ! read branching ratios
      read(fid,'(10ES11.3)') &
        (branch_ratio(i_wave,i_branch,i_sp), i_wave=1,n_wave)
    end do
  end do
  
  close(unit=fid)
  deallocate(t_crs)

end subroutine read_photo


subroutine read_photoB(sp_list, wave, & ! input variables
  im_photo_all, n_branch, ptitle, enrgI, is_ionization, charge_state, &
  im_photo_product_all, wave_lres, crs, branch_ratio)
! This subroutine reads the cross sections and branching ratios in spectral
! region B for photo reactions. Most cross sections use the low-resolution
! wavelength scale, but N2 and CO uses a high-resolution scale to account for
! the pre-dissociating states. 

  use types, only: wp => dp
  use constants
  use utils, only : find_name, intrp, locate

  implicit none

  !----------------------------------------------------------------------------
  !  External variables
  !----------------------------------------------------------------------------
  
  ! list of all species
  character(len=*), intent(in), dimension(:) :: sp_list
  ! wavelength scale (high resolution) for cross sections
  real(wp), intent(in), dimension(:) :: wave
  
  ! index mapping from list of photo species -> list of all species
  integer, intent(out), allocatable, dimension(:) :: im_photo_all
  ! number of branches
  integer, intent(out), allocatable, dimension(:) :: n_branch
  ! equation for photo reaction
  character(len=*), intent(out), allocatable, dimension(:,:) :: ptitle
  ! threshold energy in angstroms
  real(wp), intent(out), allocatable, dimension(:,:) :: enrgI
  ! is reaction an ionization reaction?
  logical, intent(out), allocatable, dimension(:,:) :: is_ionization
  ! number of ions produced in reaction
  real(wp), intent(out), allocatable, dimension(:,:) :: charge_state
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, intent(out), allocatable, dimension(:,:,:) :: im_photo_product_all
  ! wavelength scale (low resolution) for cross sections
  real(wp), intent(out), allocatable, dimension(:) :: wave_lres
  ! cross sections
  real(wp), intent(out), allocatable, dimension(:,:) :: crs
  ! branch ratio
  real(wp), intent(out), allocatable, dimension(:,:,:) :: branch_ratio
  
  !----------------------------------------------------------------------------
  !  Internal variables
  !----------------------------------------------------------------------------
  
  ! total number of species in model
  integer :: n_sp
  ! number of wavelength bins
  integer :: n_wave
  ! number of species with listed cross sections
  integer :: n_sp_photo
  ! maximum branches for each species with cross sections
  integer :: n_branch_max
  ! maximum number of products in reaction formula
  integer, parameter :: n_product_max = 4
  
  ! --- Main cross section file -----------------------------------------------
  ! file unit number
  integer :: fid
  ! number of wavelength bins
  integer :: n_wave_lres
  ! number of species with listed cross sections
  integer :: n_sp_photo_lres
  ! cross sections (low resolution)
  real(wp), allocatable, dimension(:) :: crs_lres
  ! branch ratios (low resolution)
  real(wp), allocatable, dimension(:) :: branch_ratio_lres

  ! --- N2 cross section files ------------------------------------------------
  ! file unit numbers
  integer :: fid_28N2, fid_29N2
  ! number of wavelength bins
  integer :: n_wave_N2
  ! wavelength scale for cross sections
  real(wp), allocatable, dimension(:) :: wave_N2
  ! cross sections
  real(wp), allocatable, dimension(:) :: crs_N2

  ! --- CO cross section files (a = 75-91.2 nm; b = 91.2-108 nm) --------------
  ! file unit numbers
  integer :: fid_COa, fid_COb
  ! number of wavelength bins
  integer :: n_wave_COa, n_wave_COb
  ! wavelength scale for cross sections
  real(wp), allocatable, dimension(:) :: wave_CO
  ! cross sections (total, dissociation)
  real(wp), allocatable, dimension(:) ::  crs_CO, crs_CO_diss
  ! branch ratios
  real(wp), allocatable, dimension(:,:) :: branch_ratio_CO
    
  ! --- Loop variables --------------------------------------------------------
  integer :: i_sp, i_branch, i_wave, i_wavenum, i, np
  ! --- Temporary / dummy variables -------------------------------------------
  character(len=256) :: ts_line
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  real(wp), allocatable, dimension(:) :: t_crs
  integer :: tn_isp
  
  !----------------------------------------------------------------------------
  !  Read reactions, cross sections and branching ratios from main cross
  !  section file
  !----------------------------------------------------------------------------

  n_sp = size(sp_list)-1 ! since 0th index = total
  n_wave = size(wave)
  
  open(newunit=fid, file='../data/photons/photoB.dat', &
    status='old', action='read')
    
  read(fid,*) n_wave_lres, n_sp_photo_lres, n_branch_max
  n_sp_photo = n_sp_photo_lres + 3
  
  allocate(crs(n_wave,n_sp_photo), n_branch(n_sp_photo))
  allocate(ptitle(n_branch_max,n_sp_photo))
  allocate(enrgI(n_branch_max,n_sp_photo))
  allocate(is_ionization(n_branch_max,n_sp_photo))
  allocate(charge_state(n_branch_max,n_sp_photo))
  allocate(im_photo_all(n_sp_photo))
  allocate(im_photo_product_all(n_product_max,n_branch_max,n_sp_photo))
  allocate(branch_ratio(n_wave,n_branch_max,n_sp_photo))

  ! initialize variables
  enrgI = zero
  is_ionization = .false.
  charge_state = 0
  im_photo_all = 0
  im_photo_product_all = 0

  allocate(wave_lres(n_wave_lres), crs_lres(n_wave_lres))
  allocate(branch_ratio_lres(n_wave_lres))
  allocate(t_crs(n_wave_lres))
  
  ! read wavelength scale
  read(fid,'(A)') ts_line
  read(fid,'(10ES11.3)') wave_lres
  
  do i_sp = 1, n_sp_photo_lres
    read(fid,'(A12, I4)') ts_species1, n_branch(i_sp) ! header for species
    im_photo_all(i_sp) = find_name(ts_species1, sp_list)
    
    read(fid,'(A)') ts_line ! "Total Absorption Cross Section (cm^2)"
    read(fid,*) (crs_lres(i_wave), i_wave=1,n_wave_lres)
    ! interpolate cross sections onto high resolution scale
    call intrp(wave_lres, crs_lres, wave, crs(:,i_sp))
    
    read(fid,'(A)') ts_line ! "Total Ion Cross Section (cm^2)"
    read(fid,*) (t_crs(i_wave), i_wave=1,n_wave_lres)
    
    do i_branch = 1, n_branch(i_sp)
      ! read reaction formula
      read (fid,'(A)') ts_line
      ptitle(i_branch,i_sp) = ts_line
      read (ts_line,"(5(A12,3X), A12, F7.1)") &
        (ts_species2(i),i=1,6), enrgI(i_branch,i_sp)

      ! link products to species list
      np = 0
      do i = 3, 6 ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          np = np + 1
          im_photo_product_all(np,i_branch,i_sp) = tn_isp
        end if
        ! determine if ionization reaction (is electron produced?)
        if (tn_isp == n_sp) then
          is_ionization(i_branch,i_sp) = .true.
          charge_state(i_branch,i_sp) = charge_state(i_branch,i_sp) + one
        end if
      end do
      ! read branching ratios
      read(fid,'(10ES11.3)') (branch_ratio_lres(i_wave), i_wave=1,n_wave_lres)
      ! interpolate branching ratios onto high resolution scale
      call intrp(wave_lres, branch_ratio_lres, &
        wave, branch_ratio(:,i_branch,i_sp))
    end do
  end do
     
  close(unit=fid)
  deallocate(crs_lres, branch_ratio_lres, t_crs)
  
  !----------------------------------------------------------------------------
  !  Read N2 cross sections over 800-1000 Angstroms from Lewis et al. Then
  !  interpolate to high-resolution wavelength grid. 
  !----------------------------------------------------------------------------

  n_wave_N2 = 540000
  allocate(wave_N2(n_wave_N2), crs_N2(n_wave_N2))
  
  ! --- 28N2 ------------------------------------------------------------------
  i_sp = n_sp_photo_lres + 1
  open(newunit=fid_28N2, file='../data/photons/photoB-28N2.dat', &
    status='old', action='read')
  read(fid_28N2,"(A)") ts_line ! "## Processed by..."
  read(fid_28N2,"(A)") ts_line ! "## CROSS_SECTION=..."
  read(fid_28N2,"(A)") ts_line ! "Wavenumber(cm-1) Cross-section(cm2)..."
  ! convert wavenumber scale (cm-1) to wavelength (angstrom) scale
  do i_wavenum = 1, n_wave_N2
    ! need to flip indices since wavelength is reciprocal of wavenumber
    i_wave = n_wave_N2 - i_wavenum + 1
    read(fid_28N2,*) wave_N2(i_wave), crs_N2(i_wave)
    ! convert cm-1 to angstrom
    wave_N2(i_wave) = 1.E8_wp / wave_N2(i_wave)
  end do
  close(unit=fid_28N2)

  im_photo_all(i_sp) = find_name('N2          ', sp_list)
  n_branch(i_sp) = 1
  i_branch = 1
  ptitle(i_branch,i_sp) = 'N2           + hv           = '// &
    'N2D          + N            +              +             '
  enrgI(i_branch,i_sp) = zero
  is_ionization(i_branch,i_sp) = .false.
  charge_state(i_branch,i_sp) = zero
  branch_ratio(:,i_branch,i_sp) = one
  im_photo_product_all(1,i_branch,i_sp) = find_name('N2D         ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('N           ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0
  
  ! interpolate cross sections onto high resolution scale
  call intrp(wave_N2, crs_N2, wave, crs(:,i_sp))

  ! --- 29N2 ------------------------------------------------------------------
  i_sp = n_sp_photo_lres + 2
  open(newunit=fid_29N2, file='../data/photons/photoB-29N2.dat', &
    status='old', action='read')
  read(fid_29N2,"(A)") ts_line ! "## Processed by..."
  read(fid_29N2,"(A)") ts_line ! "## CROSS_SECTION=..."
  read(fid_29N2,"(A)") ts_line ! "Wavenumber(cm-1) Cross-section(cm2)..."
  ! convert wavenumber scale (cm-1) to wavelength (angstrom) scale
  do i_wavenum = 1, n_wave_N2
    ! need to flip indices since wavelength is reciprocal of wavenumber
    i_wave = n_wave_N2 - i_wavenum + 1
    read(fid_29N2,*) wave_N2(i_wave), crs_N2(i_wave)
    ! convert cm-1 to angstrom
    wave_N2(i_wave) = 1.E8_wp / wave_N2(i_wave)
  end do
  close(unit=fid_29N2)

  im_photo_all(i_sp) = find_name('N2I         ', sp_list)
  n_branch(i_sp) = 2
  
  i_branch = 1
  ptitle(i_branch,i_sp) = 'N2I          + hv           = '// &
    'N2DI         + N            +              +             '
  enrgI(i_branch,i_sp) = zero
  is_ionization(i_branch,i_sp) = .false.
  charge_state(i_branch,i_sp) = zero
  branch_ratio(:,i_branch,i_sp) = half
  im_photo_product_all(1,i_branch,i_sp) = find_name('N2DI        ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('N           ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0

  i_branch = 2
  ptitle(i_branch,i_sp) = 'N2I          + hv           = '// &
    'N2D          + NI           +              +             '
  enrgI(i_branch,i_sp) = zero
  is_ionization(i_branch,i_sp) = .false.
  charge_state(i_branch,i_sp) = zero
  branch_ratio(:,i_branch,i_sp) = half
  im_photo_product_all(1,i_branch,i_sp) = find_name('N2D         ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('NI          ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0

  ! interpolate cross sections onto high resolution scale
  call intrp(wave_N2, crs_N2, wave, crs(:,i_sp))
  deallocate(wave_N2, crs_N2)
  
  !----------------------------------------------------------------------------
  !  Read CO cross sections from Alan Heays. The first file spans 75-91.2 nm,
  !  and the second file spans 91.2-108 nm. Then interpolate to
  !  high-resolution wavelength grid. 
  !----------------------------------------------------------------------------

  n_wave_COa = 162000
  n_wave_COb = 840000

  i_sp = n_sp_photo_lres + 3
  n_branch(i_sp) = 2
  
  allocate(wave_CO(n_wave_COa+n_wave_COb))
  allocate(crs_CO(n_wave_COa+n_wave_COb))
  allocate(crs_CO_diss(n_wave_COa+n_wave_COb))
  allocate(branch_ratio_CO(n_wave_COa+n_wave_COb,n_branch(i_sp)))
  
  ! --- 75-91.2 nm ------------------------------------------------------------
  open(newunit=fid_COa, &
    file='../data/photons/photoB-12C16O_300K_75-91.2.dat', &
    status='old', action='read')
  read(fid_COa,"(A)") ts_line ! "# 12C16O photoabsorption..."
  read(fid_COa,"(A)") ts_line ! "# wavelength(nm)..."
  do i_wave = 1, n_wave_COa
    read(fid_COa,*) wave_CO(i_wave), &
      crs_CO(i_wave), crs_CO_diss(i_wave)
    ! convert nm to angstroms
    wave_CO(i_wave) = wave_CO(i_wave) * ten
  end do
  close(unit=fid_COa)
  
  ! --- 91.2-108 nm -----------------------------------------------------------
  open(newunit=fid_COb, &
    file='../data/photons/photoB-12C16O_300K_91.2-108.dat', &
    status='old', action='read')
  read(fid_COb,"(A)") ts_line ! "# Photoabsorption..."
  read(fid_COb,"(A)") ts_line ! "# Based on the..."
  read(fid_COb,"(A)") ts_line ! "# Assembled by..."
  read(fid_COb,"(A)") ts_line ! "# wavelength(nm)..."
  do i_wave = n_wave_COa+1, n_wave_COa+n_wave_COb
    read(fid_COb,*) wave_CO(i_wave), &
      crs_CO(i_wave), crs_CO_diss(i_wave)
    ! convert nm to angstroms
    wave_CO(i_wave) = wave_CO(i_wave) * ten
  end do
  close(unit=fid_COb)
  
  im_photo_all(i_sp) = find_name('CO          ', sp_list)
  
  ! interpolate cross sections onto high resolution scale
  call intrp(wave_CO, crs_CO, wave, crs(:,i_sp))
  
  ! --- Photodissociation branch ----------------------------------------------
  i_branch = 1
  ptitle(i_branch,i_sp) = 'CO           + hv           = '// &
    'C            + O            +              +             '
  enrgI(i_branch,i_sp) = zero
  is_ionization(i_branch,i_sp) = .false.
  charge_state(i_branch,i_sp) = zero
  do concurrent (i_wave = 1:n_wave_COa+n_wave_COb)
    branch_ratio_CO(i_wave,i_branch) = crs_CO_diss(i_wave) / crs_CO(i_wave)
  end do
  im_photo_product_all(1,i_branch,i_sp) = find_name('C           ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('O           ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0
  
  ! interpolate branching ratios onto high resolution scale
  call intrp(wave_CO, branch_ratio_CO(:,i_branch), &
    wave, branch_ratio(:,i_branch,i_sp))
  
  ! --- Photoionization branch ------------------------------------------------
  i_branch = 2
  ptitle(i_branch,i_sp) = 'CO           + hv           = '// &
    'COP          + E            +              +             '
  enrgI(i_branch,i_sp) = 884.7_wp
  is_ionization(i_branch,i_sp) = .true.
  charge_state(i_branch,i_sp) = one
  do concurrent (i_wave = 1:n_wave_COa+n_wave_COb)
     branch_ratio_CO(i_wave,i_branch) = one - branch_ratio_CO(i_wave,1)
  end do
  im_photo_product_all(1,i_branch,i_sp) = find_name('COP         ',sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('E           ',sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0

  ! interpolate branching ratios onto high resolution scale
  call intrp(wave_CO, branch_ratio_CO(:,i_branch), &
    wave, branch_ratio(:,i_branch,i_sp))
  
  deallocate(wave_CO, crs_CO, crs_co_diss, branch_ratio_co)
  
end subroutine read_photoB


subroutine read_Jvals(file_jvals, sp_list, &
  im_photo_all, n_branch, ptitle, im_photo_product_all, srateJ)
! This subroutine reads in specified rates from species whose optical effects
! are not calculated in the model. 

  use types, only: wp => dp
  use constants
  use utils, only : find_name

  implicit none

  !----------------------------------------------------------------------------
  !  External variables
  !----------------------------------------------------------------------------

  ! input file name
  character(len=*), intent(in) :: file_jvals
  ! list of all species
  character(len=*), intent(in), dimension(:) :: sp_list
  
  ! index mapping from list of photo species -> list of all species
  integer, intent(out), allocatable, dimension(:) :: im_photo_all
  ! number of branches
  integer, intent(out), allocatable, dimension(:) :: n_branch
  ! equation for photo reaction
  character(len=*), intent(out), allocatable, dimension(:,:) :: ptitle
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, intent(out), allocatable, dimension(:,:,:) :: im_photo_product_all
  ! specified reaction rate (scaled to 1 AU)
  real(wp), intent(out), allocatable, dimension(:,:) :: srateJ

  !----------------------------------------------------------------------------
  !  Internal variables
  !----------------------------------------------------------------------------
  
  ! file unit number for input file
  integer :: fid
  ! total number of species in model
  integer :: n_sp
  ! number of species with listed cross sections
  integer :: n_sp_photo
  ! maximum branches for each species with cross sections
  integer :: n_branch_max
  ! maximum number of products in reaction formula
  integer, parameter :: n_product_max = 4
  
  ! loop variables
  integer :: i_sp, i_branch, i, np
  ! temporary / dummy variables
  character(len=256) :: ts_line
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  integer :: tn_isp
  
  !----------------------------------------------------------------------------
  !  Read in specified reactions and rates
  !----------------------------------------------------------------------------

  n_sp = size(sp_list)

  open(newunit=fid, file=file_jvals, status='old', action='read')
  
  read(fid,*) n_sp_photo, n_branch_max
  
  allocate(n_branch(n_sp_photo))
  allocate(ptitle(n_branch_max,n_sp_photo))
  allocate(im_photo_all(n_sp_photo))
  allocate(im_photo_product_all(n_product_max,n_branch_max,n_sp_photo))
  allocate(srateJ(n_branch_max,n_sp_photo))
  
  ! initialize variables
  im_photo_all = 0
  im_photo_product_all = 0
  
  do i_sp = 1, n_sp_photo
    read(fid,'(A12, I4)') ts_species1, n_branch(i_sp) ! header for species
    im_photo_all(i_sp) = find_name(ts_species1, sp_list)
    
    do i_branch = 1, n_branch(i_sp)
      ! read reaction formula
      read (fid,'(A)') ts_line
      ptitle(i_branch,i_sp) = ts_line
      read (ts_line,"(5(A12,3X), A12, ES11.3)") &
        (ts_species2(i), i=1,6), srateJ(i_branch,i_sp)
      
      ! link products to species list
      np = 0
      do i = 3, 6 ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          np = np + 1
          im_photo_product_all(np,i_branch,i_sp) = tn_isp
        end if
      end do
    end do
  end do
    
  close(unit=fid)

end subroutine read_Jvals


subroutine read_sol(file_sol, & ! input variable
  sol_wave, sol_flux)
! This subroutine reads in the solar fluxes in /solar/.
  use types, only: wp => dp
  use constants
  
  implicit none
  
  !----------------------------------------------------------------------------
  !  External variables
  !----------------------------------------------------------------------------
  
  character(len=*), intent(in) :: file_sol
  real(wp), intent(out), allocatable, dimension(:) :: sol_wave
  real(wp), intent(out), allocatable, dimension(:) :: sol_flux
  
  !----------------------------------------------------------------------------
  !  Internal variables
  !----------------------------------------------------------------------------
  
  ! file unit number for input file
  integer :: fid
  ! number of wavelength bins
  integer :: n_wave
  
  ! loop variables
  integer :: i_wave
  ! temporary / dummy variables
  character(len=128) :: ts_line

  !----------------------------------------------------------------------------
  !  Read solar fluxes
  !----------------------------------------------------------------------------
  
  open(newunit=fid, file=file_sol, status='old', action='read')
  
  read(fid,*) n_wave
  allocate(sol_wave(n_wave), sol_flux(n_wave))
  
  read(fid,"(A)") ts_line ! e.g., "; High solar activity reference..."
  read(fid,"(A)") ts_line ! "; Values are..."
  read(fid,"(A)") ts_line ! "; wave (A)..."
  
  do i_wave = 1, n_wave
    read(fid,*) sol_wave(i_wave), sol_flux(i_wave)
  end do
  
  close(unit=fid)

end subroutine read_sol


subroutine norm_sol(sol_wave_lres, sol_flux_lres, sol_wave_hres, sol_flux_hres)
! This subroutine normalizes the high-resolution solar fluxes in sol_flux_hres
! over wavelength scale sol_wave_hres against the low-resolution fluxes
! sol_wave_lres, which has a wavelength scale sol_flux_lres. 

  use types, only: wp => dp
  use constants
  use utils, only : locate
  
  implicit none

  !----------------------------------------------------------------------------
  !  External variables
  !----------------------------------------------------------------------------
  
  real(wp), intent(in), dimension(:) :: sol_wave_lres
  real(wp), intent(in), dimension(:) :: sol_flux_lres
  real(wp), intent(in), dimension(:) :: sol_wave_hres
  real(wp), intent(inout), dimension(:) :: sol_flux_hres
   
  !----------------------------------------------------------------------------
  !  Internal variables
  !----------------------------------------------------------------------------
  ! number of wavelength bins in low and high resolution
  integer :: n_sol_wave_lres, n_sol_wave_hres
  ! wavelength bin widths in low and high resolution
  real(wp), allocatable, dimension(:) :: dwave_lres, dwave_hres

  ! loop variables
  integer :: i_wave
  ! temporary / dummy variables
  integer :: ti_low, ti_high
  real(wp) :: t_bin_low, t_bin_high, t_flux, t_ratio

  !----------------------------------------------------------------------------
  ! Initialize wavelength scale variables
  !----------------------------------------------------------------------------
  
  n_sol_wave_lres = size(sol_wave_lres)
  allocate(dwave_lres(n_sol_wave_lres))
  do concurrent (i_wave = 2:n_sol_wave_lres-1)
    dwave_lres(i_wave) = half * &
    (sol_wave_lres(i_wave+1) - sol_wave_lres(i_wave-1))
  end do
  dwave_lres(1) = sol_wave_lres(2) - sol_wave_lres(1)
  dwave_lres(n_sol_wave_lres) = &
    sol_wave_lres(n_sol_wave_lres) - sol_wave_lres(n_sol_wave_lres-1)
  
  n_sol_wave_hres = size(sol_wave_hres)
  allocate(dwave_hres(n_sol_wave_hres))
  do concurrent (i_wave = 2:n_sol_wave_hres-1)
    dwave_hres(i_wave) = half * &
    (sol_wave_hres(i_wave+1) - sol_wave_hres(i_wave-1))
  end do
  dwave_hres(1) = sol_wave_hres(2) - sol_wave_hres(1)
  dwave_hres(n_sol_wave_hres) = &
    sol_wave_hres(n_sol_wave_hres) - sol_wave_hres(n_sol_wave_hres-1)
  
  !----------------------------------------------------------------------------
  ! Scale high-resolution solar fluxes to low-resolution solar fluxes
  !----------------------------------------------------------------------------
  
  do concurrent (i_wave = 1:n_sol_wave_lres)
    ! determine wavelength of low-wavelength bin boundaries
    t_bin_low = sol_wave_lres(i_wave) - half * dwave_lres(i_wave)
    t_bin_high = sol_wave_lres(i_wave) + half * dwave_lres(i_wave)
    
    ! check if current low-resolution bin has any high-resolution fluxes
    if((t_bin_low >= sol_wave_hres(1)) .and. &
      (t_bin_high <= sol_wave_hres(n_sol_wave_hres))) then
      
      ! determine indices of high-resolution wavelength scale corresponding to
      ! bin boundaries
      ti_low = locate(t_bin_low, sol_wave_hres)
      ti_high = locate(t_bin_high, sol_wave_hres)-1 
      ! -1 so the high end bin would not be doubly normalized when it becomes
      ! low end bin in the next iteration

      ! integrate high-resolution flux over interval spanned by low-resolution
      ! bin
      t_flux = dot_product(sol_flux_hres(ti_low:ti_high), &
        dwave_hres(ti_low:ti_high))
      
      t_ratio = sol_flux_lres(i_wave) * dwave_lres(i_wave) / t_flux
      sol_flux_hres(ti_low:ti_high) = t_ratio * sol_flux_hres(ti_low:ti_high)
      
    end if
  end do

  deallocate(dwave_lres, dwave_hres)
end subroutine norm_sol

end module
