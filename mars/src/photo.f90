module photo

contains

subroutine photo_setup
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
  character(len=128) :: file_solA, file_solB, file_solC, file_sol_hres
  
  ! high resolution wavelength bin size for region B
  real(wp), allocatable, dimension(:) :: dwaveB
  ! maximum number of branches for each species
  integer :: n_branch_maxA, n_branch_maxB, n_branch_maxC, n_branch_maxJ
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
  integer :: i_sp, i_branch, i_wave, i_r
  ! temporary / dummy variables
  integer :: ti1
  
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
    im_photoA_all, n_branchA, ptitleA, enrgIA, charge_stateA, & 
    im_photoA_product_all, waveA, csA, branch_ratioA)
  
  n_waveA = size(waveA)
  n_sp_photoA = size(csA, 2)
  n_branch_maxA = size(branch_ratioA, 2)
  
  ! Region B (750-1080 A)
  call read_photoB(sp_list, waveB, & ! input variables
    im_photoB_all, n_branchB, ptitleB, enrgIB, charge_stateB, &
    im_photoB_product_all, waveB_lres, csB, branch_ratioB)
  
  n_waveB_lres = size(waveB_lres)
  n_sp_photoB = size(csB, 2)
  n_branch_maxB = size(branch_ratioB, 2)
  
  ! Region C (1080-2000 A)
  call read_photo(file_photoC, sp_list, & ! input variables
    im_photoC_all, n_branchC, ptitleC, enrgIC, charge_stateC, &
    im_photoC_product_all, waveC, csC, branch_ratioC)
  
  n_waveC = size(waveC)
  n_sp_photoC = size(csC, 2)
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
  do concurrent (i_wave = 1:n_waveB)
    sol_fluxB(i_wave) = sol_fluxB(i_wave) * dwaveB(i_wave) / ten
  end do
  
  ! --- Region C (1080-2000 A) ------------------------------------------------
  call read_sol(file_solC, sol_waveC, sol_fluxC)
  
  ! --- Scale to Mars orbital distance (fluxes are 1-AU values) ---------------
  sol_fluxA = sol_fluxA / DH**2
  sol_fluxB = sol_fluxB / DH**2
  sol_fluxC = sol_fluxC / DH**2
  srateJ = srateJ / DH**2

  !----------------------------------------------------------------------------
  !  Combine regions
  !----------------------------------------------------------------------------
  
  ! count number of photo reactions
  i_r = 0
  do i_sp = 1, n_sp_photoA
    do i_branch = 1, n_branchA(i_sp)
      i_r = i_r + 1
    end do
  end do
  do i_sp = 1, n_sp_photoB
    do i_branch = 1, n_branchB(i_sp)
      i_r = i_r + 1
    end do
  end do
  do i_sp = 1, n_sp_photoC
    do i_branch = 1, n_branchC(i_sp)
      i_r = i_r + 1
    end do
  end do
  do i_sp = 1, n_sp_photoJ
    do i_branch = 1, n_branchJ(i_sp)
      i_r = i_r + 1
    end do
  end do
  n_prct = i_r
  
  allocate(ptitle(n_prct), im_photo_all(5,n_prct)) 
  i_r = 0
  do i_sp = 1, n_sp_photoA
    do i_branch = 1, n_branchA(i_sp)
      i_r = i_r + 1
      ptitle(i_r) = ptitleA(i_branch,i_sp)
      im_photo_all(1,i_r) = im_photoA_all(i_sp)
      im_photo_all(2:5,i_r) = im_photoA_product_all(1:4,i_branch,i_sp)
    end do
  end do

  do i_sp = 1, n_sp_photoB
    do i_branch = 1, n_branchB(i_sp)
      i_r = i_r + 1
      ptitle(i_r) = ptitleB(i_branch,i_sp)
      im_photo_all(1,i_r) = im_photoB_all(i_sp)
      im_photo_all(2:5,i_r) = im_photoB_product_all(1:4,i_branch,i_sp)
    end do
  end do

  do i_sp = 1, n_sp_photoC
    do i_branch = 1, n_branchC(i_sp)
      i_r = i_r + 1
      ptitle(i_r) = ptitleC(i_branch,i_sp)
      im_photo_all(1,i_r) = im_photoC_all(i_sp)
      im_photo_all(2:5,i_r) = im_photoC_product_all(1:4,i_branch,i_sp)
    end do
  end do

  do i_sp = 1, n_sp_photoJ
    do i_branch = 1, n_branchJ(i_sp)
      i_r = i_r + 1
      ptitle(i_r) = ptitleJ(i_branch,i_sp)
      im_photo_all(1,i_r) = im_photoJ_all(i_sp)
      im_photo_all(2:5,i_r) = im_photoJ_product_all(1:4,i_branch,i_sp)
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
  
  ! calculate Rayleigh optical depth
  call read_rayleigh

  ! calculate path lengths for solar absorption
  call paths1D

end subroutine photo_setup


subroutine read_photo(file_photo, sp_list, & ! input variables
  im_photo_all, n_branch, ptitle, enrgI, charge_state, &
  im_photo_product_all, wave, cs, branch_ratio)
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
  ! number of ions produced in reaction
  real(wp), intent(out), allocatable, dimension(:,:) :: charge_state
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, intent(out), allocatable, dimension(:,:,:) :: im_photo_product_all
  ! wavelength scale for cross sections
  real(wp), intent(out), allocatable, dimension(:) :: wave
  ! cross sections
  real(wp), intent(out), allocatable, dimension(:,:) :: cs
  ! branch ratio
  real(wp), intent(out), allocatable, dimension(:,:,:) :: branch_ratio

  !----------------------------------------------------------------------------
  !  Local variables
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
  integer :: i_sp, i_branch, i_wave, i, i_p
  ! temporary / dummy variables
  character(len=256) :: ts_line
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  real(wp), allocatable, dimension(:) :: t_cs
  integer :: tn_isp
  
  !----------------------------------------------------------------------------
  !  Read reactions, cross sections and branching ratios
  !----------------------------------------------------------------------------
  
  n_sp = size(sp_list)-1 ! since 0th index = total

  open(newunit=fid, file=file_photo, status='old', action='read')

  read(fid,*) n_wave, n_sp_photo, n_branch_max

  allocate(wave(n_wave), cs(n_wave,n_sp_photo), n_branch(n_sp_photo))
  allocate(ptitle(n_branch_max,n_sp_photo))
  allocate(enrgI(n_branch_max,n_sp_photo))
  allocate(charge_state(n_branch_max,n_sp_photo))
  allocate(im_photo_all(n_sp_photo))
  allocate(im_photo_product_all(n_product_max,n_branch_max,n_sp_photo))
  allocate(branch_ratio(n_wave,n_branch_max,n_sp_photo))

  allocate(t_cs(n_wave))
  
  ! initialize variables
  enrgI = zero
  charge_state = zero
  im_photo_all = 0
  im_photo_product_all = 0
  
  ! read wavelength scale
  read(fid,'(A)') ts_line
  read(fid,'(10ES11.3)') wave
  
  do i_sp = 1, n_sp_photo
    read(fid,'(A12, I4)') ts_species1, n_branch(i_sp) ! header for species
    im_photo_all(i_sp) = find_name(ts_species1, sp_list)
    
    if (im_photo_all(i_sp) <= 0) then
      write(*,'("Error: Photo species not found: ",A12)') ts_species1
      stop
    end if
    
    read(fid,'(A)') ts_line ! "Total Absorption Cross Section (cm^2)"
    read(fid,*) (cs(i_wave,i_sp), i_wave=1,n_wave)
    
    read(fid,'(A)') ts_line ! "Total Ionization Cross Section (cm^2)"
    read(fid,*) (t_cs(i_wave), i_wave=1,n_wave)
    
    do i_branch = 1, n_branch(i_sp)
      ! read reaction formula
      read(fid,'(A)') ts_line
      ptitle(i_branch,i_sp) = ts_line
      read(ts_line,'(5(A12, 3X), A12, F7.1)') &
        (ts_species2(i), i=1,6), enrgI(i_branch,i_sp)
      
      ! link products to species list
      i_p = 0
      do i = 3, 6 ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          i_p = i_p + 1
          im_photo_product_all(i_p,i_branch,i_sp) = tn_isp
        end if
        ! determine if ionization reaction (is electron produced?)
        if (tn_isp == n_sp) then
          charge_state(i_branch,i_sp) = charge_state(i_branch,i_sp) + one
        end if
      end do
      ! read branching ratios
      read(fid,'(10ES11.3)') &
        (branch_ratio(i_wave,i_branch,i_sp), i_wave=1,n_wave)
    end do
  end do
  
  close(unit=fid)
  deallocate(t_cs)

end subroutine read_photo


subroutine read_photoB(sp_list, wave, & ! input variables
  im_photo_all, n_branch, ptitle, enrgI, charge_state, &
  im_photo_product_all, wave_lres, cs, branch_ratio)
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
  ! number of ions produced in reaction
  real(wp), intent(out), allocatable, dimension(:,:) :: charge_state
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, intent(out), allocatable, dimension(:,:,:) :: im_photo_product_all
  ! wavelength scale (low resolution) for cross sections
  real(wp), intent(out), allocatable, dimension(:) :: wave_lres
  ! cross sections
  real(wp), intent(out), allocatable, dimension(:,:) :: cs
  ! branch ratio
  real(wp), intent(out), allocatable, dimension(:,:,:) :: branch_ratio
  
  !----------------------------------------------------------------------------
  !  Local variables
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
  integer :: i_sp, i_branch, i_wave, i_wavenum, i, i_p
  ! --- Temporary / dummy variables -------------------------------------------
  character(len=256) :: ts_line
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  real(wp), allocatable, dimension(:) :: t_cs
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
  
  allocate(cs(n_wave,n_sp_photo), n_branch(n_sp_photo))
  allocate(ptitle(n_branch_max,n_sp_photo))
  allocate(enrgI(n_branch_max,n_sp_photo))
  allocate(charge_state(n_branch_max,n_sp_photo))
  allocate(im_photo_all(n_sp_photo))
  allocate(im_photo_product_all(n_product_max,n_branch_max,n_sp_photo))
  allocate(branch_ratio(n_wave,n_branch_max,n_sp_photo))

  ! initialize variables
  enrgI = zero
  charge_state = 0
  im_photo_all = 0
  im_photo_product_all = 0

  allocate(wave_lres(n_wave_lres), crs_lres(n_wave_lres))
  allocate(branch_ratio_lres(n_wave_lres))
  allocate(t_cs(n_wave_lres))
  
  ! read wavelength scale
  read(fid,'(A)') ts_line
  read(fid,'(10ES11.3)') wave_lres
  
  do i_sp = 1, n_sp_photo_lres
    read(fid,'(A12, I4)') ts_species1, n_branch(i_sp) ! header for species
    im_photo_all(i_sp) = find_name(ts_species1, sp_list)
    
    if (im_photo_all(i_sp) <= 0) then
      write(*,'("Error: Photo species not found: ",A12)') ts_species1
      stop
    end if
    
    read(fid,'(A)') ts_line ! "Total Absorption Cross Section (cm^2)"
    read(fid,*) (crs_lres(i_wave), i_wave=1,n_wave_lres)
    ! interpolate cross sections onto high resolution scale
    call intrp(wave_lres, crs_lres, wave, cs(:,i_sp))
    
    read(fid,'(A)') ts_line ! "Total Ion Cross Section (cm^2)"
    read(fid,*) (t_cs(i_wave), i_wave=1,n_wave_lres)
    
    do i_branch = 1, n_branch(i_sp)
      ! read reaction formula
      read(fid,'(A)') ts_line
      ptitle(i_branch,i_sp) = ts_line
      read(ts_line,'(5(A12, 3X), A12, F7.1)') &
        (ts_species2(i), i=1,6), enrgI(i_branch,i_sp)

      ! link products to species list
      i_p = 0
      do i = 3, 6 ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          i_p = i_p + 1
          im_photo_product_all(i_p,i_branch,i_sp) = tn_isp
        end if
        ! determine if ionization reaction (is electron produced?)
        if (tn_isp == n_sp) then
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
  deallocate(crs_lres, branch_ratio_lres, t_cs)
  
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
  read(fid_28N2,'(A)') ts_line ! "## Processed by..."
  read(fid_28N2,'(A)') ts_line ! "## CROSS_SECTION=..."
  read(fid_28N2,'(A)') ts_line ! "Wavenumber(cm-1) Cross-section(cm2)..."
  ! convert wavenumber scale (cm-1) to wavelength (angstrom) scale
  do i_wavenum = 1, n_wave_N2
    ! need to flip indices since wavelength is reciprocal of wavenumber
    i_wave = n_wave_N2 - i_wavenum + 1
    read(fid_28N2,*) wave_N2(i_wave), crs_N2(i_wave)
    ! convert cm-1 to angstrom
    wave_N2(i_wave) = cm_to_aa / wave_N2(i_wave)
  end do
  close(unit=fid_28N2)

  im_photo_all(i_sp) = find_name('N2          ', sp_list)
  n_branch(i_sp) = 1
  i_branch = 1
  ptitle(i_branch,i_sp) = 'N2           + hv           = '// &
    'N2D          + N            +              +             '
  enrgI(i_branch,i_sp) = zero
  charge_state(i_branch,i_sp) = zero
  branch_ratio(:,i_branch,i_sp) = one
  im_photo_product_all(1,i_branch,i_sp) = find_name('N2D         ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('N           ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0
  
  ! interpolate cross sections onto high resolution scale
  call intrp(wave_N2, crs_N2, wave, cs(:,i_sp))

  ! --- 29N2 ------------------------------------------------------------------
  i_sp = n_sp_photo_lres + 2
  open(newunit=fid_29N2, file='../data/photons/photoB-29N2.dat', &
    status='old', action='read')
  read(fid_29N2,'(A)') ts_line ! "## Processed by..."
  read(fid_29N2,'(A)') ts_line ! "## CROSS_SECTION=..."
  read(fid_29N2,'(A)') ts_line ! "Wavenumber(cm-1) Cross-section(cm2)..."
  ! convert wavenumber scale (cm-1) to wavelength (angstrom) scale
  do i_wavenum = 1, n_wave_N2
    ! need to flip indices since wavelength is reciprocal of wavenumber
    i_wave = n_wave_N2 - i_wavenum + 1
    read(fid_29N2,*) wave_N2(i_wave), crs_N2(i_wave)
    ! convert cm-1 to angstrom
    wave_N2(i_wave) = cm_to_aa / wave_N2(i_wave)
  end do
  close(unit=fid_29N2)

  im_photo_all(i_sp) = find_name('N2I         ', sp_list)
  n_branch(i_sp) = 2
  
  i_branch = 1
  ptitle(i_branch,i_sp) = 'N2I          + hv           = '// &
    'N2DI         + N            +              +             '
  enrgI(i_branch,i_sp) = zero
  charge_state(i_branch,i_sp) = zero
  branch_ratio(:,i_branch,i_sp) = half
  im_photo_product_all(1,i_branch,i_sp) = find_name('N2DI        ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('N           ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0

  i_branch = 2
  ptitle(i_branch,i_sp) = 'N2I          + hv           = '// &
    'N2D          + NI           +              +             '
  enrgI(i_branch,i_sp) = zero
  charge_state(i_branch,i_sp) = zero
  branch_ratio(:,i_branch,i_sp) = half
  im_photo_product_all(1,i_branch,i_sp) = find_name('N2D         ', sp_list)
  im_photo_product_all(2,i_branch,i_sp) = find_name('NI          ', sp_list)
  im_photo_product_all(3:4,i_branch,i_sp) = 0

  ! interpolate cross sections onto high resolution scale
  call intrp(wave_N2, crs_N2, wave, cs(:,i_sp))
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
  read(fid_COa,'(A)') ts_line ! "# 12C16O photoabsorption..."
  read(fid_COa,'(A)') ts_line ! "# wavelength(nm)..."
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
  read(fid_COb,'(A)') ts_line ! "# Photoabsorption..."
  read(fid_COb,'(A)') ts_line ! "# Based on the..."
  read(fid_COb,'(A)') ts_line ! "# Assembled by..."
  read(fid_COb,'(A)') ts_line ! "# wavelength(nm)..."
  do i_wave = n_wave_COa+1, n_wave_COa+n_wave_COb
    read(fid_COb,*) wave_CO(i_wave), &
      crs_CO(i_wave), crs_CO_diss(i_wave)
    ! convert nm to angstroms
    wave_CO(i_wave) = wave_CO(i_wave) * ten
  end do
  close(unit=fid_COb)
  
  im_photo_all(i_sp) = find_name('CO          ', sp_list)
  
  ! interpolate cross sections onto high resolution scale
  call intrp(wave_CO, crs_CO, wave, cs(:,i_sp))
  
  ! --- Photodissociation branch ----------------------------------------------
  i_branch = 1
  ptitle(i_branch,i_sp) = 'CO           + hv           = '// &
    'C            + O            +              +             '
  enrgI(i_branch,i_sp) = zero
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
  !  Local variables
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
  integer :: i_sp, i_branch, i, i_p
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
    
    if (im_photo_all(i_sp) <= 0) then
      write(*,'("Error: Photo species not found: ",A12)') ts_species1
      stop
    end if
    
    do i_branch = 1, n_branch(i_sp)
      ! read reaction formula
      read(fid,'(A)') ts_line
      ptitle(i_branch,i_sp) = ts_line
      read(ts_line,'(5(A12,3X), A12, ES11.3)') &
        (ts_species2(i), i=1,6), srateJ(i_branch,i_sp)
      
      ! link products to species list
      i_p = 0
      do i = 3, 6 ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          i_p = i_p + 1
          im_photo_product_all(i_p,i_branch,i_sp) = tn_isp
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
  !  Local variables
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
  
  read(fid,'(A)') ts_line ! e.g., "; High solar activity reference..."
  read(fid,'(A)') ts_line ! "; Values are..."
  read(fid,'(A)') ts_line ! "; wave (A)..."
  
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
  !  Local variables
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
  
  do i_wave = 1, n_sol_wave_lres
    ! determine wavelength of low-wavelength bin boundaries
    t_bin_low = sol_wave_lres(i_wave) - half * dwave_lres(i_wave)
    t_bin_high = sol_wave_lres(i_wave) + half * dwave_lres(i_wave)
    
    ! check if current low-resolution bin has any high-resolution fluxes
    if ((t_bin_low >= sol_wave_hres(1)) .and. &
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


subroutine read_rayleigh
! This subroutine sets up the calculation for Rayleigh scattering by CO2. 
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : intrp
  
  implicit none

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  ! file unit number for input file
  integer :: fid
  ! number of wavelengths for Rayleigh scattering
  integer :: n_wave
  ! wavelength scale for Rayleigh scattering
  real(wp), allocatable, dimension(:) :: wave
  ! cross sections for Rayleigh scattering
  real(wp), allocatable, dimension(:) :: cs
  
  ! loop variables
  integer :: i_wave, i_z

  !----------------------------------------------------------------------------
  !  Read Rayleigh scattering cross sections
  !----------------------------------------------------------------------------
  
  open(newunit=fid, file='../data/photons/CO2_Rayleigh.dat', &
    status='old', action='read')
    
  read(fid,*) n_wave ! number of entries
  allocate(wave(n_wave), cs(n_wave))
  do i_wave = 1, n_wave
    read(fid,*) wave(i_wave), cs(i_wave)
  end do
  
  close(unit=fid)
  
  !----------------------------------------------------------------------------
  !  Interpolate cross sections onto photo region C
  !----------------------------------------------------------------------------
  
  allocate(cs_ray(n_waveC), tau_ray(n_z,n_waveC))
  call intrp(wave, cs, waveC, cs_ray)

  do concurrent (i_wave = 1:n_waveC, i_z = 1:n_z)
    if (waveC(i_wave) > 500._wp) then
      tau_ray(i_z,i_wave) = cs_ray(i_wave) * col(i_z)
    else
      tau_ray(i_z,i_wave) = zero
    end if
  end do

  deallocate(wave, cs)

end subroutine read_rayleigh


subroutine paths1D
! This subroutine calculates the physical path lengths for radiative transfer
! through the atmosphere for various illumination conditions:
! 1: daylight
! 0: twilight (past terminator but illuminated at altitude)
! -1: night
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : locate
  
  implicit none
  
  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  ! radius for tangent point of light ray
  real(wp), allocatable, dimension(:) :: r_tan
  ! radius to opaque aerosol layer
  real(wp) :: rShadow
  
  ! loop variables
  integer :: i_z1, i_z2
  
  
  allocate(itan(n_z), ds(n_z,n_z))
  allocate(r_tan(n_z))
  
  rShadow = rPlanet + z_bot
  ds = 1.0e+33_wp ! high optical path by default

  if (cos_sza >= zero) then
  ! day
    illum = 1
    do concurrent (i_z1 = 1:n_z-1)
      r_tan(i_z1) = rz(i_z1) * sin_sza
      itan(i_z1) = i_z1
      do concurrent (i_z2 = i_z1:n_z-1)
        ds(i_z2,i_z1) = sqrt(rz(i_z2+1)**two - r_tan(i_z1)**two) &
          - sqrt(rz(i_z2)**two - r_tan(i_z1)**two)
      end do
    end do
  else if ((cos_sza < zero) .and. (rz(n_z)*sin_sza > rShadow)) then
  ! twilight (sin sza = cos solar elev ang)
    illum = 0
    ibot = locate(rShadow, rz*sin_sza) + 1
    do concurrent (i_z1 = ibot:n_z-1)
      r_tan(i_z1) = rz(i_z1) * sin_sza
      
      itan(i_z1) = locate(r_tan(i_z1), rz)
      if (rz(itan(i_z1)) < r_tan(i_z1)) then
        itan(i_z1) = itan(i_z1) + 1
      end if
      ds(itan(i_z1),i_z1) = sqrt(rz(itan(i_z1)+1)**two - r_tan(i_z1)**two)
      do concurrent (i_z2 = itan(i_z1)+1:n_z-1)
        ds(i_z2,i_z1) = sqrt(rz(i_z2+1)**two - r_tan(i_z1)**two) &
          - sqrt(rz(i_z2)**two - r_tan(i_z1)**two)
      end do
    end do
    r_tan(i_z1) = rz(n_z) * sin_sza
    itan(n_z) = locate(r_tan(i_z1),rz) + 1
  else if ((cos_sza < zero) .and. (rz(n_z)*sin_sza <= rShadow)) then
  ! night
    illum = -1
    ! just use default 1E33 value
  end if
  
end subroutine paths1D


subroutine cal_photo_rates
! This subroutine 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : find_bin

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  implicit none
  ! Column density: unit=cm-2; dim=(altitude level, species #)
  real(wp), allocatable, dimension(:,:) :: clm
  ! Scale height at top of atmosphere: unit=cm
  real(wp) :: Htop
  ! Chapman function parameters
  real(wp) :: ch_sin_a, ch_cos_a, ch_x
  ! Energies (ionization, photon, electron): unit=eV
  real(wp) :: E_ion, E_hv, E_e
!  real(wp), allocatable, dimension(:) :: Selsum
  
  ! loop variables
  integer :: i_sp1, i_sp2, i_z1, i_z2, i_branch, i_wave, i_r
  ! temporary / dummy variables
  integer :: ienrg
  real(wp) :: t

  !----------------------------------------------------------------------------
  ! Photolysis Rates at Top of Atmosphere
  !----------------------------------------------------------------------------
  if (cal_photoA) then
    do concurrent (i_wave = 1:n_waveA, i_sp1 = 1:n_sp_photoA)
      do i_branch = 1, n_branchA(i_sp1) 
        prtA(i_wave,i_branch,i_sp1) = branch_ratioA(i_wave,i_branch,i_sp1) &
          * csA(i_wave,i_sp1) * sol_fluxA(i_wave)
      end do
    end do
  end if

  if (cal_photoB) then
    do concurrent (i_wave = 1:n_waveB, i_sp1 = 1:n_sp_photoB)
      do i_branch = 1, n_branchB(i_sp1)
        prtB(i_wave,i_branch,i_sp1) = branch_ratioB(i_wave,i_branch,i_sp1) &
          * csB(i_wave,i_sp1) * sol_fluxB(i_wave)
      end do
    end do
  end if

  if (cal_photoC) then
    do concurrent (i_wave = 1:n_waveC, i_sp1 = 1:n_sp_photoC)
      do i_branch = 1, n_branchC(i_sp1)
        prtC(i_wave,i_branch,i_sp1) = branch_ratioC(i_wave,i_branch,i_sp1) &
          * csC(i_wave,i_sp1) * sol_fluxC(i_wave)
      end do
    end do
  end if

  !----------------------------------------------------------------------------
  ! Column Densities
  !----------------------------------------------------------------------------

  allocate(clm(n_z,n_sp))

  if (illum == 1) then    ! dayside
    do i_sp1 = 1, n_sp-1
      do i_z1 = 1, n_z
        ! element at top of atmosphere (includes contribution above top bin)

        ! Chapman function 
        ! Ch(a, x) = sqrt(pi*x/2) * exp(x/2*cos(a)^2) * erfc(sqrt(x/2)*cos(a))
        ! a is SZA at the point of interest
        ! x is radius measured from planetary center divided by scale height
        ch_sin_a = rz(i_z1) * sin_sza / rz(n_z)
        ch_cos_a = sqrt(one - ch_sin_a * ch_sin_a)
        ch_x = rz(n_z) / Ht(n_z,i_sp1)
        if (ch_x < 100._wp) then    ! curvature is important
          Htop = Ht(n_z,i_sp1) * sqrt(pi * half * ch_x) &
            * exp(half * ch_x * ch_cos_a * ch_cos_a) &
            * erfc(sqrt(half * ch_x) * ch_cos_a)
        else
          Htop = Ht(n_z,i_sp1) / cos_sza
        end if
        t = Htop * den(n_z,i_sp1)
        
        ! elements over rest of atmosphere
        do i_z2 = n_z-1, i_z1, -1
          t = t + half * (den(i_z2,i_sp1) + den(i_z2+1,i_sp1)) * ds(i_z2,i_z1)
        end do
        clm(i_z1,i_sp1) = t
      end do
   end do

  else if (illum == 0) then    ! twilight
    do i_sp1 = 1, n_sp-1
      clm(1:ibot-1,i_sp1) = 1.E33_wp
      do i_z1 = ibot, n_z
        ! element at top of atmosphere (includes contribution above top bin)
        ch_sin_a = rz(i_z1) * sin_sza / rz(n_z)
        ch_cos_a = sqrt(one - ch_sin_a * ch_sin_a)
        ch_x = rz(n_z) / Ht(n_z,i_sp1)
        if (ch_x < 100._wp) then    ! curvature is important
          Htop = Ht(n_z,i_sp1) * sqrt(pi * half * ch_x) &
            * exp(half * ch_x * ch_cos_a * ch_cos_a) &
            * erfc(sqrt(half * ch_x) * ch_cos_a)
        else
          Htop = Ht(n_z,i_sp1) / cos_sza
        end if
        t = Htop * den(n_z,i_sp1)
        
        do i_z2 = i_z1, n_z-1                  ! exterior to rz
          t = t + half * (den(i_z2,i_sp1) + den(i_z2+1,i_sp1)) * ds(i_z2,i_z1)
        end do
        do i_z2 = max(itan(i_z1),1), i_z1-1    ! interior to rz
          t = t + (den(i_z2,i_sp1)+den(i_z2+1,i_sp1)) * ds(i_z2,i_z1)
        end do
        clm(i_z1,i_sp1) = t
      end do
    end do

  else if (illum == -1) then    ! nightside
    clm(i_z1,i_sp1) = 1.E33_wp
  
  end if

  !----------------------------------------------------------------------------
  ! Optical Depths
  !----------------------------------------------------------------------------

  if (cal_photoA) then
    do i_z1 = 1, n_z
      do i_wave = 1, n_waveA
        t = zero
        do i_sp1 = 1, n_sp_photoA
          i_sp2 = im_photoA_all(i_sp1)
          t = t + csA(i_wave,i_sp1) * clm(i_z1,i_sp2)
        end do
        trnA(i_wave,i_z1) = exp(-t)
      end do
    end do
  end if

  if (cal_photoB) then
    do i_z1 = 1, n_z
      do i_wave = 1, n_waveB
        t = zero
        do i_sp1 = 1, n_sp_photoB
          i_sp2 = im_photoB_all(i_sp1)
          t = t + csB(i_wave,i_sp1) * clm(i_z1,i_sp2)
        end do
        trnB(i_wave,i_z1) = exp(-t)
      end do
    end do
  end if

  if (cal_photoC) then
    do i_z1 = 1, n_z
      do i_wave = 1, n_waveC
        t = zero
        do i_sp1 = 1, n_sp_photoC
          i_sp2 = im_photoC_all(i_sp1)
          t = t + csC(i_wave,i_sp1) * clm(i_z1,i_sp2)
        end do
        trnC(i_wave,i_z1) = exp(-(t + tau_ray(i_z1,i_wave) / cos_sza))
        ! with addition absorption from CO2 Rayleigh scattering
      end do
    end do
  end if

  !----------------------------------------------------------------------------
  ! Absorption Rates
  !----------------------------------------------------------------------------

  i_r = 0
  
  do i_sp1 = 1, n_sp_photoA
    do i_branch = 1, n_branchA(i_sp1)
      i_r = i_r + 1
      if (cal_photoA) then
        do i_z1 = 1, n_z
          t = zero
          do i_wave = 1, n_waveA
            t = t + prtA(i_wave,i_branch,i_sp1) * trnA(i_wave,i_z1)
          end do
          rph(i_r,i_z1) = diurnal_average * t
        end do
      end if
    end do
  end do

  do i_sp1 = 1, n_sp_photoB
    do i_branch = 1, n_branchB(i_sp1)
      i_r = i_r + 1
      if (cal_photoB) then  
        do i_z1 = 1, n_z
          t = zero
          do i_wave = 1, n_waveB
            t = t + prtB(i_wave,i_branch,i_sp1) * trnB(i_wave,i_z1)
          end do
          rph(i_r,i_z1) = diurnal_average * t
        end do
      end if
    end do
  end do
  
  do i_sp1 = 1, n_sp_photoC
    do i_branch = 1, n_branchC(i_sp1)
      i_r = i_r + 1
      if (cal_photoC) then
        do i_z1 = 1, n_z
          t = zero
          do i_wave = 1, n_waveC
            t = t + prtC(i_wave,i_branch,i_sp1) * trnC(i_wave,i_z1)
          end do
          rph(i_r,i_z1) = diurnal_average * t
        end do
      end if
    end do
  end do

  do i_sp1 = 1, n_sp_photoJ
    do i_branch = 1, n_branchJ(i_sp1)
      i_r = i_r + 1
      if ((illum >= 0) .and. cal_photoJ) then
        rph(i_r,1:n_z) = diurnal_average * srateJ(i_branch,i_sp1)
      else
        rph(i_r,1:n_z) = zero
      end if
    end do
  end do

  !----------------------------------------------------------------------------
  ! Photoelectron Production
  !----------------------------------------------------------------------------

  esrc = zero
  i_r = 0
  
  if (cal_photoA) then
    do i_sp1 = 1, n_sp_photoA
      i_sp2 = im_photoA_all(i_sp1)
      do i_branch = 1, n_branchA(i_sp1)
        i_r = i_r + 1
        if (charge_stateA(i_branch,i_sp1) > zero) then
          E_ion = hc / (enrgIA(i_branch,i_sp1) * aa_to_cm) * erg_to_ev
          do i_wave = 1, n_waveA
            E_hv = hc / (waveA(i_wave) * aa_to_cm) * erg_to_ev
            E_e = (E_hv - E_ion) / charge_stateA(i_branch,i_sp1)
            if (E_e > zero) then
              ienrg = find_bin(E_e, e_enrg, e_denrg) 
              do i_z1 = iz_bot, n_z
                esrc(i_r,i_z1,ienrg) = esrc(i_r,i_z1,ienrg) &
                  + (E_e / e_enrg(ienrg) / e_denrg(ienrg)) &
                  * prtA(i_wave,i_branch,i_sp1) &
                  * charge_stateA(i_branch,i_sp1) &
                  * den(i_z1,i_sp2) * trnA(i_wave,i_z1)
              end do
            end if
          end do
        end if
      end do
    end do
  end if

  if (cal_photoB) then
    do i_sp1 = 1, n_sp_photoB
      i_sp2 = im_photoB_all(i_sp1)
      do i_branch = 1, n_branchB(i_sp1)
        i_r = i_r + 1
        if (charge_stateB(i_branch,i_sp1) > zero) then
          E_ion = hc / (enrgIB(i_branch,i_sp1) * aa_to_cm) * erg_to_ev
          do i_wave = 1, n_waveB
            E_hv = hc / (waveB(i_wave) * aa_to_cm) * erg_to_ev
            E_e = (E_hv - E_ion) / charge_stateB(i_branch,i_sp1)
            if (E_e > zero) then
              ienrg = find_bin(E_e, e_enrg, e_denrg)
              do i_z1 = iz_bot, n_z
                esrc(i_r,i_z1,ienrg) = esrc(i_r,i_z1,ienrg) &
                  + (E_e / e_enrg(ienrg) / e_denrg(ienrg)) &
                  * prtB(i_wave,i_branch,i_sp1) &
                  * charge_stateB(i_branch,i_sp1) &
                  * den(i_z1,i_sp2) * trnB(i_wave,i_z1)
              end do
            end if
          end do
        end if
      end do
    end do
  end if

  if (cal_photoC) then
    do i_sp1 = 1, n_sp_photoC
      i_sp2 = im_photoC_all(i_sp1)
      do i_branch = 1, n_branchC(i_sp1)
        i_r = i_r + 1
        if (charge_stateC(i_branch,i_sp1) > zero) then
          E_ion = hc / (enrgIC(i_branch,i_sp1) * aa_to_cm) * erg_to_ev
          do i_wave = 1, n_waveC
            E_hv = hc / (waveC(i_wave) * aa_to_cm) * erg_to_ev
            E_e = (E_hv - E_ion) / charge_stateC(i_branch,i_sp1)
            if (E_e > zero) then
              ienrg = find_bin(E_e, e_enrg, e_denrg)
              do i_z1 = iz_bot, n_z
                esrc(i_r,i_z1,ienrg) = esrc(i_r,i_z1,ienrg) &
                  + (E_e / e_enrg(ienrg) / e_denrg(ienrg)) &
                  * prtC(i_wave,i_branch,i_sp1) &
                  * charge_stateC(i_branch,i_sp1) &
                  * den(i_z1,i_sp2) * trnC(i_wave,i_z1)
              end do
            end if
          end do
        end if
      end do
    end do
  end if

  esrc = diurnal_average * esrc

  Sel = zero
  do concurrent (ienrg = 1:n_e_enrg, i_z1 = 1:n_z)
    Sel(i_z1,ienrg) = sum(esrc(:,i_z1,ienrg))
  end do

!  ALLOCATE(Selsum(n_z))
!  Selsum = zero ! total electrons produced in altitude bin
!  do i_z1 = iz_bot, n_z
!    t = zero
!    do ienrg = 1, n_e_enrg
!      t = t + e_denrg(ienrg)*Sel(i_z1,ienrg)
!    end do
!    Selsum(i_z1) = t
!  end do

  deallocate(clm)

end subroutine cal_photo_rates


end module
