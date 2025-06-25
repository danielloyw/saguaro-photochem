! This subroutine reads in the settings for neutral and ion species 
! (nmolecules.settings, imolecules.settings), then allocate the relevant
! variables. 
subroutine read_species

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use types, only: wp => dp
  use constants
  use global_variables
  use subs, only : find_name

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  implicit none
  ! file unit numbers for neutral and ion settings
  integer :: fid_neu, fid_ion
  ! number of ions in chemical equilibrium
  integer :: ichk
  ! number of species in chemical equilibrium + electrons
  integer :: n_chem_e
  
  ! loop variables
  integer :: i_sp, nc, nd
  ! temporary / dummy variables
  integer :: tn1
  character(len=140) :: ts_header

  !----------------------------------------------------------------------------
  !  Read settings for neutral and ionic species
  !----------------------------------------------------------------------------

  open(newunit=fid_neu, file='nmolecules.settings', &
    status='old', action='read')
  open(newunit=fid_ion, file='imolecules.settings', &
    status='old', action='read')

  read(fid_neu,*) n_neu ! number of neutral species
  read(fid_ion,*) n_ion ! number of ion species

  n_sp = n_neu + n_ion + 1 ! total number of species (neutrals+ions+electrons)

  ! allocate variables in preparation for population with file read
  allocate(sp_list(0:n_sp), istat(n_sp), chrg(n_sp), mmw(n_sp), &
    nhyd(n_sp), ncar(n_sp), n14n(n_sp), n15n(n_sp), noxy(n_sp), &
    dtype(n_sp), ad(n_sp), sd(n_sp), phi(n_sp), sd_2(n_sp), sd_3(n_sp), &
    ibnd(n_sp,2), bval(n_sp,2))
  
  ! read in settings for neutrals
  read(fid_neu,"(A)") ts_header ! column headers
  do i_sp = 1, n_neu
    read(fid_neu,961) tn1, &
      sp_list(i_sp), istat(i_sp), chrg(i_sp), mmw(i_sp), &
      nhyd(i_sp), ncar(i_sp), n14n(i_sp), n15n(i_sp), noxy(i_sp), &
      dtype(i_sp), ad(i_sp), sd(i_sp), phi(i_sp), sd_2(i_sp), sd_3(i_sp), &
      ibnd(i_sp,1), bval(i_sp,1), ibnd(i_sp,2), bval(i_sp,2)
  end do

961 format(I4, &
      1X, A12, 1X, I1, 1X, I2, 1X, F7.1, &
      5I4, 2X, &
      I4, 1X, ES9.2, 1X, F6.3, 1X, ES9.2, 1X, F5.1, 1X, F6.1, 6X, &
      I1, 2X, ES10.3, 3X, I1, 2X, ES10.3)
  
  ! read in settings for ions
  read(fid_ion,"(A)") ts_header
  ichk = 0
  do i_sp = n_neu+1, n_neu+n_ion
     read(fid_ion,962) tn1, &
      sp_list(i_sp), istat(i_sp), chrg(i_sp), mmw(i_sp), &
      nhyd(i_sp), ncar(i_sp), n14n(i_sp), n15n(i_sp), noxy(i_sp), &
      ibnd(i_sp,1), bval(i_sp,1), ibnd(i_sp,2), bval(i_sp,2)
    if(istat(i_sp) == 1) ichk = ichk + 1
  end do
  
962 format(I4, &
      1X, A12, 1X, I1, 1X, I2, 1X, F7.1, &
      5I4, 4X, &
      I1, 2X, ES10.3, 3X, I1, 2X, ES10.3)

  close(unit=fid_neu)
  close(unit=fid_ion)

  sp_list(0)    = '            ' ! for total
  sp_list(n_sp) = 'E           ' ! electrons
  mmw(n_sp) = zero

  !----------------------------------------------------------------------------
  !  Create lists of species in chemical and diffusive equilibrium
  !----------------------------------------------------------------------------

  ! count number of species in chemical and diffusive equilibrium
  nc = 0; nd = 0
  do i_sp = 1, n_sp-1
    if(istat(i_sp) == 1) nc = nc + 1 ! species in chemical equilibrium
    if(istat(i_sp) == 2) nd = nd + 1 ! species diffusive equilibrium
  end do
  n_chem = nc
  n_diff = nd
  

  if(ichk > 0) then
    have_ions = .true.
    n_chem_e = n_chem + 1 ! with electrons
  else
    have_ions = .false.
    n_chem_e = n_chem
  end if

  ! index mapping between lists
  ! list of chemical species -> list of all species
  allocate(im_chem_all(n_chem_e))
  im_chem_all = 0
  ! list of all species -> list of chemical species
  allocate(im_all_chem(0:n_sp))
  im_all_chem = 0
  ! list of diffusing species -> list of all species
  allocate(im_diff_all(n_diff))
  im_diff_all = 0
  ! list of all species -> list of diffusing species
  allocate(im_all_diff(0:n_sp))
  im_all_diff = 0


  ! calculate mapping between list of diffusing species and list of all species
  nd = 0
  do i_sp = 1, n_sp-1
    if(istat(i_sp) == 2) then
      nd=nd+1
      im_diff_all(nd) = i_sp
      im_all_diff(i_sp) = nd
    end if
  end do

  write(*,"('DIFFUSING SPECIES')")
  if(size(im_diff_all) > 0) then
    write(*,"(10(1X,A12))") (sp_list(im_diff_all(nd)), nd=1, n_diff)
  else
    write(*,"('  NONE')")
  end if


  ! calculate mapping between list of chemical species and list of all species
  nc = 0
  do i_sp = 1, n_sp-1
    if(istat(i_sp) == 1) then ! species in chemical equilibrium
      nc=nc+1
      im_chem_all(nc) = i_sp
      im_all_chem(i_sp) = nc
    end if
  end do
  
  ! entry for electrons
  if(have_ions) then
    im_chem_all(n_chem_e) = n_sp
    im_all_chem(n_sp) = n_chem_e
  end if

  ! print list of chemical species to screen
  write(*,"('CHEMICAL SPECIES')")
  if(size(im_chem_all) > 0) then
    write(*,"(10(1X,A12))") (sp_list(im_chem_all(nc)), nc=1, n_chem)
  else
    write(*,"('  NONE')")
  end if

  !----------------------------------------------------------------------------
  !  Determine indices of key species
  !----------------------------------------------------------------------------
  
  iN2  = find_name('N2          ', sp_list)
  iCO2 = find_name('CO2         ', sp_list)
  iELE = find_name('E           ', sp_list)

end subroutine read_species
