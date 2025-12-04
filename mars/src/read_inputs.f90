module read_inputs

contains

subroutine read_species
! This subroutine reads in the settings for neutral and ion species 
! (nmolecules.settings, imolecules.settings), then allocate the relevant
! variables. 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use constants
  use global_variables
  use utils, only : find_name

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  implicit none
  ! file unit numbers for neutral and ion settings
  integer :: fid_neu, fid_ion
  ! number of ions
  integer :: ichk
  ! number of species in chemistry loop + electrons
  integer :: n_chem_e
  
  ! loop variables
  integer :: i_sp, n_c, n_d
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

  read(fid_neu,*) n_neu    ! number of neutral species
  read(fid_ion,*) n_ion    ! number of ion species

  ! total number of species (neutrals+ions+electrons)
  n_sp = n_neu + n_ion + 1

  ! allocate variables in preparation for population with file read
  allocate(sp_list(0:n_sp), istat(n_sp), chrg(n_sp), mmw(n_sp), &
    nhyd(n_sp), ncar(n_sp), n14N(n_sp), n15N(n_sp), noxy(n_sp), &
    dtype(n_sp), ad(n_sp), sd1(n_sp), phi(n_sp), sd2(n_sp), sd3(n_sp), &
    ibnd(n_sp,2), bval(n_sp,2))
  
  ! read in settings for neutrals
  read(fid_neu,'(A)') ts_header    ! column headers
  do i_sp = 1, n_neu
    read(fid_neu,961) tn1, &
      sp_list(i_sp), istat(i_sp), chrg(i_sp), mmw(i_sp), &
      nhyd(i_sp), ncar(i_sp), n14N(i_sp), n15N(i_sp), noxy(i_sp), &
      dtype(i_sp), ad(i_sp), sd1(i_sp), phi(i_sp), sd2(i_sp), sd3(i_sp), &
      ibnd(i_sp,1), bval(i_sp,1), ibnd(i_sp,2), bval(i_sp,2)
  end do

961 format(I4, &
      1X, A12, 1X, I1, 1X, I2, 1X, F7.1, &
      5I4, 2X, &
      I4, 1X, ES9.2, 1X, F6.3, 1X, ES9.2, 1X, F5.1, 1X, F6.1, 6X, &
      I1, 2X, ES10.3, 3X, I1, 2X, ES10.3)
  
  ! read in settings for ions
  read(fid_ion,'(A)') ts_header
  ichk = 0
  do i_sp = n_neu+1, n_neu+n_ion
     read(fid_ion,962) tn1, &
      sp_list(i_sp), istat(i_sp), chrg(i_sp), mmw(i_sp), &
      nhyd(i_sp), ncar(i_sp), n14N(i_sp), n15N(i_sp), noxy(i_sp), &
      ibnd(i_sp,1), bval(i_sp,1), ibnd(i_sp,2), bval(i_sp,2)
    if (istat(i_sp) == 1) ichk = ichk + 1
  end do
  
962 format(I4, &
      1X, A12, 1X, I1, 1X, I2, 1X, F7.1, &
      5I4, 4X, &
      I1, 2X, ES10.3, 3X, I1, 2X, ES10.3)

  close(unit=fid_neu)
  close(unit=fid_ion)

  sp_list(0)    = '            '    ! for total
  sp_list(n_sp) = 'E           '    ! electrons
  mmw(n_sp) = zero

  !----------------------------------------------------------------------------
  !  Create lists of species in chemistry and diffusion loops
  !----------------------------------------------------------------------------

  ! count number of species in chemistry and diffusion loops
  n_c = 0; n_d = 0
  do i_sp = 1, n_sp-1
    if (istat(i_sp) == 1) n_c = n_c + 1    ! species in chemistry loop
    if (istat(i_sp) == 2) n_d = n_d + 1    ! species in diffusion loop
  end do
  n_chem = n_c
  n_diff = n_d
  

  if (ichk > 0) then
    have_ions = .true.
    n_chem_e = n_chem + 1    ! with electrons
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
  n_d = 0
  do i_sp = 1, n_sp-1
    if (istat(i_sp) == 2) then
      n_d = n_d + 1
      im_diff_all(n_d) = i_sp
      im_all_diff(i_sp) = n_d
    end if
  end do

  write(*,'("DIFFUSING SPECIES")')
  if (size(im_diff_all) > 0) then
    write(*,'(10(1X,A12))') (sp_list(im_diff_all(n_d)), n_d=1, n_diff)
  else
    write(*,'("  NONE")')
  end if


  ! calculate mapping between list of chemical species and list of all species
  n_c = 0
  do i_sp = 1, n_sp-1
    if (istat(i_sp) == 1) then    ! species in chemistry loop
      n_c = n_c + 1
      im_chem_all(n_c) = i_sp
      im_all_chem(i_sp) = n_c
    end if
  end do
  
  ! entry for electrons
  if (have_ions) then
    im_chem_all(n_chem_e) = n_sp
    im_all_chem(n_sp) = n_chem_e
  end if

  ! print list of chemical species to screen
  write(*,'("CHEMISTRY SPECIES")')
  if (size(im_chem_all) > 0) then
    write(*,'(10(1X,A12))') (sp_list(im_chem_all(n_c)), n_c=1, n_chem)
  else
    write(*,'("  NONE")')
  end if

  !----------------------------------------------------------------------------
  !  Determine indices of key species
  !----------------------------------------------------------------------------
  
  iN2  = find_name('N2          ', sp_list)
  iCO2 = find_name('CO2         ', sp_list)
  iELE = find_name('E           ', sp_list)

end subroutine read_species


subroutine read_atmos
! This subroutine reads in the model atmosphere from atm1D.in, and calculates 
! other atmospheric state quantities. 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : find_name, locate

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  implicit none
  ! file unit numbers for atm1D.in
  integer :: fid
  
  ! species names in atm1D.in
  character(len=12), allocatable, dimension(:) :: sp_atm
  ! index mapping from sp_atm to sp_list
  integer, allocatable, dimension(:) :: im_atm_list
  ! does species from settings have density defined in atmosphere?
  logical, allocatable, dimension(:) :: has_den
  
  ! loop variables
  integer :: i_sp, i_z
  ! temporary / dummy variables
  character(len=128) :: ts_header
  character(len=12) :: ts_name
  integer :: tn_sp1, tn_sp2
  
  !----------------------------------------------------------------------------
  !  Read in model atmosphere from atm1D
  !----------------------------------------------------------------------------
  
  ! start with all species undefined, then flip when density found
  allocate(has_den(n_sp))
  has_den = .false.
  
  ! read in model atmosphere file atm1D.in
  open(newunit=fid, file='atm1D.in', status='old', action='read')
    
    read(fid,*) n_z, tn_sp1
    ! allocate variables
    allocate(sp_atm(tn_sp1), im_atm_list(tn_sp1))
    allocate(z(n_z), rz(0:n_z+1), grv(n_z), eK(n_z))
    allocate(Tn(n_z), Te(n_z), prs(n_z), mass(n_z), rho(n_z))
    allocate(den(n_z,0:n_sp), den_old(n_z,0:n_sp), vmr(n_z,n_sp))

    ! read species and associate them with sp_list from species settings files
    ! i.e., nmolecules.settings and imolecules.settings
    read(fid,'(A)') ts_header
    read(fid,'(10(2X,A12,1X))') (sp_atm(i_sp), i_sp=1, tn_sp1)
    do i_sp = 1, tn_sp1
      tn_sp2 = find_name(sp_atm(i_sp),sp_list)
      im_atm_list(i_sp) = tn_sp2
      if (tn_sp2 /= 0) has_den(tn_sp2) = .true.
    end do

    ! read in model atmosphere and populate relevant variables
    read(fid,'(A)') ts_header
    read(fid,*) (z(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (rz(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (grv(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (Tn(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (Te(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (prs(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (rho(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (mass(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (eK(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (den(i_z,0), i_z=1, n_z)    ! den(*,0) is for overall density
    do i_sp = 1, tn_sp1
      tn_sp2 = im_atm_list(i_sp)
      read(fid,'(A)') ts_name
      read(fid,*) (den(i_z,tn_sp2), i_z=1, n_z)
    end do
    
  close(unit=fid)
  
  ! altitude to opaque atmosphere layer (e.g., aerosols)
  z_bot = 0.E5_wp
  iz_bot = locate(z_bot, z)
  
  !----------------------------------------------------------------------------
  !  Recalculate quanitites to ensure consistency
  !----------------------------------------------------------------------------

  ! apply bottom boundary condition of stipulated mole ratio
  do i_sp = 1, n_diff
    tn_sp2 = im_diff_all(i_sp)
    if (ibnd(tn_sp2,1) == 3) then
      den(1,tn_sp2) = bval(tn_sp2,1) * den(1,0)
    end if
  end do

  ! set density to min value if species does not have density defined
  do concurrent (i_sp = 1:n_sp-1)
    if (.not. has_den(i_sp)) then
      den(:,i_sp) = eps
      write(*,'("Density not defined for ", A, &
        "... Assigning minimum density of ", ES8.1)') sp_list(i_sp), eps
    end if
  end do

  ! set density to min value if below minimum
  where (den < eps) den = eps
  
  ! reset electron density to sum of ion densities
  do concurrent (i_z = 1:n_z)
    den(i_z,n_sp) = sum(den(i_z, n_neu+1:n_sp-1))
  end do
  
  ! reset total density
  do concurrent (i_z = 1:n_z)
    den(i_z,0) = sum(den(i_z,1:n_sp-1))
  end do
  
  ! set density to fixed mole fraction according to settings
  do concurrent (i_sp = 1:n_sp-1)
    if (istat(i_sp) == 3) then
      den(:,i_sp) = bval(i_sp,1) * den(:,0)
    end if
  end do

  ! calculate mole fraction
  do concurrent (i_z = 1:n_z, i_sp = 1:n_sp)
    vmr(i_z,i_sp) = den(i_z,i_sp) / den(i_z,0)
  end do

  ! reset gravity
  do concurrent (i_z = 1:n_z)
    grv(i_z) = GM / rz(i_z)**two
  end do

  ! reset pressure to be consistent with N & T
  do concurrent (i_z = 1:n_z)
    prs(i_z) = kB * Tn(i_z) * den(i_z,0)
  end do

  ! calculate mass density, mean molecular mass
  do concurrent (i_z = 1:n_z)
    rho(i_z) = dot_product(mmw, den(i_z,1:n_sp)) * amu
    mass(i_z) = rho(i_z) / den(i_z,0) / amu
  end do

  !----------------------------------------------------------------------------
  !  Auxillary altitude arrays
  !----------------------------------------------------------------------------

  allocate(rz_mid(0:n_z), dr_mid(0:n_z), dr(n_z))
  allocate(Ht(n_z,0:n_sp), col(n_z))

  rz(0) = rz(1) - (rz(2) - rz(1))
  rz(n_z+1) = rz(n_z) + (rz(n_z) - rz(n_z-1))
  do concurrent (i_z = 0:n_z)
    rz_mid(i_z) = half * (rz(i_z+1) + rz(i_z))
    dr_mid(i_z) = rz(i_z+1) - rz(i_z)
  end do

  do concurrent (i_z = 1:n_z)
    dr(i_z) = rz_mid(i_z) - rz_mid(i_z-1)
  end do

  !----------------------------------------------------------------------------
  !  Scale heights
  !----------------------------------------------------------------------------
  
  ! for each species
  do concurrent (i_z = 1:n_z, i_sp = 1:n_sp-1)
    Ht(i_z,i_sp) = kB * tn(i_z) / grv(i_z) / mmw(i_sp) / amu
  end do
  
  ! mean scale height
  do concurrent (i_z = 1:n_z)
    Ht(i_z,0) = kB * tn(i_z) / grv(i_z) / mass(i_z) / amu
  end do

  !----------------------------------------------------------------------------
  !  Column density above
  !----------------------------------------------------------------------------

  col(n_z) = den(n_z,0) * Ht(n_z,0)
  do i_z = n_z-1, 1, -1
    col(i_z) = col(i_z+1) + half * (den(i_z,0) + den(i_z+1,0)) * dr_mid(i_z)
  end do

  deallocate(sp_atm, im_atm_list, has_den)
end subroutine read_atmos


subroutine read_reactions
! This subroutine reads in the list of reactions and their rate constants in 
! nreactions.csv and ireactions.csv. 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : find_name

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  implicit none
  ! file unit numbers for neutral and ion reactions
  integer :: fid_nrct, fid_irct
  ! number of neutral and ion reactions
  integer :: n_nrct, n_irct
  ! constants for calculating reaction rate coefficients:
  ! dim=(constant 1-10, reaction #)
  real(wp), allocatable, dimension(:,:) :: rkn
  ! rate constants (k0, kinf, krad)
  real(wp) :: rk0, rkInf, rkRad
  ! reduced pressure (k0 * [M] / kinf)
  real(wp) :: P_red
  ! constants for Baulch et al. [2005] broadening factor
  real(wp) :: FC, NF, CF, FF
  
  ! loop variables
  integer :: i_r, i, i_r2, i_z
  ! temporary / dummy variables
  character(len=256) :: ts_line, ts1 
  integer :: tn1, tm_chem_type
  character(len=12), dimension(5) :: ts_species
  real(wp), dimension(10) :: t_rk
  
  ! Tabulated rates variables (deactivated)
  ! integer :: ntab_rct
  ! integer :: nt
  ! integer, allocatable, dimension(:,:) :: itab_rct
  ! logical :: ltab

  !----------------------------------------------------------------------------
  !  Read Tabulated Reaction Rate Coefficients (deactivated for now)
  !----------------------------------------------------------------------------

  ! OPEN(unit=40,file='../data/reactions/nreactions.tab',status='old')
     ! READ(40,*) ntab_rct
     ! ALLOCATE(itab_rct(5,ntab_rct),nprs_rct(ntab_rct),ntmp_rct(ntab_rct), &
          ! plog_rct(nprs_tab_max,ntab_rct),tmp_rct(ntmp_tab_max,ntab_rct), &
          ! rct_tab(nprs_tab_max,ntmp_tab_max,ntab_rct))
     ! DO i_r2 = 1, ntab_rct
        ! READ(40,"(A12,3X,A12,3X,A12,3X,A12,3X,A12)") (fn(i),i=1,5)
        ! DO i = 1, 5
           ! itab_rct(i,i_r2) = FIND_NAME(fn(i),name)
        ! END DO
        ! READ(40,*) nprs_rct(i_r2),ntmp_rct(i_r2)
        ! READ(40,*) (tmp_rct(nt,i_r2),nt=1,ntmp_rct(i_r2))
        ! DO np = 1, nprs_rct(i_r2)
           ! READ(40,*) rdum,(rct_tab(np,nt,i_r2),nt=1,ntmp_rct(i_r2))
           ! plog_rct(np,i_r2) = LOG(rdum)
        ! END DO
     ! END DO
  ! CLOSE(unit=40)
  ! rct_tab = LOG(rct_tab)
  
  !----------------------------------------------------------------------------
  !  Read reaction rate constants
  !----------------------------------------------------------------------------
   
  open(newunit=fid_nrct, file='../data/reactions/nreactions.csv', & 
    status='old', action='read')
  open(newunit=fid_irct, file='../data/reactions/ireactions.csv', &
    status='old', action='read')

  ! determine number of entries in neutral reactions list
  read(fid_nrct,'(A)') ts_line
  read(ts_line,*) n_nrct, ts1
  read(fid_nrct,'(A)') ts_line
  read(fid_nrct,'(A)') ts_line

  ! count how many of the entries are actually active (type /= 0)
  n_rct = 0
  do i_r = 1, n_nrct
    read(fid_nrct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) n_rct = n_rct + 1
  end do

  rewind(unit=fid_nrct)

  ! determine number of entries in ion reactions list
  read(fid_irct,'(A)') ts_line
  read(ts_line,*) n_irct, ts1
  read(fid_irct,'(A)') ts_line
  read(fid_irct,'(A)') ts_line

  ! count how many of the entries are actually active (type /= 0)
  do i_r = 1, n_irct
    read(fid_irct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) n_rct = n_rct + 1
  end do

  rewind(unit=fid_irct)

  ! now allocate global variables and actually populate them
  allocate(ctitle(n_rct), chem_type(n_rct), rkn(10,n_rct), ireactant(5,n_rct))
  ! allocate(ntab(n_rct))

  ! --- Neutral reactions -----------------------------------------------------
  read(fid_nrct,'(A)') ts_line
  read(fid_nrct,'(A)') ts_line
  read(fid_nrct,'(A)') ts_line  
  i_r2 = 0
  do i_r = 1, n_nrct
    read(fid_nrct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) then
      i_r2 = i_r2 + 1
      chem_type(i_r2) = tm_chem_type
      ctitle(i_r2) = ts_species(1)//' + '//ts_species(2)//' = '// & 
        ts_species(3)//' + '//ts_species(4)//' + '//ts_species(5)
      ! determine indices for products
      do concurrent (i = 3:5)
        ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
      end do
      
      ! determine indices for reactants
      if (abs(chem_type(i_r2)) == 1) then
      ! unimolecular reaction
        rkn(:,i_r2) = t_rk
        i = 1
        ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
        if (ireactant(i,i_r2) <= 0) then
          write(*,'(I6, ":", 2X, A)') i_r2, ctitle(i_r2)
          write(*,'("Error: Reactant ", I2, " not found: ", A12)') & 
            i, ts_species(i)
          stop
        end if
        ireactant(2,i_r2) = 0

      else if ((abs(chem_type(i_r2)) > 1) .and. (abs(chem_type(i_r2)) < 6)) then
      ! bimolecular & trimolecular reaction
        rkn(:,i_r2) = t_rk
        do i = 1, 2
          ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
          if (ireactant(i,i_r2) <= 0) then
            write(*,'(I6, ":", 2X, A)') i_r2, ctitle(i_r2)
            write(*,'("Error: Reactant ", I2, " not found: ", A12)') & 
              i, ts_species(i)
            stop
          end if
        end do
         
      ! Tabulated reaction rates (deactivated for now)
      ! else if (abs(chem_type(i_r2)) == 6) then
        ! do i = 1, 2
          ! ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
          ! if (ireactant(i,i_r2) <= 0) then
            ! write(*,"(I6,':',2X,A)") i_r2, ctitle(i_r2)
            ! write(*,"('Error: Reactant ',I2,' not found: ',A12)") & 
              ! i, ts_species(i)
            ! stop
          ! end if
        ! end do
        ! do nt = 1, ntab_rct
          ! ltab = .true.
          ! do i = 1, 5
            ! if (ireactant(i,i_r2) /= itab_rct(i,nt)) ltab = .false.
        ! end do
        ! if (ltab) then
          ! ntab(i_r2) = nt
          ! EXIT
        ! end if
      ! end do
      ! if (.not.ltab) then
        ! write(*,"(I6,':',2X,A)") i_r2,ctitle(i_r2)
        ! write(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") i_r2
        ! STOP
      ! end if

      else 
        write(*,'("Error: Invalid setting for neutral reaction ", I0, &
          "! Exiting...")') i_r
        stop
      end if
    end if
  end do
  close(unit=fid_nrct)

  ! --- Ion reactions ---------------------------------------------------------
  read(fid_irct,'(A)') ts_line
  read(fid_irct,'(A)') ts_line
  read(fid_irct,'(A)') ts_line  
  do i_r = 1, n_irct
    read(fid_irct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) then
      i_r2 = i_r2 + 1
      chem_type(i_r2) = tm_chem_type
      ctitle(i_r2) = ts_species(1)//' + '//ts_species(2)//' = '// & 
        ts_species(3)//' + '//ts_species(4)//' + '//ts_species(5)
      ! determine indices for products
      do concurrent (i = 3:5)
        ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
      end do
      
      ! determine indices for reactants
      if (abs(chem_type(i_r2)) == 1) then
      ! unimolecular reaction
        rkn(:,i_r2) = t_rk
        i = 1
        ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
        if (ireactant(i,i_r2) <= 0) then
          write(*,'(I6, ":", 2X, A)') i_r2, ctitle(i_r2)
          write(*,'("Error: Reactant ", I0, " not found: ", A12)') &
            i, ts_species(i)
          stop
        end if
        ireactant(2,i_r2) = 0

      else if ((abs(chem_type(i_r2)) > 1) .and. (abs(chem_type(i_r2)) < 6)) then
      ! bimolecular & trimolecular reaction
      rkn(:,i_r2) = t_rk
        do i = 1, 2
          ireactant(i,i_r2) = find_name(ts_species(i), sp_list)
          if (ireactant(i,i_r2) <= 0) then
            write(*,'(I6, ":", 2X, A)') i_r2, ctitle(i_r2)
            write(*,'("Error: Reactant ", I0, " not found: ", A12)') & 
              i, ts_species(i)
            stop
          end if
        end do

      ! Tabulated reaction rates (deactivated for now)
      ! else if (ABS(chem_type(i_r2)) == 6) then
        ! do i = 1, 2
          ! ireactant(i,i_r2) = FIND_NAME(ts_species(i),sp_list)
            ! if (ireactant(i,i_r2) <= 0) then
              ! write(*,"(I6,':',2X,A)") i_r2,ctitle(i_r2)
              ! write(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)") i_r2,i,ts_species(i)
              ! STOP
            ! end if
          ! end do
           
          ! do nt = 1, ntab_rct
            ! ltab = .true.
            ! do i = 1, 5
              ! if (ireactant(i,i_r2) /= itab_rct(i,nt)) ltab = .false.
            ! end do
            ! if (ltab) then
              ! ntab(i_r2) = nt
              ! EXIT
            ! end if
          ! end do
          ! if (.not.ltab) then
            ! write(*,"(I6,':',2X,A)") i_r2,ctitle(i_r2)
            ! write(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") i_r2
            ! STOP
          ! end if

      else 
        write(*,'("Error: Invalid setting for ion reaction ", I0, &
          "! Exiting...")') i_r
        stop
      end if
    end if
  end do
  close(unit=fid_irct)

  ! deallocate(itab_rct)

  !----------------------------------------------------------------------------
  !  Calculate reaction coefficients
  !----------------------------------------------------------------------------
  
  allocate(rk(n_rct,n_z), rate_rct(n_rct,n_z)) 
  
  do i_r = 1, n_rct
    
    ! --- Neutral reactions ---------------------------------------------------
    if (chem_type(i_r) == 1) then
    ! unimolecular reactions
    ! rate coefficient is simply rate constant
      do concurrent (i_z = 1:n_z)
        rk(i_r,i_z) = rkn(1,i_r)
      end do

    else if (chem_type(i_r) == 2) then
    ! bimolecular reactions
    ! k1 * T^k2 * exp(k3 * T)
      do concurrent (i_z = 1:n_z)
        rk(i_r,i_z) = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * exp(rkn(3,i_r) / Tn(i_z))
      end do

    else if (chem_type(i_r) == 3) then
    ! association reactions
    ! k0 * [M] * kinf / (k0 * [M] + kinf)
    ! k0 and kinf are in the form k1 * T^k2 * exp(k3 * T)
    ! F is the Troe broadening factor
    ! log F = log Fc / (1 + ((log(P)+C) / (N - 0.14*(log(P)+C)))^2)
    ! where
    ! reduced pressure P = k0 * [M] / kinf
    ! N = 0.75 - 1.27 * log Fc
    ! C = -0.4 - 0.67 * log Fc
      do i_z = 1, n_z
        rkInf = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * exp(rkn(3,i_r) / Tn(i_z))
        rk0 = rkn(4,i_r) * (Tn(i_z) ** rkn(5,i_r)) &
          * exp(rkn(6,i_r) / Tn(i_z))
        P_red = rk0 * den(i_z,0) / rkInf
        if (rkn(10,i_r) == zero) then
          rk(i_r,i_z) = rk0 * den(i_z,0) / (one + P_red)
        else
          FC = rkn(10,i_r)
          NF = 0.75_wp - 1.27_wp * log10(FC)
          CF = -0.4_wp - 0.67_wp * log10(FC)
          FF = 10._wp ** (log10(FC) / (one + &
            ((log10(P_red) + CF) / (NF - 0.14_wp * (log10(P_red) + CF))) ** 2))
          rk(i_r,i_z) = FF * rk0 * den(i_z,0) / (one + P_red)
          end if
      end do
!? check if reaction types should be 3 rather than 4
    else if (chem_type(i_r) == 4) then
    ! association & radiative association reactions
      do i_z = 1, n_z
        rkInf = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * exp(rkn(3,i_r) / Tn(i_z))
        rk0 = rkn(4,i_r) * (Tn(i_z) ** rkn(5,i_r)) &
          * exp(rkn(6,i_r) / Tn(i_z))
        rkRad = rkn(7,i_r) * (Tn(i_z) ** rkn(8,i_r)) &
          * exp(rkn(9,i_r) / Tn(i_z))
        P_red = rk0 * den(i_z,0) / rkInf
        if (rkn(10,i_r) == zero) then
          rk(i_r,i_z) = rkRad + rk0 * den(i_z,0) / (one + P_red)
        else
          FC = rkn(10,i_r)
          NF = 0.75_wp - 1.27_wp * log10(FC)
          CF = -0.4_wp - 0.67_wp * log10(FC)
          FF = 10._wp ** (log10(FC) / (one + &
            ((log10(P_red) + CF) / (NF - 0.14_wp * (log10(P_red) + CF))) ** 2))
          rk(i_r,i_z) = min(rkInf, &
            rkRad + FF * rk0 * den(i_z,0) / (one + P_red))
        end if
      end do

    else if (chem_type(i_r) == 5) then
    ! association reactions (Sander's formula)
    ! k0 * [M] * kinf / (k0 * [M] + kinf)
    !   * 0.6 ^ (1 / (1 + (log10(k0 * [M] / kinf))^2)
      do i_z = 1, n_z
        rkInf = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * exp(rkn(3,i_r) / Tn(i_z))
        rk0 = rkn(4,i_r) * (Tn(i_z) ** rkn(5,i_r)) &
          * exp(rkn(6,i_r) / Tn(i_z))
        P_red = rk0 * den(i_z,0) / rkInf
        rk(i_r,i_z) = rk0 * den(i_z,0) / (one + P_red) &
          * 0.6_wp ** (one / (one + log10(P_red) ** two))
      end do

    ! Tabulated reaction rates (deactivated for now)
    ! else if (chem_type(i_r) == 6) then
      ! nt = ntab(i_r)
      ! do i_z = 1, n_z
        ! if (Tn(i_z) < tmp_rct(1,nt)) then
          ! nt1 = 1
          ! nt2 = 2
          ! rtmp = zero
        ! else if (Tn(i_z) > tmp_rct(ntmp_rct(nt),nt)) then
          ! nt1 = ntmp_rct(nt)-1
          ! nt2 = nt1 + 1
          ! rtmp = one
        ! else
          ! nt1 = LOCATE(Tn(i_z),tmp_rct(1:ntmp_rct(nt),nt),)
          ! nt2 = nt1 + 1
          ! rtmp = (Tn(i_z)-tmp_rct(nt1,nt))/(tmp_rct(nt2,nt)-tmp_rct(nt1,nt))
        ! end if
        ! if (log(prs(i_z)) < plog_rct(1,nt)) then
          ! rprs = zero
          ! np1 = 1
          ! np2 = 2
        ! else if (log(prs(i_z)) > plog_rct(nprs_rct(nt),nt)) then
          ! rprs = one
          ! np1 = nprs_rct(nt)-1
          ! np2 = np1 + 1
        ! else
          ! np1 = LOCATE(log(prs(i_z)),plog_rct(1:nprs_rct(nt),nt))
          ! np2 = np1 + 1
          ! rprs = (log(prs(i_z))-plog_rct(np1,nt))/(plog_rct(np2,nt)-plog_rct(np1,nt))
        ! end if
        ! rk(i_r,i_z) = (one-rtmp)*(one-rprs)*rct_tab(np1,nt1,nt) &
                 ! + rtmp*(one-rprs)*rct_tab(np2,nt1,nt) &
                 ! + (one-rtmp)*rprs*rct_tab(np1,nt2,nt) &
                 ! + rtmp*rprs*rct_tab(np2,nt2,nt) 
        ! rk(i_r,i_z) = exp(rk(i_r,i_z))
      ! end do

    ! --- Ion reactions -------------------------------------------------------
    else if (chem_type(i_r) == -1) then
!? why is there a factor of two?
    ! unimolecular reaction
      do concurrent (i_z = 1:n_z)
        rk(i_r,i_z) = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * two * exp(rkn(3,i_r) / Tn(i_z))
      end do

    else if (chem_type(i_r) == -2) then
    ! normal two-body reaction
      do concurrent (i_z = 1:n_z)
        rk(i_r,i_z) = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * two * exp(rkn(3,i_r) / Tn(i_z))
      end do

    else if (chem_type(i_r) == -3) then
    ! 3-body reaction
      do concurrent (i_z = 1:n_z)
        rk(i_r,i_z) = rkn(1,i_r) * (Tn(i_z) ** rkn(2,i_r)) &
          * two * exp(rkn(3,i_r) / Tn(i_z))
      end do

    else if (chem_type(i_r) == -4) then
    ! electron recombination
      do concurrent (i_z = 1:n_z)
        rk(i_r,i_z) = rkn(1,i_r) * (Te(i_z) ** rkn(2,i_r)) &
          * two * exp(rkn(3,i_r) / Te(i_z))
      end do
    
    else if (chem_type(i_r) .ne. 0) then
      write(*,'("Error! Invalid setting for reaction ", I0, ": ", A)') &
        i_r, ctitle(i_r)
      write(*,'("Exiting...")')
      stop
    end if

  end do

end subroutine read_reactions

end module
