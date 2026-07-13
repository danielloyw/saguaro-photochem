module output

contains

subroutine write_output
! This subroutine writes the output files based on the settings in 
! output.settings.

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only: heapsort
  use photo, only: cal_photo_rates
  use electrons, only: cal_electrons

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  implicit none
  ! settings for what output files to produce
  character(len=3), dimension(15) :: out_settings
  ! output directory path
  character(len=20) :: output_dir
  
  ! --- File unit numbers -----------------------------------------------------
  ! output settings
  integer :: fid_settings
  ! atm1D.out, atm1D.csv
  integer :: fid_atm1, fid_atm2
  ! diff.csv
  integer :: fid_diff
  ! solar_flux.out
  integer :: fid_solar
  ! eflux.out
  integer :: fid_eflux
  ! ratecoeff.out
  integer :: fid_ratecoeff
  ! chemrates.out
  integer :: fid_chemrates
  ! colrates.out
  integer :: fid_colrates
  ! photorates.out
  integer :: fid_photorates
  ! COphotofreqB.out
  integer :: fid_COphotofreqB
  ! pcolrates.out
  integer :: fid_pcolrates
  ! elerates.out
  integer :: fid_elerates
  ! ecolrates.out
  integer :: fid_ecolrates
  ! balance.out
  integer :: fid_balance
  ! species-specific files
  integer :: fid_species

  ! column headings
  character(len=15), allocatable, dimension(:) :: heads
  
  ! total number of wavelength bins for spectral regions A, B (low-res), C
  integer :: n_wave
  ! combined wavelength scale with all 3 spectral regions
  real(wp), allocatable, dimension(:) :: wave
  ! combined solar flux for spectral ranges A, B, C: 
  ! unit=photons cm-2 s-1 bin-1; dim=(wavelength bin, altitude level)
  real(wp), allocatable, dimension(:,:) :: solar_flux
  ! CO photodissociation rates in region B: unit=cm-3 s-1 (bin width in nm)-1
  ! dim=(wavelength bin, altitude level)
  real(wp), allocatable, dimension(:,:) :: COphotorate
  
  ! --- Column quantities -----------------------------------------------------
  ! reaction rates: unit=cm-2 s-1;
  ! dim=(reaction #, above aerosols / below aerosols / total)
  real(wp), allocatable, dimension(:,:) :: rate_col
  ! density: unit=cm-2; dim=(species #)
  real(wp), allocatable, dimension(:) :: den_col
  ! production and loss from chemical processes: unit=cm-2 s-1; dim=(species #)
  real(wp), allocatable, dimension(:) :: pr_c_col, ls_c_col
  ! production and loss from photo processes: unit=cm-2 s-1; dim=(species #)
  real(wp), allocatable, dimension(:) :: pr_p_col, ls_p_col
  ! production and loss from electron processes: unit=cm-2 s-1; dim=(species #)
  real(wp), allocatable, dimension(:) :: pr_e_col, ls_e_col
  ! net production and loss: unit=cm-2 s-1; dim=(species #)
  real(wp), allocatable, dimension(:) :: pr_col, ls_col
  ! diffusive fluxes at top and bottom: unit=cm-2 s-1; dim=(species #)
  real(wp), allocatable, dimension(:) :: flux_top_col, flux_bot_col
  ! balance: unit=cm-2 s-1; dim=(species #)
  real(wp), allocatable, dimension(:) :: bal_col
  ! time constant to achieve balance = 0: unit=s; dim(species #)
  real(wp), allocatable, dimension(:) :: bal_time

  ! --- Oxygen (atoms) quantities ---------------------------------------------
  ! density: unit=cm-3; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: den_oxy
  ! production and loss from chemical processes: unit=cm-3 s-1;
  ! dim=(altitude level)
  real(wp), allocatable, dimension(:) :: pr_c_oxy, ls_c_oxy
  ! production and loss from photo processes: unit=cm-3 s-1;
  ! dim=(altitude level)
  real(wp), allocatable, dimension(:) :: pr_p_oxy, ls_p_oxy
  ! production and loss from electron processes: unit=cm-3 s-1;
  ! dim=(altitude level)
  real(wp), allocatable, dimension(:) :: pr_e_oxy, ls_e_oxy
  ! net production and loss: unit=cm-3 s-1; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: pr_oxy, ls_oxy
  ! column production and loss from chemical processes: unit=cm-2 s-1
  real(wp) :: pr_c_col_oxy, ls_c_col_oxy
  ! column production and loss from photo processes: unit=cm-2 s-1
  real(wp) :: pr_p_col_oxy, ls_p_col_oxy
  ! column production and loss from electron processes: unit=cm-2 s-1
  real(wp) :: pr_e_col_oxy, ls_e_col_oxy
  ! diffusive fluxes at top and bottom: unit=cm-2 s-1
  real(wp) :: flux_top_col_oxy, flux_bot_col_oxy
  ! column net production and loss: unit=cm-2 s-1
  real(wp) :: pr_col_oxy, ls_col_oxy
  ! flux: unit=cm-2 s-1; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: flux_oxy
  ! flux divergence: unit=cm-3 s-1; dim=(altitude level)
  real(wp), allocatable, dimension(:) :: div_flux_oxy
  ! column balance: unit=cm-2 s-1
  real(wp) :: bal_col_oxy
  
  ! --- Hydrogen (atoms) quantities -----------------------------------------------
  ! density: unit=cm-3; dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: den_hyd
  ! production and loss from chemical processes: unit=cm-3 s-1;
  ! dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: pr_c_hyd, ls_c_hyd
  ! production and loss from photo processes: unit=cm-3 s-1;
  ! dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: pr_p_hyd, ls_p_hyd
  ! production and loss from electron processes: unit=cm-3 s-1;
  ! dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: pr_e_hyd, ls_e_hyd
  ! net production and loss: unit=cm-3 s-1; dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: pr_hyd, ls_hyd
  ! flux: unit=cm-2 s-1; dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: flux_hyd
  ! flux divergence: unit=cm-3 s-1; dim=(altitude level, odd/even)
  real(wp), allocatable, dimension(:,:) :: div_flux_hyd
  ! column production and loss from chemical processes: unit=cm-2 s-1;
  ! dim=(odd/even)
  real(wp), dimension(2) :: pr_c_col_hyd, ls_c_col_hyd
  ! column production and loss from photo processes: unit=cm-2 s-1;
  ! dim=(odd/even)
  real(wp), dimension(2) :: pr_p_col_hyd, ls_p_col_hyd
  ! column production and loss from electron processes: unit=cm-2 s-1;
  ! dim=(odd/even)
  real(wp), dimension(2) :: pr_e_col_hyd, ls_e_col_hyd
  ! column net production and loss: unit=cm-2 s-1; dim=(odd/even)
  real(wp), dimension(2) :: pr_col_hyd, ls_col_hyd
  ! diffusive fluxes at top and bottom: unit=cm-2 s-1; dim=(odd/even)
  real(wp), dimension(2) :: flux_top_col_hyd, flux_bot_col_hyd
  ! column balance: unit=cm-2 s-1; dim=(odd/even)
  real(wp), dimension(2) :: bal_col_hyd
  
  ! loop variables
  integer :: i, i_sp, i_z, i_wave, i_waveB, i_enrg, i_r, i_h
  ! temporary / dummy variables
  character(len=20) :: ts_header
  integer, allocatable, dimension(:) :: t_sort_order
  integer :: ti_sp1, ti_sp2, ti_sp3, ti_sp4, ti_sp5
  integer :: ti_sp_c1, ti_sp_c2, ti_sp_c3, ti_sp_c4, ti_sp_c5
  
  !----------------------------------------------------------------------------
  !  Recalculate rates
  !----------------------------------------------------------------------------

  cal_photoA = .true.
  cal_photoB = .true.
  cal_photoC = .true.
  cal_photoJ = .true.
  call cal_photo_rates
  call cal_electrons
  
  pr_p = zero; ls_p = zero; rate_prct = zero
  pr_e = zero; ls_e = zero; rate_erct = zero
  pr_c = zero; ls_c = zero; rate_rct = zero
  where (rk < 1.0E-99_wp) rk = zero
  do i_z = 1, n_z
    ! Chemical reactions
    do i_r = 1, n_rct
      if (chem_type(i_r) == 1) then    ! unimolecular
        ti_sp1 = ireactant(1,i_r); ti_sp_c1 = im_all_chem(ti_sp1)
        ti_sp3 = ireactant(3,i_r); ti_sp_c3 = im_all_chem(ti_sp3)
        ti_sp4 = ireactant(4,i_r); ti_sp_c4 = im_all_chem(ti_sp4)
        ti_sp5 = ireactant(5,i_r); ti_sp_c5 = im_all_chem(ti_sp5)
        rate_rct(i_r,i_z) = rk(i_r,i_z) * den(i_z,ti_sp1)
        ls_c(i_z,ti_sp1) = ls_c(i_z,ti_sp1) + rate_rct(i_r,i_z)
        pr_c(i_z,ti_sp3) = pr_c(i_z,ti_sp3) + rate_rct(i_r,i_z)
        pr_c(i_z,ti_sp4) = pr_c(i_z,ti_sp4) + rate_rct(i_r,i_z)
        pr_c(i_z,ti_sp5) = pr_c(i_z,ti_sp5) + rate_rct(i_r,i_z)
      else
        ti_sp1 = ireactant(1,i_r); ti_sp_c1 = im_all_chem(ti_sp1)
        ti_sp2 = ireactant(2,i_r); ti_sp_c2 = im_all_chem(ti_sp2)
        ti_sp3 = ireactant(3,i_r); ti_sp_c3 = im_all_chem(ti_sp3)
        ti_sp4 = ireactant(4,i_r); ti_sp_c4 = im_all_chem(ti_sp4)
        ti_sp5 = ireactant(5,i_r); ti_sp_c5 = im_all_chem(ti_sp5)
        rate_rct(i_r,i_z) = rk(i_r,i_z) * den(i_z,ti_sp1) * den(i_z,ti_sp2)
        ls_c(i_z,ti_sp1) = ls_c(i_z,ti_sp1) + rate_rct(i_r,i_z)
        ls_c(i_z,ti_sp2) = ls_c(i_z,ti_sp2) + rate_rct(i_r,i_z)
        pr_c(i_z,ti_sp3) = pr_c(i_z,ti_sp3) + rate_rct(i_r,i_z)
        pr_c(i_z,ti_sp4) = pr_c(i_z,ti_sp4) + rate_rct(i_r,i_z)
        pr_c(i_z,ti_sp5) = pr_c(i_z,ti_sp5) + rate_rct(i_r,i_z)
      end if
    end do
    
    ! Photo reactions
    do i_r = 1, n_prct
      ti_sp1 = im_photo_all(1,i_r); ti_sp_c1 = im_all_chem(ti_sp1)
      ti_sp2 = im_photo_all(2,i_r); ti_sp_c2 = im_all_chem(ti_sp2)
      ti_sp3 = im_photo_all(3,i_r); ti_sp_c3 = im_all_chem(ti_sp3)
      ti_sp4 = im_photo_all(4,i_r); ti_sp_c4 = im_all_chem(ti_sp4)
      ti_sp5 = im_photo_all(5,i_r); ti_sp_c5 = im_all_chem(ti_sp5)
      rate_prct(i_r,i_z) = freq_prct(i_r,i_z) * den(i_z,ti_sp1)
      ls_p(i_z,ti_sp1) = ls_p(i_z,ti_sp1) + rate_prct(i_r,i_z)
      pr_p(i_z,ti_sp2) = pr_p(i_z,ti_sp2) + rate_prct(i_r,i_z)
      pr_p(i_z,ti_sp3) = pr_p(i_z,ti_sp3) + rate_prct(i_r,i_z)
      pr_p(i_z,ti_sp4) = pr_p(i_z,ti_sp4) + rate_prct(i_r,i_z)
      pr_p(i_z,ti_sp5) = pr_p(i_z,ti_sp5) + rate_prct(i_r,i_z)
    end do
    
    ! Electron reactions
    do i_r = 1, n_erct
      ti_sp1 = im_e_all(1,i_r)
      ti_sp2 = im_e_all(2,i_r)
      ti_sp3 = im_e_all(3,i_r)
      ti_sp4 = im_e_all(4,i_r)
      ti_sp5 = im_e_all(5,i_r)
      rate_erct(i_r,i_z) = freq_erct(i_r,i_z) * den(i_z,ti_sp1)
      ls_e(i_z,ti_sp1) = ls_e(i_z,ti_sp1) + rate_erct(i_r,i_z)
      pr_e(i_z,ti_sp2) = pr_e(i_z,ti_sp2) + rate_erct(i_r,i_z)
      pr_e(i_z,ti_sp3) = pr_e(i_z,ti_sp3) + rate_erct(i_r,i_z)
      pr_e(i_z,ti_sp4) = pr_e(i_z,ti_sp4) + rate_erct(i_r,i_z)
      pr_e(i_z,ti_sp5) = pr_e(i_z,ti_sp5) + rate_erct(i_r,i_z)
    end do
  end do
  pr = pr_c + pr_p + pr_e
  ls = ls_c + ls_p + ls_e
  
  ! set very small values to zero
  where (rate_prct < 1.0E-99_wp) rate_prct = zero
  where (rate_erct < 1.0E-99_wp) rate_erct = zero
  where (rate_rct < 1.0E-99_wp) rate_rct = zero
  where (pr_c < 1.E-99_wp) pr_c = zero
  where (ls_c < 1.E-99_wp) ls_c = zero
  where (pr_p < 1.E-99_wp) pr_p = zero
  where (ls_p < 1.E-99_wp) ls_p = zero
  where (pr_e < 1.E-99_wp) pr_e = zero
  where (ls_e < 1.E-99_wp) ls_e = zero
  where (pr < 1.E-99_wp) pr = zero
  where (ls < 1.E-99_wp) ls = zero
  
  !----------------------------------------------------------------------------
  !  Read output settings
  !----------------------------------------------------------------------------

  open(newunit=fid_settings, file='output.settings', &
    status='old', action='read')
  
  read(fid_settings, *) ts_header, out_settings(1) ! atm1D.csv
  read(fid_settings, *) ts_header, out_settings(2) ! diff.csv
  read(fid_settings, *) ts_header, out_settings(3) ! solar_flux.out
  read(fid_settings, *) ts_header, out_settings(4) ! eflux.out
  read(fid_settings, *) ts_header, out_settings(5) ! balance.out
  read(fid_settings, *) ts_header, out_settings(6) ! ratecoeff.out
  read(fid_settings, *) ts_header, out_settings(7) ! chemrates.out
  read(fid_settings, *) ts_header, out_settings(8) ! colrates.out
  read(fid_settings, *) ts_header, out_settings(9) ! photorates.out
  read(fid_settings, *) ts_header, out_settings(10) ! COphotofreqB.out
  read(fid_settings, *) ts_header, out_settings(11) ! pcolrates.out
  read(fid_settings, *) ts_header, out_settings(12) ! elerates.out
  read(fid_settings, *) ts_header, out_settings(13) ! ecolrates.out
  read(fid_settings, *) ts_header, out_settings(14) ! species-specific files
  
  output_dir = '../runs/'//trim(runID)//'/output/'
  
  !----------------------------------------------------------------------------
  !  atm1D.out
  !----------------------------------------------------------------------------

  open(newunit=fid_atm1, file=trim(output_dir)//'atm1D.out', &
    status='unknown', action='write')

  write(fid_atm1, '(2I6)') n_z, n_sp
  write(fid_atm1, '("SPECIES")')
  write(fid_atm1, '(10(2X, A12, 1X))') (sp_list(i_sp), i_sp=1, n_sp)
  write(fid_atm1, '("ALTITUDE (km)")')
  write(fid_atm1, '(10ES15.7)') (cm_to_km * z(i_z), i_z=1, n_z)
  write(fid_atm1, '("RADIUS (km)")')
  write(fid_atm1, '(10ES15.7)') (cm_to_km * rz(i_z), i_z=1, n_z)
  write(fid_atm1, '("GRAVITY (cm s^-2)")')
  write(fid_atm1, '(10ES15.7)') (grv(i_z), i_z=1, n_z)
  write(fid_atm1, '("NEUTRAL TEMPERATURE (K)")')
  write(fid_atm1, '(10ES15.7)') (Tn(i_z), i_z=1, n_z)
  write(fid_atm1, '("ELECTRON TEMPERATURE (K)")')
  write(fid_atm1, '(10ES15.7)') (Te(i_z), i_z=1, n_z)
  write(fid_atm1, '("PRESSURE (dyne cm^-2 = 0.1 Pa)")')
  write(fid_atm1, '(10ES15.7)') (prs(i_z), i_z=1, n_z)
  write(fid_atm1, '("MASS DENSITY (g cm^-3)")')
  write(fid_atm1, '(10ES15.7)') (rho(i_z), i_z=1, n_z)
  write(fid_atm1, '("MEAN MOLECULAR WEIGHT (amu)")')
  write(fid_atm1, '(10ES15.7)') (mass(i_z), i_z=1, n_z)
  write(fid_atm1, '("EDDY COEFFICIENT (cm^2 s^-1)")')
  write(fid_atm1, '(10ES15.7)') (eK(i_z), i_z=1, n_z)
  write(fid_atm1, '("TOTAL DENSITY (cm^-3)")')
  write(fid_atm1, '(10ES15.7)') (den(i_z,0), i_z=1, n_z)
  do i_sp = 1, n_sp
    write(fid_atm1, '(A12,F5.2)') sp_list(i_sp)
    write(fid_atm1, '(10ES15.7)') (den(i_z,i_sp), i_z=1, n_z)
  end do

  close(unit=fid_atm1)

  !----------------------------------------------------------------------------
  !  atm1D.csv
  !----------------------------------------------------------------------------

  if (out_settings(1) == 'on ') then
    ! set column headings
    allocate(heads(n_sp+10))
    heads(1) = 'Alt (km)'
    heads(2) = 'Radius (km)'
    heads(3) = 'Grav (cm/s2)'
    heads(4) = 'Tn (K)'
    heads(5) = 'Te (K)'
    heads(6) = 'P (ubar)'
    heads(7) = 'rho (gm/cm3)'
    heads(8) = 'MMW (amu)'
    heads(9) = 'Kzz (cm2/s)'
    heads(10) = 'Ntot (cm-3)'
    do concurrent (i_sp = 1:n_sp)
      heads(10+i_sp) = sp_list(i_sp)
    end do
    
    ! write atmosphere to file
    open(newunit=fid_atm2, file=trim(output_dir)//'atm1D.csv', &
      status='unknown', action='write')
    write(fid_atm2, '(A12, 500(", ", A12))') (heads(i_h), i_h=1, 10+n_sp)
    do i_z = 1, n_z
      write(fid_atm2, '(ES15.7, 500(", ", ES15.7))') &
        cm_to_km * z(i_z), cm_to_km * rz(i_z), grv(i_z), Tn(i_z), Te(i_z), &
        prs(i_z), rho(i_z), mass(i_z), eK(i_z), den(i_z,0), &
        (vmr(i_z,i_sp), i_sp=1, n_sp)
    end do
    close(unit=fid_atm2)
    
    deallocate(heads)
  end if

  !----------------------------------------------------------------------------
  !  Diffusion coefficients
  !----------------------------------------------------------------------------

  if (out_settings(2) == 'on ') then
    ! set column headings
    allocate(heads(2*n_sp+3))
    heads(1) = 'Alt (km)'
    heads(2) = 'Kzz (cm^2/s)'
    do concurrent (i_sp = 1:n_sp)
      heads(2+i_sp) = sp_list(i_sp)
    end do
    heads(3+n_sp) = 'Ht (km)'
    do concurrent (i_sp = 1:n_sp)
      heads(3+n_sp+i_sp) = sp_list(i_sp)
    end do

    ! write diffusion coefficients and scale heights to file
    open(newunit=fid_diff, file=trim(output_dir)//'diff.csv', &
      status='unknown', action='write')
    write(fid_diff, '(A12, 500(", ", A12))') &
      (heads(i_h), i_h=1, 2*n_sp+3)
    ! 500 for being >> n_sp
    do i_z = 1, n_z
      write(fid_diff, '(ES15.7, 500(", ", ES15.7))') &
        cm_to_km * z(i_z), eK(i_z), &
        (Df(i_z,i_sp), i_sp=1, n_sp), &
        (cm_to_km * Ht(i_z,i_sp), i_sp=0, n_sp)
    end do
    close(unit=fid_diff)
    
    deallocate(heads)
  end if
  
  !----------------------------------------------------------------------------
  !  Solar fluxes
  !----------------------------------------------------------------------------

  if (out_settings(3) == 'on ') then
    open(newunit=fid_solar, file=trim(output_dir)//'solar_flux.out', &
      status='unknown', action='write')

    n_wave = n_waveA + n_waveB_lres + n_waveC - 2
    ! -2 because there is overlap at edge bins
    ! write header
    write(fid_solar, '(2I6)') n_wave, n_z
    
    ! write wavelength scale
    write(fid_solar, '("Wavelength (Angstroms)")')
    ! combine wavelength scales
    allocate(wave(n_wave))
    wave(1:n_waveA) = waveA
    wave(n_waveA+1:n_waveA+n_waveB_lres-2) = waveB_lres(2:n_waveB_lres-1)
    wave(n_waveA+n_waveB_lres-1:n_wave) = waveC
    write(fid_solar, '(10ES11.3)') (wave(i_wave), i_wave=1, n_wave)
    
    ! --- Write solar fluxes -----------------------------------------
    allocate(solar_flux(n_wave, n_z))
    solar_flux = zero
    do i_z = 1, n_z
      write(fid_solar, '(F9.1)') cm_to_km * z(i_z)
      ! region A
      do i_wave = 1, n_waveA
         solar_flux(i_wave,i_z) = sol_fluxA(i_wave) * trnA(i_wave,i_z)
      end do
      ! region B (sum over high-resolution bins)
      do i_wave = 1, n_waveB_lres-2
         do i_waveB = (i_wave-1)*200+1, i_wave*200
            solar_flux(n_waveA+i_wave,i_z) = solar_flux(n_waveA+i_wave,i_z) &
              + sol_fluxB(i_waveB) * trnB(i_waveB,i_z)
         end do
      end do
      ! region C
      do i_wave = 1, n_waveC
         solar_flux(n_waveA+n_waveB_lres-2+i_wave,i_z) = &
          sol_fluxC(i_wave) * trnC(i_wave,i_z)
      end do
      
      write(fid_solar, '(10ES11.3)') (solar_flux(i_wave,i_z), i_wave=1, n_wave)
    end do
    close(unit=fid_solar)
    
    deallocate(wave, solar_flux)
  end if
  
  !----------------------------------------------------------------------------
  !  Electron fluxes
  !----------------------------------------------------------------------------

  if (out_settings(4) == 'on ') then
    open(newunit=fid_eflux, file=trim(output_dir)//'eflux.out', &
      status='unknown', action='write')
    ! write header
    write(fid_eflux, '(2I6)') n_e_enrg, n_z
    ! write energy grid
    write(fid_eflux, '("ENERGY GRID (eV)")')
    write(fid_eflux, '(10ES11.3)') (e_enrg(i_enrg), i_enrg=1, n_e_enrg)
    ! write energy flux
    do i_z = 1, n_z
      write(fid_eflux, '(F11.3)') cm_to_km * z(i_z)
      write(fid_eflux, '(10ES11.3)') (eflux(i_z,i_enrg), i_enrg=1, n_e_enrg)
    end do
    close(unit=fid_eflux)
  end if

  !----------------------------------------------------------------------------
  !  Rate coefficients
  !----------------------------------------------------------------------------

  if (out_settings(6) == 'on ') then
    open(newunit=fid_ratecoeff, file=trim(output_dir)//'ratecoeff.out', &
      status='unknown', action='write')
    ! write header
    write(fid_ratecoeff, '(2I6)') n_rct, n_z
    ! write altitudes
    write(fid_ratecoeff, '("ALTITUDE (km)")')
    write(fid_ratecoeff, '(10ES11.3)') (cm_to_km * z(i_z), i_z=1, n_z)
    ! write rate coefficients
    do i_r = 1, n_rct
      write(fid_ratecoeff, '(A)') ctitle(i_r)
      write(fid_ratecoeff, '(10ES11.3)') (rk(i_r,i_z), i_z=1, n_z)
    end do
    close(unit=fid_ratecoeff)
  end if

  !----------------------------------------------------------------------------
  !  Reaction rates
  !----------------------------------------------------------------------------

  if (out_settings(7) == 'on ') then
    open(newunit=fid_chemrates, file=trim(output_dir)//'chemrates.out', &
      status='unknown', action='write')
    ! write header
    write(fid_chemrates, '(2I6)') n_rct, n_z
    ! write altitudes
    write(fid_chemrates, '("ALTITUDE (km)")')
    write(fid_chemrates, '(10ES11.3)') (cm_to_km * z(i_z), i_z=1, n_z)
    ! write rates
    do i_r = 1, n_rct
      write(fid_chemrates, '(A)') ctitle(i_r)
      write(fid_chemrates, '(10ES11.3)') (rate_rct(i_r,i_z), i_z=1, n_z)
    end do
    close(unit=fid_chemrates)
  end if

  !----------------------------------------------------------------------------
  !  Column-integrated reaction rates
  !----------------------------------------------------------------------------
  
  if (out_settings(8) == 'on ') then
    ! set column headings
    allocate(heads(4))
    heads(1) = 'Above aerosols'
    heads(2) = 'Below aerosols'
    heads(3) = 'Total'
    heads(4) = 'Reaction'
    
    allocate(rate_col(n_rct,3))
    open(newunit=fid_colrates, file=trim(output_dir)//'colrates.out', &
      status='unknown', action='write')
    do i_r = 1, n_rct
      rate_col(i_r,1) = col_integrate(rate_rct(i_r,:), iz_bot, n_z)
      rate_col(i_r,2) = col_integrate(rate_rct(i_r,:), 1, iz_bot)
      rate_col(i_r,3) = rate_col(i_r,1) + rate_col(i_r,2)
    end do
    
    allocate(t_sort_order(n_rct))
    call heapsort(rate_col(:,3), t_sort_order)
    write(fid_colrates, '(4A15)') (heads(i_h), i_h=1, 4)
    do i_r = n_rct, 1, -1
      write(fid_colrates, '(I4, 3(2X, ES11.3), 2X, A73)') &
        t_sort_order(i_r), rate_col(t_sort_order(i_r),1), &
        rate_col(t_sort_order(i_r),2), rate_col(t_sort_order(i_r),3), &
        ctitle(t_sort_order(i_r))
    end do
    close(unit=fid_colrates)
    
    deallocate(heads, rate_col, t_sort_order)
  end if
  
  !----------------------------------------------------------------------------
  !  Photo reaction rates
  !----------------------------------------------------------------------------

  if (out_settings(9) == 'on ') then
    open(newunit=fid_photorates, file=trim(output_dir)//'photorates.out', &
      status='unknown', action='write')

    write(fid_photorates, '(2I6)') n_prct, n_z
    ! write altitudes
    write(fid_photorates, '("ALTITUDE (km)")')
    write(fid_photorates, '(10ES11.3)') (cm_to_km * z(i_z), i_z=1, n_z)
    ! write rates
    do i_r = 1, n_prct
      write(fid_photorates, '(A)') ptitle(i_r)
      write(fid_photorates, '(10ES11.3)') (rate_prct(i_r, i_z), i_z=1, n_z)
     end do
    close(unit=fid_photorates)

  end if

  !----------------------------------------------------------------------------
  !  CO photodissociation rates in spectral region B
  !----------------------------------------------------------------------------

  if (out_settings(10) == 'on ') then
    open(newunit=fid_COphotofreqB, file=trim(output_dir)//'COphotofreqB.out', &
      status='unknown', action='write')

    write(fid_COphotofreqB, '(1I6)') n_waveB/100, n_z
    ! write wavelength scale
    write(fid_COphotofreqB, '("Wavelength (Angstroms)")')
    write(fid_COphotofreqB, '(10ES11.3)') &
      ((waveB((i_wave-1)*100+50) + waveB((i_wave-1)*100+51)) / 2, &
      i_wave=1, n_waveB/100)
    
    allocate(COphotorate(n_waveB/100, n_z))
    do i_z = 1, n_z
      write(fid_COphotofreqB, '(F11.3)') cm_to_km * z(i_z)
      do concurrent (i_wave = 1:n_waveB/100)
        COphotorate(i_wave,i_z) = dot_product( &
          prtB((i_wave-1)*100+1:i_wave*100,1,22), &
          trnB((i_wave-1)*100+1:i_wave*100,i_z))
      end do
      write(fid_COphotofreqB, '(10ES11.3)') &
        (COphotorate(i_wave,i_z), i_wave=1, n_waveB/100)
    end do

    close(unit=fid_COphotofreqB)
    
    deallocate(COphotorate)
  end if
  
  !----------------------------------------------------------------------------
  !  Column-integrated photo reaction rates
  !----------------------------------------------------------------------------
  
  if (out_settings(11) == 'on ') then
    ! set column headings
    allocate(heads(4))
    heads(1) = 'Above aerosols'
    heads(2) = 'Below aerosols'
    heads(3) = 'Total'
    heads(4) = 'Reaction'
    
    allocate(rate_col(n_prct,3))
    open(newunit=fid_pcolrates, file=trim(output_dir)//'pcolrates.out', &
      status='unknown', action='write')
    do i_r = 1, n_prct
      rate_col(i_r,1) = col_integrate(rate_prct(i_r,:), iz_bot, n_z)
      rate_col(i_r,2) = col_integrate(rate_prct(i_r,:), 1, iz_bot)
      rate_col(i_r,3) = rate_col(i_r,1) + rate_col(i_r,2)
    end do
    
    allocate(t_sort_order(n_prct))
    call heapsort(rate_col(:,3), t_sort_order)
    write(fid_pcolrates, '(4A15)') (heads(i_h), i_h=1, 4)
    do i_r = n_prct, 1, -1
      write(fid_pcolrates,'(I4, 3(2X, ES11.3), 2X, A87)') &
        t_sort_order(i_r), rate_col(t_sort_order(i_r),1), &
        rate_col(t_sort_order(i_r),2), rate_col(t_sort_order(i_r),3), &
        ptitle(t_sort_order(i_r))
    end do
    close(unit=fid_pcolrates)
    
    deallocate(heads, rate_col, t_sort_order)
  end if
  
  !----------------------------------------------------------------------------
  !  Electron reaction rates
  !----------------------------------------------------------------------------
  
  if (out_settings(12) == 'on ') then
    open(newunit=fid_elerates, file=trim(output_dir)//'elerates.out', &
      status='unknown', action='write')
    ! write header
    write(fid_elerates, '(2I6)') n_erct, n_z
    ! write altitudes
    write(fid_elerates, '("ALTITUDE (km)")')
    write(fid_elerates, '(10ES11.3)') (cm_to_km * z(i_z), i_z=1, n_z)
    do i_r = 1, n_erct
      write(fid_elerates, '(A)') etitle(i_r)
      write(fid_elerates, '(10ES11.3)') (rate_erct(i_r,i_z), i_z=1, n_z)
    end do
    close(unit=fid_elerates)
  end if

  !----------------------------------------------------------------------------
  !  Column-integrated electron reaction rates
  !----------------------------------------------------------------------------
  
  if (out_settings(13) == 'on ') then
    ! set column headings
    allocate(heads(2))
    heads(1) = 'Total'
    heads(2) = 'Reaction'
    
    allocate(rate_col(n_erct,1))
    open(newunit=fid_ecolrates, file=trim(output_dir)//'ecolrates.out', &
      status='unknown', action='write')
    do i_r = 1, n_erct
      rate_col(i_r, 1) = col_integrate(rate_erct(i_r,:), 1, n_z)
    end do

    allocate(t_sort_order(n_erct))
    call heapsort(rate_col(:,1), t_sort_order)
    write(fid_ecolrates, '(8X, 2A15)') (heads(i_h), i_h=1, 2)
    do i_r = n_erct, 1, -1
      write(fid_ecolrates, '(I4, 2X, ES11.3, 2X, A87)') &
        t_sort_order(i_r), rate_col(t_sort_order(i_r), 1), &
        etitle(t_sort_order(i_r))
    end do
    close(unit=fid_ecolrates)
    
    deallocate(heads, rate_col, t_sort_order)
  end if

  !----------------------------------------------------------------------------
  !  Column-integrated quantities for individual species
  !----------------------------------------------------------------------------
  
  if ((out_settings(5) == 'on ') .or. (out_settings(14) == 'on ')) then
    allocate(pr_c_col(n_sp), ls_c_col(n_sp))
    allocate(pr_p_col(n_sp), ls_p_col(n_sp))
    allocate(pr_e_col(n_sp), ls_e_col(n_sp))
    allocate(pr_col(n_sp), ls_col(n_sp))
    allocate(flux_top_col(n_sp), flux_bot_col(n_sp))
    allocate(bal_col(n_sp), den_col(n_sp), bal_time(n_sp))
    do i_sp = 1, n_sp
      pr_c_col(i_sp) = col_integrate(pr_c(:, i_sp), 1, n_z)
      ls_c_col(i_sp) = col_integrate(ls_c(:, i_sp), 1, n_z)
      pr_p_col(i_sp) = col_integrate(pr_p(:, i_sp), 1, n_z)
      ls_p_col(i_sp) = col_integrate(ls_p(:, i_sp), 1, n_z)
      pr_e_col(i_sp) = col_integrate(pr_e(:, i_sp), 1, n_z)
      ls_e_col(i_sp) = col_integrate(ls_e(:, i_sp), 1, n_z)
      pr_col(i_sp) = pr_p_col(i_sp) + pr_e_col(i_sp) + pr_c_col(i_sp)
      ls_col(i_sp) = ls_p_col(i_sp) + ls_e_col(i_sp) + ls_c_col(i_sp)
      den_col(i_sp) = col_integrate(den(:,i_sp), 1, n_z)
      ! normalize fluxes to surface
      flux_top_col(i_sp) = flux(n_z,i_sp) * (rz(n_z) / rPlanet)**two
      flux_bot_col(i_sp) = flux(1,i_sp) * (rz(1) / rPlanet)**two
      bal_col(i_sp) = pr_col(i_sp) - ls_col(i_sp) &
        + flux_bot_col(i_sp) - flux_top_col(i_sp)
      bal_time(i_sp) = den_col(i_sp) / abs(bal_col(i_sp))
    end do
  
    ! --- Column-integrated balances for each species -------------------------
    if (out_settings(5) == 'on ') then
      open(newunit=fid_balance, file=trim(output_dir)//'balance.out', &
        status='unknown', action='write')
      
      ! headings
      write(fid_balance, '("Species      Status  Density      ", &
        "Bot Flux     Top Flux     ", &
        "Production   Loss         Balance      Time Const")')
      
      allocate(t_sort_order(n_sp))
      call heapsort(abs(bal_col), t_sort_order)
      do i = n_sp, 1, -1    ! descending order
        i_sp = t_sort_order(i)
        if (istat(i_sp) > 0) then ! write only active species
          write(fid_balance, '(A12, I3, 7(2X, ES11.3))') &
            sp_list(i_sp), istat(i_sp), den_col(i_sp), &
            flux_bot_col(i_sp), flux_top_col(i_sp), &
            pr_col(i_sp), ls_col(i_sp), bal_col(i_sp), bal_time(i_sp)
        end if
      end do
      do i = n_sp, 1, -1    ! descending order
        i_sp = t_sort_order(i)
        if (istat(i_sp) <= 0) then ! write the other species
          write(fid_balance, '(A12, I3, 7(2X, ES11.3))') &
            sp_list(i_sp), istat(i_sp), den_col(i_sp), &
            flux_bot_col(i_sp), flux_top_col(i_sp), &
            pr_col(i_sp), ls_col(i_sp), bal_col(i_sp), bal_time(i_sp)
        end if
      end do
      close(unit=fid_balance)
      deallocate(t_sort_order)
    end if
  
    ! --- Species specific files ----------------------------------------------
    if (out_settings(14) == 'on ') then
      
      ! set headings
      allocate(heads(14))
      heads(1) = 'Alt'
      heads(2) = 'Density'
      heads(3) = 'VMR'
      heads(4) = 'Prod (Chem)'
      heads(5) = 'Loss (Chem)'
      heads(6) = 'Prod (Photo)'
      heads(7) = 'Loss (Photo)'
      heads(8) = 'Prod (Electron)'
      heads(9) = 'Loss (Electron)'
      heads(10) = 'Net Prod'
      heads(11) = 'Net Loss'
      heads(12) = 'Flux'
      heads(13) = 'div(flux)'
      heads(14) = 'Balance'
      
      do i = 1, size(im_chem_all)
        i_sp = im_chem_all(i)
        open(newunit=fid_species, &
          file=trim(output_dir)//'molecules/'//trim(sp_list(i_sp))//'.out', &
          status='unknown', action='write')
        
        ! column quantities
        write(fid_species, '(ES11.3," | Column Production (Chemical)")') &
          pr_c_col(i_sp)
        write(fid_species, '(ES11.3," | Column Loss (Chemical)")') &
          ls_c_col(i_sp)
        write(fid_species, '(ES11.3," | Column Production (Photo)")') &
          pr_p_col(i_sp)
        write(fid_species, '(ES11.3," | Column Loss (Photo)")') &
          ls_p_col(i_sp)
        write(fid_species, '(ES11.3," | Column Production (Electron)")') &
          pr_e_col(i_sp)
        write(fid_species, '(ES11.3," | Column Loss (Electron)")') &
          ls_e_col(i_sp)
        write(fid_species, '(ES11.3," | Column Net Production")') &
          pr_col(i_sp)
        write(fid_species, '(ES11.3," | Column Net Loss")') &
          ls_col(i_sp)
        write(fid_species, '(ES11.3," | Escape Flux")') &
          flux_top_col(i_sp)
        write(fid_species, '(ES11.3," | Surface Flux")') &
          flux_bot_col(i_sp)
        write(fid_species, '(ES11.3," | Net Balance")') &
          bal_col(i_sp)
        
        ! per altitude
        write(fid_species, '(14A15)') (heads(i_h), i_h=1, 14)
        do i_z = 1, n_z
          write(fid_species, '(F15.1, 13(4X, ES11.3))') &
            cm_to_km * z(i_z), den(i_z,i_sp), den(i_z,i_sp) / den(i_z, 0), &
            pr_c(i_z,i_sp), ls_c(i_z,i_sp), pr_p(i_z,i_sp), ls_p(i_z,i_sp), &
            pr_e(i_z,i_sp), ls_e(i_z,i_sp), pr(i_z,i_sp), ls(i_z,i_sp), &
            flux(i_z,i_sp), div_flux(i_z,i_sp),  &
            pr(i_z,i_sp) - ls(i_z,i_sp) - div_flux(i_z,i_sp)
        end do
        close(unit=fid_species)
      end do
      
      do i = 1, n_diff
        i_sp = im_diff_all(i)
        open(newunit=fid_species, &
          file=trim(output_dir)//'molecules/'//trim(sp_list(i_sp))//'.out', &
          status='unknown', action='write')
        
        ! column quantities
        write(fid_species, '(ES11.3," | Column Production (Chemical)")') &
          pr_c_col(i_sp)
        write(fid_species, '(ES11.3," | Column Loss (Chemical)")') &
          ls_c_col(i_sp)
        write(fid_species, '(ES11.3," | Column Production (Photo)")') &
          pr_p_col(i_sp)
        write(fid_species, '(ES11.3," | Column Loss (Photo)")') &
          ls_p_col(i_sp)
        write(fid_species, '(ES11.3," | Column Production (Electron)")') &
          pr_e_col(i_sp)
        write(fid_species, '(ES11.3," | Column Loss (Electron)")') &
          ls_e_col(i_sp)
        write(fid_species, '(ES11.3," | Column Net Production")') &
          pr_col(i_sp)
        write(fid_species, '(ES11.3," | Column Net Loss")') &
          ls_col(i_sp)
        write(fid_species, '(ES11.3," | Escape Flux")') &
          flux_top_col(i_sp)
        write(fid_species, '(ES11.3," | Surface Flux")') &
          flux_bot_col(i_sp)
        write(fid_species, '(ES11.3," | Net Balance")') &
          bal_col(i_sp)
        
        ! per altitude
        write(fid_species, '(14A15)') (heads(i_h), i_h=1, 14)
        do i_z = 1, n_z
          write(fid_species, '(F15.1, 13(4X, ES11.3))') &
            cm_to_km * z(i_z), den(i_z,i_sp), den(i_z,i_sp) / den(i_z, 0), &
            pr_c(i_z,i_sp), ls_c(i_z,i_sp), pr_p(i_z,i_sp), ls_p(i_z,i_sp), &
            pr_e(i_z,i_sp), ls_e(i_z,i_sp), pr(i_z,i_sp), ls(i_z,i_sp), &
            flux(i_z,i_sp), div_flux(i_z,i_sp),  &
            pr(i_z,i_sp) - ls(i_z,i_sp) - div_flux(i_z,i_sp)
        end do
        close(unit=fid_species)
      end do
      
      ! --- Oxygen ------------------------------------------------------------
      allocate(den_oxy(n_z))
      allocate(pr_c_oxy(n_z), ls_c_oxy(n_z))
      allocate(pr_p_oxy(n_z), ls_p_oxy(n_z))
      allocate(pr_e_oxy(n_z), ls_e_oxy(n_z))
      allocate(pr_oxy(n_z), ls_oxy(n_z))
      allocate(flux_oxy(n_z), div_flux_oxy(n_z))
      do concurrent (i_z = 1:n_z)
        den_oxy(i_z) = dot_product(noxy, den(i_z, 1:n_sp))
        pr_p_oxy(i_z) = dot_product(noxy, pr_p(i_z, :))
        ls_p_oxy(i_z) = dot_product(noxy, ls_p(i_z, :))
        pr_e_oxy(i_z) = dot_product(noxy, pr_e(i_z, :))
        ls_e_oxy(i_z) = dot_product(noxy, ls_e(i_z, :))
        pr_c_oxy(i_z) = dot_product(noxy, pr_c(i_z, :))
        ls_c_oxy(i_z) = dot_product(noxy, ls_c(i_z, :))
        pr_oxy(i_z) = dot_product(noxy, pr(i_z, :))
        ls_oxy(i_z) = dot_product(noxy, ls(i_z, :))
        flux_oxy(i_z) = dot_product(noxy, flux(i_z, :))
        div_flux_oxy(i_z) = dot_product(noxy, div_flux(i_z, :))
      end do
      
      pr_c_col_oxy = col_integrate(pr_c_oxy, 1, n_z)
      ls_c_col_oxy = col_integrate(ls_c_oxy, 1, n_z)
      pr_p_col_oxy = col_integrate(pr_p_oxy, 1, n_z)
      ls_p_col_oxy = col_integrate(ls_p_oxy, 1, n_z)
      pr_e_col_oxy = col_integrate(pr_e_oxy, 1, n_z)
      ls_e_col_oxy = col_integrate(ls_e_oxy, 1, n_z)
      pr_col_oxy = pr_c_col_oxy + pr_p_col_oxy + pr_e_col_oxy
      ls_col_oxy = ls_c_col_oxy + ls_p_col_oxy + ls_e_col_oxy
      flux_top_col_oxy = flux_oxy(n_z) * (rz(n_z) / rPlanet)**two
      flux_bot_col_oxy = flux_oxy(1) * (rz(1) / rPlanet)**two
      bal_col_oxy = pr_col_oxy - ls_col_oxy &
        + flux_bot_col_oxy - flux_top_col_oxy

      open(newunit=fid_species, &
        file=trim(output_dir)//'molecules/O_atoms.out', &
        status='unknown', action='write')
      ! column quantities
      write(fid_species, '(ES11.3," | Column Production (Chemical)")') &
        pr_c_col_oxy
      write(fid_species, '(ES11.3," | Column Loss (Chemical)")') &
        ls_c_col_oxy
      write(fid_species, '(ES11.3," | Column Production (Photo)")') &
        pr_p_col_oxy
      write(fid_species, '(ES11.3," | Column Loss (Photo)")') &
        ls_p_col_oxy
      write(fid_species, '(ES11.3," | Column Production (Electron)")') &
        pr_e_col_oxy
      write(fid_species, '(ES11.3," | Column Loss (Electron)")') &
        ls_e_col_oxy
      write(fid_species, '(ES11.3," | Column Net Production")') &
        pr_col_oxy
      write(fid_species, '(ES11.3," | Column Net Loss")') &
        ls_col_oxy
      write(fid_species, '(ES11.3," | Escape Flux")') &
        flux_top_col_oxy
      write(fid_species, '(ES11.3," | Surface Flux")') &
        flux_bot_col_oxy
      write(fid_species, '(ES11.3," | Net Balance")') &
        bal_col_oxy
      ! per altitude
      write(fid_species, '(14A15)') (heads(i_h), i_h=1, 14)
      do i_z = 1, n_z
        write(fid_species, '(F15.1, 13(4X, ES11.3))') &
          cm_to_km * z(i_z), den_oxy(i_z), den_oxy(i_z) / den(i_z, 0), &
          pr_c_oxy(i_z), ls_c_oxy(i_z), pr_p_oxy(i_z), ls_p_oxy(i_z), &
          pr_e_oxy(i_z), ls_e_oxy(i_z), pr_oxy(i_z), ls_oxy(i_z), &
          flux_oxy(i_z), div_flux_oxy(i_z), &
          pr_oxy(i_z) - ls_oxy(i_z) - div_flux_oxy(i_z)
      end do
      close(unit=fid_species)
      
      deallocate(den_oxy)
      deallocate(pr_c_oxy, ls_c_oxy, pr_p_oxy, ls_p_oxy, pr_e_oxy, ls_e_oxy)
      deallocate(pr_oxy, ls_oxy)
      deallocate(flux_oxy, div_flux_oxy)

      ! --- Hydrogen ----------------------------------------------------------
      allocate(den_hyd(n_z, 2))
      allocate(pr_c_hyd(n_z, 2), ls_c_hyd(n_z, 2))
      allocate(pr_p_hyd(n_z, 2), ls_p_hyd(n_z, 2))
      allocate(pr_e_hyd(n_z, 2), ls_e_hyd(n_z, 2))
      allocate(pr_hyd(n_z, 2), ls_hyd(n_z, 2))
      allocate(flux_hyd(n_z, 2), div_flux_hyd(n_z, 2))
      den_hyd = zero
      pr_c_hyd = zero; ls_c_hyd = zero
      pr_p_hyd = zero; ls_p_hyd = zero
      pr_e_hyd = zero; ls_e_hyd = zero
      pr_hyd = zero; ls_hyd = zero
      flux_hyd = zero; div_flux_hyd = zero
      do i_sp = 1, n_sp-1
        do concurrent (i_z = 1:n_z)
          if ((nhyd(i_sp) == 1) .or. (nhyd(i_sp) == 2)) then
            den_hyd(i_z, nhyd(i_sp)) = den_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * den(i_z, i_sp)
            flux_hyd(i_z, nhyd(i_sp)) = flux_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * flux(i_z, i_sp)
            div_flux_hyd(i_z, nhyd(i_sp)) = div_flux_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * div_flux(i_z, i_sp)
            pr_c_hyd(i_z, nhyd(i_sp)) = pr_c_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * pr_c(i_z, i_sp)
            ls_c_hyd(i_z, nhyd(i_sp)) = ls_c_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * ls_c(i_z, i_sp)
            pr_p_hyd(i_z, nhyd(i_sp)) = pr_p_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * pr_p(i_z, i_sp)
            ls_p_hyd(i_z, nhyd(i_sp)) = ls_p_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * ls_p(i_z, i_sp)
            pr_e_hyd(i_z, nhyd(i_sp)) = pr_e_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * pr_e(i_z, i_sp)
            ls_e_hyd(i_z, nhyd(i_sp)) = ls_e_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * ls_e(i_z, i_sp)
            pr_hyd(i_z, nhyd(i_sp)) = pr_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * pr(i_z, i_sp)
            ls_hyd(i_z, nhyd(i_sp)) = ls_hyd(i_z, nhyd(i_sp)) &
              + nhyd(i_sp) * ls(i_z, i_sp)
          end if
        end do
      end do
      
      do i = 1, 2
        pr_c_col_hyd(i) = col_integrate(pr_c_hyd(:,i), 1, n_z)
        ls_c_col_hyd(i) = col_integrate(ls_c_hyd(:,i), 1, n_z)
        pr_p_col_hyd(i) = col_integrate(pr_p_hyd(:,i), 1, n_z)
        ls_p_col_hyd(i) = col_integrate(ls_p_hyd(:,i), 1, n_z)
        pr_e_col_hyd(i) = col_integrate(pr_e_hyd(:,i), 1, n_z)
        ls_e_col_hyd(i) = col_integrate(ls_e_hyd(:,i), 1, n_z)
        pr_col_hyd(i) = pr_p_col_hyd(i) + pr_e_col_hyd(i) + pr_c_col_hyd(i)
        ls_col_hyd(i) = ls_p_col_hyd(i) + ls_e_col_hyd(i) + ls_c_col_hyd(i)
        flux_top_col_hyd(i) = flux_hyd(n_z, i) * (rz(n_z) / rPlanet)**two
        flux_bot_col_hyd(i) = flux_hyd(1, i) * (rz(1) / rPlanet)**two
        bal_col_hyd(i) = pr_col_hyd(i) - ls_col_hyd(i) &
          + flux_bot_col_hyd(i) - flux_top_col_hyd(i)
      end do
      
      ! write odd hydrogen
      i = 1
      open(newunit=fid_species, file=trim(output_dir)//'molecules/H_odd.out', &
        status='unknown', action='write')
      ! column quantities
      write(fid_species, '(ES11.3," | Column Production (Chemical)")') &
        pr_c_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Loss (Chemical)")') &
        ls_c_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Production (Photo)")') &
        pr_p_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Loss (Photo)")') &
        ls_p_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Production (Electron)")') &
        pr_e_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Loss (Electron)")') &
        ls_e_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Net Production")') &
        pr_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Net Loss")') &
        ls_col_hyd(i)
      write(fid_species, '(ES11.3," | Escape Flux")') &
        flux_top_col_hyd(i)
      write(fid_species, '(ES11.3," | Surface Flux")') &
        flux_bot_col_hyd(i)
      write(fid_species, '(ES11.3," | Net Balance")') &
        bal_col_hyd(i)
      ! per altitude
      write(fid_species, '(14A15)') (heads(i_h), i_h=1, 14)
      do i_z = 1, n_z
        write(fid_species, '(F15.1, 13(4X,ES11.3))') &
          cm_to_km * z(i_z), den_hyd(i_z, i), den_hyd(i_z, i) / den(i_z, 0), &
          pr_c_hyd(i_z, i), ls_c_hyd(i_z, i), &
          pr_p_hyd(i_z, i), ls_p_hyd(i_z, i), &
          pr_e_hyd(i_z, i), ls_e_hyd(i_z, i), &
          pr_hyd(i_z, i), ls_hyd(i_z, i), &
          flux_hyd(i_z, i), div_flux_hyd(i_z, i), &
          pr_hyd(i_z, i) - ls_hyd(i_z, i) - div_flux_hyd(i_z, i)
      end do
      close(unit=fid_species)

      ! write even hydrogen
      i = 2
      open(newunit=fid_species, file=trim(output_dir)//'molecules/H_odd.out', &
        status='unknown', action='write')
      ! column quantities
      write(fid_species, '(ES11.3," | Column Production (Chemical)")') &
        pr_c_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Loss (Chemical)")') &
        ls_c_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Production (Photo)")') &
        pr_p_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Loss (Photo)")') &
        ls_p_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Production (Electron)")') &
        pr_e_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Loss (Electron)")') &
        ls_e_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Net Production")') &
        pr_col_hyd(i)
      write(fid_species, '(ES11.3," | Column Net Loss")') &
        ls_col_hyd(i)
      write(fid_species, '(ES11.3," | Escape Flux")') &
        flux_top_col_hyd(i)
      write(fid_species, '(ES11.3," | Surface Flux")') &
        flux_bot_col_hyd(i)
      write(fid_species, '(ES11.3," | Net Balance")') &
        bal_col_hyd(i)
      ! per altitude
      write(fid_species, '(14A15)') (heads(i_h), i_h=1, 14)
      do i_z = 1, n_z
        write(fid_species, '(F15.1, 13(4X,ES11.3))') &
          cm_to_km * z(i_z), den_hyd(i_z, i), den_hyd(i_z, i) / den(i_z, 0), &
          pr_c_hyd(i_z, i), ls_c_hyd(i_z, i), &
          pr_p_hyd(i_z, i), ls_p_hyd(i_z, i), &
          pr_e_hyd(i_z, i), ls_e_hyd(i_z, i), &
          pr_hyd(i_z, i), ls_hyd(i_z, i), &
          flux_hyd(i_z, i), div_flux_hyd(i_z, i), &
          pr_hyd(i_z, i) - ls_hyd(i_z, i) - div_flux_hyd(i_z, i)
      end do
      deallocate(den_hyd)
      deallocate(pr_c_hyd, ls_c_hyd, pr_p_hyd, ls_p_hyd, pr_e_hyd, ls_e_hyd)
      deallocate(pr_hyd, ls_hyd)
      deallocate(flux_hyd, div_flux_hyd)
      
      deallocate(heads)
    end if
    
    deallocate(pr_c_col, ls_c_col, pr_p_col, ls_p_col, pr_e_col, ls_e_col)
    deallocate(pr_col, ls_col, flux_top_col, flux_bot_col)
    deallocate(bal_col, den_col, bal_time)
  end if
end subroutine write_output


real(kind=wp) function col_integrate(var, alt0, alt1) result(var_col)
! This function does column integration on altitude profile var between the
! altitudes indices alt0 and alt1.
  use types, only: wp => dp
  use constants
  use global_variables, only: rz, dr_mid, rPlanet
  implicit none
  real(wp), dimension(:), intent(in) :: var
  integer, intent(in) :: alt0, alt1
  integer :: i_z
  
  var_col = zero
  do i_z = alt0, alt1-1
    var_col = var_col + dr_mid(i_z+1) * half &
      * (var(i_z) * rz(i_z)**two + var(i_z+1) * rz(i_z+1)**two) / rPlanet**two
  end do  
end function col_integrate


end module
