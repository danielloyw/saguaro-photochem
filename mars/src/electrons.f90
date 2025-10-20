module electrons

contains

subroutine electrons_setup
! This subroutine sets up the calculation for electron impact reactions. 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : find_name, find_bin

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  implicit none
  ! file unit numbers for eimpact.dat
  integer :: fid
  ! maximum states for each species with cross sections
  integer :: n_state_max
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, allocatable, dimension(:,:,:) :: im_e_product_all
  ! number of reactions (excitation, dissociation, ionization)
  integer :: n_excite, n_dissoc, n_ioniz
  ! number of dissociation + ionization reactions
  integer :: n_sp_e

  ! loop variables
  integer :: i, i_sp, i_branch, i_state, i_enrg, i_p, i_r
  ! --- Temporary / dummy variables -------------------------------------------
  character(len=256) :: ts_line
  character(len=10) :: ts_state
  character(len=87), allocatable, dimension(:,:) :: ts_title
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  real(wp), allocatable, dimension(:) :: t_cs
  ! energy of excited / ionized state
  real(wp) :: t_enrg_state
  ! energy after ionization
  real(wp) :: t_enrg_e
  ! discretized energy after ionization
  real(wp) :: t_enrgdis_e
  ! discretized ionization energy
  real(wp) :: t_enrgdis_ion
  integer :: tn_isp, tn_ibin
  ! integer :: tn1, tn2

  !----------------------------------------------------------------------------
  !  Read electron impact cross sections
  !----------------------------------------------------------------------------

  open(newunit=fid, file='../data/electrons/eimpact.dat', &
    status='old', action='read')

  read(fid,*) n_e_enrg, n_sp_e_thk, n_sp_e_thn, n_state_max

  n_sp_e = n_sp_e_thk + n_sp_e_thn ! number of species
  
  ! read energy bins
  allocate(e_enrg(n_e_enrg), e_denrg(n_e_enrg))
  read(fid,'(A)') ts_line ! "MEAN ENERGY IN BINS"
  read(fid,*) (e_enrg(i_enrg), i_enrg=1,n_e_enrg)
  read(fid,'(A)') ts_line ! "WIDTH OF ENERGY BINS"
  read(fid,*) (e_denrg(i_enrg), i_enrg=1,n_e_enrg)
  
  ! read in cross sections
  allocate(n_states(3,n_sp_e_thk))
  allocate(n_branch_e(n_sp_e))
  allocate(enrg_state(n_state_max,n_sp_e_thk))
  allocate(ecs(n_e_enrg,n_state_max,n_sp_e_thk))
  !allocate(branch_ratio_el(n_e_enrg,n_state_max,n_sp_e))
  !allocate(crs_tot_inel(n_e_enrg,n_sp_e))
  
  allocate(im_e_reactant_all(n_sp_e))
  allocate(im_e_product_all(4,n_state_max,n_sp_e))
  allocate(ts_title(n_state_max,n_sp_e))
  allocate(t_cs(n_e_enrg))
  
  do i_sp = 1, n_sp_e_thk
    read(fid,'(A12, 3I4)') ts_species1, n_excite, n_dissoc, n_ioniz
    im_e_reactant_all(i_sp) = find_name(ts_species1, sp_list)

    ! check if species exists in sp_list
    if (im_e_reactant_all(i_sp) <= 0) then
      write(*,'("Error: Electron impact species not found: ", A12)') & 
        ts_species1
      stop
    end if
      
    n_states(1,i_sp) = n_excite
    n_states(2,i_sp) = n_dissoc
    n_states(3,i_sp) = n_ioniz
    n_branch_e(i_sp) = n_dissoc + n_ioniz
    
    ! read total cross sections
    read(fid,'(A)') ts_line    ! "TOTAL ELASTIC CROSS SECTION"
    read(fid,*) (t_cs(i_enrg), i_enrg=1,n_e_enrg)
    read(fid,'(A)') ts_line    ! "TOTAL INELASTIC CROSS SECTION"
    read(fid,*) (t_cs(i_enrg), i_enrg=1,n_e_enrg)
    
    ! read cross sections for excited states
    do i_state = 1, n_excite
      read(fid,'(A12, F9.3)') ts_state, enrg_state(i_state,i_sp)
      read(fid,*) (ecs(i_enrg,i_state,i_sp), i_enrg=1,n_e_enrg)
    end do
    
    ! read cross sections for dissociation & ionization states 
    do i_branch = 1, n_branch_e(i_sp)
      ! read state and ionization energy
      read(fid,'(A12, F9.3)') ts_state, enrg_state(n_excite+i_branch,i_sp)
      
      ! read reaction formula
      read(fid,'(A)') ts_line
      ts_title(i_branch,i_sp) = ts_line
      read(ts_line,'(5(A12, 3X), A12)') (ts_species2(i), i=1,6)
      
      ! link products to species list
      i_p = 0
      do i = 3, 6    ! products
        tn_isp = find_name(ts_species2(i), sp_list)
        if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          i_p = i_p + 1
          im_e_product_all(i_p,i_branch,i_sp) = tn_isp
        end if
      end do
      ! read cross sections
      read(fid,*) (ecs(i_enrg,n_excite+i_branch,i_sp), i_enrg=1,n_e_enrg)
    end do
  end do
  
  ! do i_sp = 1, n_sp_e_thn
    ! i_sp = i_sp + 1
    ! read(fid,'(A12, I4)') ts_species1, n_branch_e(i_sp)
    ! im_e_reactant_all(i_sp) = find_name(ts_species1, sp_list)
    
    ! ! check if species exists in sp_list
    ! if (im_e_reactant_all(i_sp) <= 0) then
      ! write(*,'("Error: Electron impact species not found: ", A12)') & 
        ! ts_species1
      ! stop
    ! end if
  
    ! ! read total cross sections
    ! read(fid,'(A)') ts_line
    ! read(fid,'(10ES11.3)') (crs_tot_inel(i_enrg,i_sp),i_enrg=1,n_e_enrg)
    
    ! do i_branch = 1, n_branch_e(i_sp)
      ! read(fid,'(A)') ts_line
      ! ts_title(i_branch,i_sp) = ts_line
      ! read(ts_line,'((5(A12,3X),A12,F6.1))') (ts_species2(j),j=1,6), rdum
      ! i_p = 0
      ! do i = 3, 6 ! products
        ! tn_isp = find_name(ts_species2(i), sp_list)
        ! if ((tn_isp > 0) .and. (tn_isp <= n_sp)) then
          ! i_p = i_p + 1
          ! im_e_product_all(i_p,i_branch,i_sp) = tn_isp
        ! end if
      ! end do
      ! nprdcts_el(i_branch,i_sp) = i_p
      ! read(fid,'(10ES11.3)') (branch_ratio_el(i_enrg,i_branch,i_sp), i_enrg=1,n_e_enrg)
    ! end do
  ! end do

  close(unit=fid)
  
  !----------------------------------------------------------------------------
  !  Assign global variables
  !----------------------------------------------------------------------------
  
  ! count number of electron reactions
  i_r = 0
  do i_sp = 1, n_sp_e
    do i_branch = 1, n_branch_e(i_sp)
      i_r = i_r + 1
   end do
  end do
  n_erct = i_r

  allocate(etitle(n_erct), im_e_all(5,n_erct))
  i_r = 0
  do i_sp = 1, n_sp_e
    do i_branch = 1, n_branch_e(i_sp)
      i_r = i_r + 1
      etitle(i_r) = ts_title(i_branch,i_sp)
      im_e_all(1,i_r) = im_e_reactant_all(i_sp)
      im_e_all(2:5,i_r) = im_e_product_all(1:4,i_branch,i_sp)
    end do
  end do

  ! calculate total cross sections for excitation/dissociation and ionization
  ! allocate(ecs_exc(n_e_enrg,n_sp_e_thk), ecs_ion(n_e_enrg,n_sp_e_thk))
  ! do i_sp = 1, n_sp_e_thk
    ! do i_enrg = 1, n_e_enrg
      ! tn2 = n_states(1,i_sp) + n_states(2,i_sp)
      ! ecs_exc(i_enrg,i_sp) = sum(ecs(i_enrg,1:tn2,i_sp))
      ! tn1 = tn2 + 1
      ! tn2 = tn2 + n_states(3,i_sp)
      ! ecs_ion(i_enrg,i_sp) = sum(ecs(i_enrg,tn1:tn2,i_sp))
    ! end do
  ! end do
  
  ! calculate cross sections by summing differential cross sections
  allocate(cs_ion(n_e_enrg,n_state_max,n_sp_e))
  allocate(cs_sec(n_e_enrg,n_state_max,n_sp_e))
  do i_sp = 1, n_sp_e_thk
    do i_state = 1, n_states(3,i_sp)
      do i_enrg = 1, n_e_enrg
        t_enrg_state = &
          enrg_state(n_states(1,i_sp)+n_states(2,i_sp)+i_state, i_sp)
        if (e_enrg(i_enrg) > t_enrg_state) then
        ! electron energy > ionization energy
          
          ! --- Calculate discretized ionization energy -----------------------
          ! In modeling the electron cascade, electrons within the initial bin
          ! are pushed up against the higher-energy side of the bin. When a
          ! reaction results in the electron falling into the lower half of the
          ! final bin, it is likely that the final bin is actually one lower
          ! without the "pushing up" in the initial bin. Electrons falling into
          ! the upper half of the final bin will get pushed up. 
          t_enrg_e = e_enrg(i_enrg) - t_enrg_state
          tn_ibin = find_bin(t_enrg_e, e_enrg, e_denrg)
          if (t_enrg_e < e_enrg(tn_ibin)) then
          ! total energy after ionization < bin center
            t_enrgdis_e = e_enrg(tn_ibin) - half * e_denrg(tn_ibin)
            tn_ibin = tn_ibin - 1
          else    ! total energy after ionization >= bin center
            ! then push electron energy
            t_enrgdis_e = e_enrg(tn_ibin) + half * e_denrg(tn_ibin)
          end if
          t_enrgdis_ion = e_enrg(i_enrg) - t_enrgdis_e
          ! -------------------------------------------------------------------
          
          if (tn_ibin == 0) then
          ! when exiting e energy < first bin center
            ! then all electrons in single bin and normalization factor = 1
            cs_ion(i_enrg,i_state,i_sp) = one 
            cs_sec(i_enrg,i_state,i_sp) = one
          else
            do i = 1, tn_ibin
              cs_ion(i_enrg,i_state,i_sp) = cs_ion(i_enrg,i_state,i_sp) &
                + difcs_sec(e_enrg(i_enrg) - e_enrg(i) - t_enrgdis_ion, &
                e_enrg(i_enrg), t_enrgdis_ion) & 
                * e_denrg(i)
              cs_sec(i_enrg,i_state,i_sp) = cs_sec(i_enrg,i_state,i_sp) &
                + difcs_sec(e_enrg(i), e_enrg(i_enrg), t_enrgdis_ion) &
                * e_denrg(i)
            end do
          end if
        end if
      end do
    end do
  end do

  ! allocate global variables
  allocate(esrc(n_prct,n_z,n_e_enrg))
  allocate(Sel(n_z,n_e_enrg))
  allocate(eflux(n_z,n_e_enrg))
  allocate(eph(n_erct,n_z))
  allocate(rpe(n_erct,n_z))
  ! allocate(pS(n_z))
  rpe = zero

  ! deallocate local arrays
  deallocate(im_e_product_all, ts_title, t_cs)
  
end subroutine electrons_setup


subroutine cal_electrons
! This subroutine calculates the electron impact reaction rates and
! omnidirectional electron flux, including the production of secondary
! electrons and the energy cascade from the electrons colliding with neutrals
! and the background of thermal electrons. The subroutine assumes no vertical
! transport and solves backwards for the mean electron intensity at each energy
! bin assuming zero input for the highest energy bin.  

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
  ! attenuation by neutral excitation, secondary electrons, ionization, total:
  ! unit=cm-1; dim=(final energy bin, initial energy bin)
  real(wp), allocatable, dimension(:,:) :: Ael_exc, Ael_sec, Ael_ion, Ael
  ! total excitation and ionization cross sections: unit=cm2;
  ! dim=(energy bin, e species #)
  real(wp), allocatable, dimension(:,:) :: sum_cs_exc, sum_cs_ion
  ! sigma_T * N in photoelectron notes Eq. 4: unit=cm-1; dim=(energy bin)
  real(wp), allocatable, dimension(:) :: sigmaN

  ! loop variables
  integer :: i, i_z, i_sp, i_branch, i_enrg, i_state1, i_state2, i_r
  ! --- Temporary / dummy variables -------------------------------------------
  integer :: tn_ibin
  ! energy after ionization
  real(wp) :: t_enrg_e
  ! discretized energy after ionization
  real(wp) :: t_enrgdis_e
  ! discretized ionization energy
  real(wp) :: t_enrgdis_ion
  ! electron cross section scaled to discretization scheme
  real(wp) :: t_ecs_scaled
  real(wp) :: t_sc
  ! terms in photoelectron notes Eq. 4 (alpha_j, alpha_(j+1)
  real(wp) :: t_aj, t_aj1

  !----------------------------------------------------------------------------
  ! Calculate electron flux at each energy and level
  !----------------------------------------------------------------------------
  
  allocate(Ael_exc(n_e_enrg,n_e_enrg), Ael_sec(n_e_enrg,n_e_enrg))
  allocate(Ael_ion(n_e_enrg,n_e_enrg), Ael(n_e_enrg,n_e_enrg))
  allocate(sum_cs_ion(n_e_enrg,n_sp_e_thk), sum_cs_exc(n_e_enrg,n_sp_e_thk))
  allocate(sigmaN(n_e_enrg))

  eflux = zero
  do i_z = iz_bot, n_z
    Ael = zero
    Ael_exc = zero
    Ael_sec = zero
    Ael_ion = zero

    if (sum(Sel(i_z,1:n_e_enrg)) > zero) then
      sum_cs_exc = zero
      sum_cs_ion = zero

      ! --- Calculate Attenuation ---------------------------------------------
      do i_enrg = 2, n_e_enrg ! lowest bin doesn't need further attenuation
        do i_sp = 1, n_sp_e_thk

          ! Excitation to discrete excitation and dissociation states
          do i_state1 = 1, n_states(1,i_sp) + n_states(2,i_sp)
            ! exiting energy of electron
            t_enrg_e = e_enrg(i_enrg) - enrg_state(i_state1,i_sp)
            if (t_enrg_e > zero) then    ! excitation occurs
              tn_ibin = find_bin(t_enrg_e, e_enrg, e_denrg)
              if (tn_ibin < i_enrg) then
              ! scale cross sections for discretization
              ! transition out of bin (Eq. (25) of photoelectron notes)
                t_ecs_scaled = ecs(i_enrg,i_state1,i_sp) &
                  * enrg_state(i_state1,i_sp) &
                  / (e_enrg(i_enrg) - e_enrg(tn_ibin))
                Ael_exc(tn_ibin,i_enrg) = Ael_exc(tn_ibin,i_enrg) &
                  + den(i_z,im_e_reactant_all(i_sp)) * t_ecs_scaled
              else if (tn_ibin == i_enrg) then
              ! transition within bin (Eq. (27) of photoelectron notes)
                t_ecs_scaled = ecs(i_enrg,i_state1,i_sp) &
                  * enrg_state(i_state1,i_sp) &
                  / (e_enrg(tn_ibin) - e_enrg(tn_ibin-1))
                Ael_exc(tn_ibin-1,i_enrg) = Ael_exc(tn_ibin-1,i_enrg) &
                  + den(i_z,im_e_reactant_all(i_sp)) * t_ecs_scaled
              end if
              sum_cs_exc(i_enrg,i_sp) = sum_cs_exc(i_enrg,i_sp) + t_ecs_scaled
            end if
          end do

          ! Ionization channels
          do i_state1 = 1, n_states(3,i_sp)
            i_state2 = n_states(1,i_sp) + n_states(2,i_sp) + i_state1
            ! exiting energy of electron
            t_enrg_e = e_enrg(i_enrg) - enrg_state(i_state2,i_sp)
            if (t_enrg_e > zero) then    ! ionization occurs
              tn_ibin = find_bin(t_enrg_e, e_enrg, e_denrg)
              if (t_enrg_e < (e_enrg(tn_ibin)) .and. (tn_ibin /= 1)) then
                t_enrgdis_e = e_enrg(tn_ibin) - half * e_denrg(tn_ibin)
                tn_ibin = tn_ibin - 1
              else
                t_enrgdis_e = e_enrg(tn_ibin) + half * e_denrg(tn_ibin)
              end if
              t_enrgdis_ion = e_enrg(i_enrg) - t_enrgdis_e
              do i = 1, tn_ibin
                if (i == tn_ibin) then
                  ! scale to conserve energy
                  t_sc = enrg_state(i_state2,i_sp) / t_enrgdis_ion
                else
                  t_sc = one
                end if
                ! weighted averaging to account for shape of differential
                ! cross sections: 
                ! sum(den*ecs*difcs_sec*e_denrg) / sum(difcs_sec*e_denrg)
                Ael_ion(i,i_enrg) = Ael_ion(i,i_enrg) &
                  + den(i_z,im_e_reactant_all(i_sp)) &
                  * ecs(i_enrg,i_state2,i_sp) &
                  * difcs_sec(e_enrg(i_enrg) - e_enrg(i) - t_enrgdis_ion, &
                  e_enrg(i_enrg), t_enrgdis_ion) &
                  * e_denrg(i_enrg) / cs_ion(i_enrg,i_state1,i_sp) * t_sc 
                Ael_sec(i,i_enrg) = Ael_sec(i,i_enrg) &
                  + den(i_z,im_e_reactant_all(i_sp)) &
                  * ecs(i_enrg,i_state2,i_sp) &
                  * difcs_sec(e_enrg(i), e_enrg(i_enrg), t_enrgdis_ion) &
                  * e_denrg(i_enrg) / cs_sec(i_enrg,i_state1,i_sp)
              end do
              sum_cs_ion(i_enrg,i_sp) = sum_cs_ion(i_enrg,i_sp) &
                + ecs(i_enrg,i_state2,i_sp) 
            end if
          end do
        end do
      end do

      Ael = Ael_exc + Ael_ion + Ael_sec

      ! --- Calculate Electron Flux -------------------------------------------
      sigmaN = zero
      do concurrent (i_enrg = 1:n_e_enrg)
        do i_sp = 1, n_sp_e_thk
          sigmaN(i_enrg) = sigmaN(i_enrg) + den(i_z,im_e_reactant_all(i_sp)) &
            * (sum_cs_exc(i_enrg,i_sp) + sum_cs_ion(i_enrg,i_sp))
        end do
      end do
      
      ! only 1st and 4th term of Eq. 4 from photoelectron notes for highest
      ! energy bin
      eflux(i_z,n_e_enrg) = Sel(i_z,n_e_enrg) / sigmaN(n_e_enrg)
      
      do i_enrg = n_e_enrg-1, 1, -1
      ! Eq. 20 in photoelectron notes, multiplied by dE
        t_aj = den(i_z,iELE) &
          * loss(e_enrg(i_enrg), den(i_z,iELE), Te(i_z)) &
          + e_denrg(i_enrg) * sigmaN(i_enrg)
        t_aj1 = den(i_z,iELE) &
          * loss(e_enrg(i_enrg+1), den(i_z,iELE), Te(i_z)) &
          - e_denrg(i_enrg+1) * sigmaN(i_enrg+1)
        if (t_aj < 1.0E-20) then    ! if t_aj = 0
          eflux(i_z,i_enrg) = 0.57 * den(i_z,iELE)
        else
          eflux(i_z,i_enrg) = (e_denrg(i_enrg) * Sel(i_z,i_enrg) &
            + e_denrg(i_enrg+1) * Sel(i_z,i_enrg+1) &
            + t_aj1 * eflux(i_z,i_enrg+1) &
            + e_denrg(i_enrg) &
            * dot_product(Ael(i_enrg,i_enrg+1:n_e_enrg), &
              eflux(i_z,i_enrg+1:n_e_enrg)) &
            + e_denrg(i_enrg+1) &
            * dot_product(Ael(i_enrg+1,i_enrg+2:n_e_enrg), &
              eflux(i_z,i_enrg+2:n_e_enrg)) &
            ) / t_aj
        end if
      end do
    end if
  end do 
  
  !----------------------------------------------------------------------------
  ! Secondary Production
  !----------------------------------------------------------------------------

  ! pS = zero
  ! do concurrent (i_z = iz_bot:n_z, i_sp = 1:n_sp_e_thk, i_enrg = 1:n_e_enrg)
    ! pS(i_z) = pS(i_z) + den(i_z,im_e_reactant_all(i_sp)) &
      ! * ecs_ion(i_enrg,i_sp) * eflux(i_z,i_enrg) * e_denrg(i_enrg) 
  ! end do
  
  i_r = 0
  eph = zero
  do i_sp = 1, n_sp_e_thk
    do i_branch = 1, n_branch_e(i_sp)
      i_r = i_r + 1
      do i_z = iz_bot, n_z
        do i_enrg = 1, n_e_enrg
          eph(i_r,i_z) = eph(i_r,i_z) + eflux(i_z,i_enrg) &
          * ecs(i_enrg,n_states(1,i_sp)+i_branch,i_sp) * e_denrg(i_enrg) 
        end do
      end do
    end do
  end do
  
  ! do i_sp = n_sp_e_thk+1, n_sp_e_thk+n_sp_e_thn
    ! do i_branch = 1, n_branch_e(i_sp)
      ! i_r = i_r + 1
      ! do i_z = iz_bot, n_z
        ! do i_enrg = 1, n_e_enrg
          ! eph(i_r,i_z) = eph(i_r,i_z) + branch_ratio_el(i_enrg,i_branch,i_sp) &
          ! * eflux(i_z,i_enrg) * crs_tot_inel(i_enrg,i_sp) * e_denrg(i_enrg) 
        ! end do
      ! end do
    ! end do
  ! end do

  deallocate(Ael_exc, Ael_sec, Ael_ion, Ael)
  deallocate(sum_cs_ion, sum_cs_exc)
  deallocate(sigmaN)

end subroutine cal_electrons


real(kind=wp) pure function difcs_sec(E_sec, E, I)
! This function calculates the shape of the differential cross-section for
! secondary electron production, with secondary electron energy E_sec, primary
! electron energy E and ionization potential I (in eV).
! See Eq. (15) in photoelectron notes.  
  use types, only: wp => dp
  use constants
  implicit none
  real(wp), intent(in) :: E_sec, E, I
  real(wp), parameter :: E0 = 13.0 ! shape parameter for N2
  !? need to change E0 to CO2 for Mars
  difcs_sec = (one + E_sec/E0)**(-2.1_wp) / E0 / tanh((E-I) / two / E0)
end function difcs_sec


real(kind=wp) pure function loss(E, n_e, Te)
! This function calculates the degradation of suprathermal electrons with 
! energy E from collisions with ambient thermal electrons of density n_e and
! temperature Te. 
! See Eqs. (10) and (11) in photoelectron notes.
  use types, only: wp => dp
  use constants
  implicit none
  real(wp), intent(in) :: E, n_e, Te
  real(wp) :: E_e

  E_e = 8.618E-5_wp * Te
  if (E_e > E) then
    loss = zero
  else 
    loss = 3.37E-12_wp / E**0.94_wp / n_e**0.03_wp &
      * ((E - E_e) / (E - 0.53_wp * E_e))**2.36_wp
  end if
end function loss

end module
