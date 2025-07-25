module elctrn_mod

contains

subroutine elctrn
! This subroutine sets up the calculation for electron impact reactions. 

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
  ! file unit numbers for eimpact.dat
  integer :: fid
  ! maximum states for each species with cross sections
  integer :: n_state_max
  ! index mapping from list of electron impact species -> list of all species
  integer, allocatable, dimension(:) :: im_e_reactant_all
  ! index mapping from list of products in photo reaction -> list of all
  ! species
  integer, allocatable, dimension(:,:,:) :: im_e_product_all
  ! number of reactions (excitation, dissociation, ionization)
  integer :: n_excite, n_dissoc, n_ioniz
  ! number of dissociation + ionization reactions
  integer :: n_sp_e
  ! discretized ionization energy
  real(wp) :: enrg_ion

  ! loop variables
  integer :: i, i_sp, i_branch, i_state, i_enrg, i_p, i_r
  ! temporary / dummy variables
  character(len=256) :: ts_line
  character(len=10) :: ts_state
  character(len=87), allocatable, dimension(:,:) :: ts_title
  character(len=12) :: ts_species1
  character(len=12), dimension(6) :: ts_species2
  real(wp), allocatable, dimension(:) :: t_cs
  real(wp) :: t_enrg1, t_enrg2, t_enrg3
  integer :: tn1, tn2, tn_isp, tn_ibin

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
    read(fid,'(A)') ts_line ! "TOTAL ELASTIC CROSS SECTION"
    read(fid,*) (t_cs(i_enrg), i_enrg=1,n_e_enrg)
    read(fid,'(A)') ts_line ! "TOTAL INELASTIC CROSS SECTION"
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
      do i = 3, 6 ! products
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
  allocate(ecs_exc(n_e_enrg,n_sp_e_thk), ecs_ion(n_e_enrg,n_sp_e_thk))
  do i_sp = 1, n_sp_e_thk
    do i_enrg = 1, n_e_enrg
      tn2 = n_states(1,i_sp) + n_states(2,i_sp)
      ecs_exc(i_enrg,i_sp) = sum(ecs(i_enrg,1:tn2,i_sp))
      tn1 = tn2 + 1
      tn2 = tn2 + n_states(3,i_sp)
      ecs_ion(i_enrg,i_sp) = sum(ecs(i_enrg,tn1:tn2,i_sp))
    end do
  end do
  
  ! calculate cross sections by summing differential cross sections
  allocate(cs_ion(n_e_enrg,n_state_max,n_sp_e))
  allocate(cs_sec(n_e_enrg,n_state_max,n_sp_e))
  do i_sp = 1, n_sp_e_thk
    do i_state = 1, n_states(3,i_sp)
      do i_enrg = 1, n_e_enrg
        t_enrg1 = enrg_state(n_states(1,i_sp)+n_states(2,i_sp)+i_state,i_sp)
        if (e_enrg(i_enrg) > t_enrg1) then
        ! electron energy > ionization energy
          
          ! --- Calculate discretized ionization energy -----------------------
          ! In modeling the electron cascade, electrons within the initial bin
          ! are pushed up against the higher-energy side of the bin. When a
          ! reaction results in the electron falling into the lower half of the
          ! final bin, it is likely that the final bin is actually one lower
          ! without the "pushing up" in the initial bin. Electrons falling into
          ! the upper half of the final bin will get pushed up. 
          t_enrg2 = e_enrg(i_enrg) - t_enrg1   ! total energy after ionization
          tn_ibin = find_bin(t_enrg2, e_enrg, e_denrg)
          if (t_enrg2 < e_enrg(tn_ibin)) then
          ! total energy after ionization < bin center
            ! discretized total energy after ionization
            t_enrg3 = e_enrg(tn_ibin) - half * e_denrg(tn_ibin)
            tn_ibin = tn_ibin-1
          else   ! total energy after ionization >= bin center
            ! then push electron energy
            t_enrg3 = e_enrg(tn_ibin) + half * e_denrg(tn_ibin)
          end if
          enrg_ion = e_enrg(i_enrg) - t_enrg3
          ! -------------------------------------------------------------------
          
          if (tn_ibin == 0) then
          ! when exiting e energy < first bin center
            ! then remove electron from model with large 1-cm2 cross section
            cs_ion(i_enrg,i_state,i_sp) = one 
            cs_sec(i_enrg,i_state,i_sp) = one
          else
            do i = 1, tn_ibin
              ! exiting electron cross-section
              cs_ion(i_enrg,i_state,i_sp) = cs_ion(i_enrg,i_state,i_sp) &
                + e_denrg(i) * difcs_sec(e_enrg(i_enrg)-e_enrg(i)-enrg_ion, &
                e_enrg(i_enrg), enrg_ion)
              ! secondary electron cross-section
              cs_sec(i_enrg,i_state,i_sp) = cs_sec(i_enrg,i_state,i_sp) &
                + e_denrg(i) * difcs_sec(e_enrg(i), e_enrg(i_enrg), enrg_ion)
            end do
          end if
        end if
      end do
    end do
  end do

  ! allocate external variables

  allocate(esrc(n_prct,n_z,n_e_enrg))
  allocate(Sel(n_z,n_e_enrg))
  allocate(eflux(n_z,n_e_enrg))
  allocate(pS(n_z))
  allocate(eph(n_erct,n_z))
  allocate(rpe(n_erct,n_z))
  esrc = zero
  Sel = zero
  eflux = zero
  pS = zero
  eph = zero
  rpe = zero

  ! deallocate local arrays
  deallocate(im_e_reactant_all, im_e_product_all, ts_title, t_cs)
  
end subroutine elctrn


integer pure function find_bin(enrg, e_enrg, e_denrg)
! This function finds the index of the bin containing the value enrg in an
! monotically increasing vector e_enrg with bin sizes e_denrg.
! Returns extreme bins if outside of the range spanned by e_enrg. 
  use types, only: wp => dp
  use constants
  implicit none
  real(wp), intent(in) :: enrg
  real(wp), intent(in), dimension(:) :: e_enrg, e_denrg
  integer :: i, n_e_enrg
  
  n_e_enrg = size(e_enrg)
  
  ! initialize result to first bin
  find_bin = 1
  
  do i = 2, n_e_enrg-1
    if ((e_enrg(i) - half * e_denrg(i) <= enrg) .and. &
      (enrg <= e_enrg(i+1) - half * e_denrg(i+1))) then
      find_bin = i
      return
     end if
  end do
  
  ! case for last bin
  if (enrg > e_enrg(n_e_enrg) - half * e_denrg(n_e_enrg)) find_bin = n_e_enrg
  
end function find_bin


real(kind=wp) pure function difcs_sec(E_sec, E, I)
! This function calculates the differential cross-section for secondary
! electron production, with secondary electron energy E_sec, primary electron
! energy E and ionization potential I (in eV).
! See Eq. (15) in photoelectron notes.  
  use types, only: wp => dp
  use constants
  implicit none
  real(wp), intent(in) :: E_sec, E, I
  real(wp), parameter :: E0 = 13.0 ! shape parameter for N2
  !? need to change E0 to CO2 for Mars
  difcs_sec = (one + E_sec/E0)**(-2.1_wp) / E0 / tanh((E-I) / two / E0)
end function difcs_sec

end module
