module reactions_mod

contains

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
  ! constants for Baulch et al. [2005] broadening factor
  real(wp) :: FC, PFL, FCL, NF, CF, FF, XF, rexp
  
  ! loop variables
  integer :: i_rct, i, i_r, i_z
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
     ! DO i_r = 1, ntab_rct
        ! READ(40,"(A12,3X,A12,3X,A12,3X,A12,3X,A12)") (fn(i),i=1,5)
        ! DO i = 1, 5
           ! itab_rct(i,i_r) = FIND_NAME(fn(i),name)
        ! END DO
        ! READ(40,*) nprs_rct(i_r),ntmp_rct(i_r)
        ! READ(40,*) (tmp_rct(nt,i_r),nt=1,ntmp_rct(i_r))
        ! DO np = 1, nprs_rct(i_r)
           ! READ(40,*) rdum,(rct_tab(np,nt,i_r),nt=1,ntmp_rct(i_r))
           ! plog_rct(np,i_r) = LOG(rdum)
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
  read (fid_nrct,'(A)') ts_line
  read (ts_line,*) n_nrct, ts1
  read (fid_nrct,'(A)') ts_line
  read (fid_nrct,'(A)') ts_line

  ! count how many of the entries are actually active (type /= 0)
  i_r = 0
  do i_rct = 1, n_nrct
    read(fid_nrct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) i_r = i_r + 1
  end do

  rewind(unit=fid_nrct)

  ! determine number of entries in ion reactions list
  read (fid_irct,'(A)') ts_line
  read (ts_line,*) n_irct, ts1
  read (fid_irct,'(A)') ts_line
  read (fid_irct,'(A)') ts_line

  ! count how many of the entries are actually active (type /= 0)
  do i_rct = 1, n_irct
    read(fid_irct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) i_r = i_r + 1
  end do
  n_rct = i_r ! assign to global variable

  rewind(unit=fid_irct)

  ! now allocate global variables and actually populate them
  allocate(ctitle(n_rct), chem_type(n_rct), rkn(10,n_rct), ireactant(5,n_rct))
  ! allocate(ntab(n_rct))

  ! --- Neutral reactions -----------------------------------------------------
  read (fid_nrct,'(A)') ts_line
  read (fid_nrct,'(A)') ts_line
  read (fid_nrct,'(A)') ts_line  
  i_r = 0
  do i_rct = 1, n_nrct
    read(fid_nrct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) then
      i_r = i_r + 1
      chem_type(i_r) = tm_chem_type
      ctitle(i_r) = ts_species(1)//' + '//ts_species(2)//' = '// & 
        ts_species(3)//' + '//ts_species(4)//' + '//ts_species(5)
      ! determine indices for products
      do concurrent (i = 3:5)
        ireactant(i,i_r) = find_name(ts_species(i), sp_list)
      end do
      
      ! determine indices for reactants
      if (abs(chem_type(i_r)) == 1) then
      ! unimolecular reaction
        rkn(:,i_r) = t_rk
        i = 1
        ireactant(i,i_r) = find_name(ts_species(i), sp_list)
        if (ireactant(i,i_r) <= 0) then
          write(*,'(I6, ":", 2X, A)') i_r, ctitle(i_r)
          write(*,'("Error: Reactant ", I2, " not found: ", A12)') & 
            i, ts_species(i)
          stop
        end if
        ireactant(2,i_r) = 0

      else if ((abs(chem_type(i_r)) > 1) .and. (abs(chem_type(i_r)) < 6)) then
      ! bimolecular & trimolecular reaction
        rkn(:,i_r) = t_rk
        do i = 1, 2
          ireactant(i,i_r) = find_name(ts_species(i), sp_list)
          if (ireactant(i,i_r) <= 0) then
            write(*,'(I6, ":", 2X, A)') i_r, ctitle(i_r)
            write(*,'("Error: Reactant ", I2, " not found: ", A12)') & 
              i, ts_species(i)
            stop
          end if
        end do
         
      ! Tabulated reaction rates (deactivated for now)
      ! else if(abs(chem_type(i_r)) == 6) then
        ! do i = 1, 2
          ! ireactant(i,i_r) = find_name(ts_species(i), sp_list)
          ! if (ireactant(i,i_r) <= 0) then
            ! write(*,"(I6,':',2X,A)") i_r, ctitle(i_r)
            ! write(*,"('Error: Reactant ',I2,' not found: ',A12)") & 
              ! i, ts_species(i)
            ! stop
          ! end if
        ! end do
        ! do nt = 1, ntab_rct
          ! ltab = .true.
          ! do i = 1, 5
            ! if(ireactant(i,i_r) /= itab_rct(i,nt)) ltab = .false.
        ! end do
        ! if(ltab) then
          ! ntab(i_r) = nt
          ! EXIT
        ! end if
      ! end do
      ! if(.not.ltab) then
        ! write(*,"(I6,':',2X,A)") i_r,ctitle(i_r)
        ! write(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") i_r
        ! STOP
      ! end if

      else 
        write(*,'("Error: Invalid setting for neutral reaction ", I0, &
          "! Exiting...")') i_rct
        stop
      end if
    end if
  end do
  close(unit=fid_nrct)

  ! --- Ion reactions ---------------------------------------------------------
  read (fid_irct,'(A)') ts_line
  read (fid_irct,'(A)') ts_line
  read (fid_irct,'(A)') ts_line  
  do i_rct = 1, n_irct
    read(fid_irct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) then
      i_r = i_r + 1
      chem_type(i_r) = tm_chem_type
      ctitle(i_r) = ts_species(1)//' + '//ts_species(2)//' = '// & 
        ts_species(3)//' + '//ts_species(4)//' + '//ts_species(5)
      ! determine indices for products
      do concurrent (i = 3:5)
        ireactant(i,i_r) = find_name(ts_species(i), sp_list)
      end do
      
      ! determine indices for reactants
      if (abs(chem_type(i_r)) == 1) then
      ! unimolecular reaction
        rkn(:,i_r) = t_rk
        i = 1
        ireactant(i,i_r) = find_name(ts_species(i), sp_list)
        if (ireactant(i,i_r) <= 0) then
          write(*,'(I6, ":", 2X, A)') i_r, ctitle(i_r)
          write(*,'("Error: Reactant ", I0, " not found: ", A12)') &
            i, ts_species(i)
          stop
        end if
        ireactant(2,i_r) = 0

      else if ((abs(chem_type(i_r)) > 1) .and. (abs(chem_type(i_r)) < 6)) then
      ! bimolecular & trimolecular reaction
      rkn(:,i_r) = t_rk
        do i = 1, 2
          ireactant(i,i_r) = find_name(ts_species(i), sp_list)
          if (ireactant(i,i_r) <= 0) then
            write(*,'(I6, ":", 2X, A)') i_r, ctitle(i_r)
            write(*,'("Error: Reactant ", I0, " not found: ", A12)') & 
              i, ts_species(i)
            stop
          end if
        end do

      ! Tabulated reaction rates (deactivated for now)
      ! else if(ABS(chem_type(i_r)) == 6) then
        ! do i = 1, 2
          ! ireactant(i,i_r) = FIND_NAME(ts_species(i),sp_list)
            ! if (ireactant(i,i_r) <= 0) then
              ! write(*,"(I6,':',2X,A)") i_r,ctitle(i_r)
              ! write(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)") i_r,i,ts_species(i)
              ! STOP
            ! end if
          ! end do
           
          ! do nt = 1, ntab_rct
            ! ltab = .true.
            ! do i = 1, 5
              ! if(ireactant(i,i_r) /= itab_rct(i,nt)) ltab = .false.
            ! end do
            ! if(ltab) then
              ! ntab(i_r) = nt
              ! EXIT
            ! end if
          ! end do
          ! if(.not.ltab) then
            ! write(*,"(I6,':',2X,A)") i_r,ctitle(i_r)
            ! write(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") i_r
            ! STOP
          ! end if

      else 
        write(*,'("Error: Invalid setting for ion reaction ", I0, &
          "! Exiting...")') i_rct
        stop
      end if
    end if
  end do
  close(unit=fid_irct)

  ! deallocate(itab_rct)

  !----------------------------------------------------------------------------
  !  Calculate reaction coefficients
  !----------------------------------------------------------------------------
  
  allocate(rk(n_rct,n_z), rct_rate(n_rct,n_z)) 
  
  do i_rct = 1, n_rct
    
    ! --- Neutral reactions ---------------------------------------------------
    if (chem_type(i_rct) == 1) then
    ! unimolecular reactions
    ! rate coefficient is simply rate constant
      do concurrent (i_z = 1:n_z)
        rk(i_rct,i_z) = rkn(1,i_rct)
      end do

    else if(chem_type(i_rct) == 2) then
    ! bimolecular reactions
    ! k1 * T^k2 * exp(k3 * T)
      do concurrent (i_z = 1:n_z)
        rk(i_rct,i_z) = rkn(1,i_rct) * (Tn(i_z) ** rkn(2,i_rct)) &
          * exp(rkn(3,i_rct) / Tn(i_z))
      end do

    else if (chem_type(i_rct) == 3) then
    ! association reactions
    ! k0 * [M] * kinf / (k0 * [M] + kinf)
    ! k0 and kinf are in the form k1 * T^k2 * exp(k3 * T)
    ! F is the Troe broadening factor
    ! log F = log Fc / (1 + (log(P)+C) / (N - 0.14*(log(P)+C))^2)
    ! where
    ! reduced pressure P = k0 * [M] / kinf
    ! N = 0.75 - 1.27 * log Fc
    ! C = -0.4 - 0.67 * log Fc
      do i_z = 1, n_z
        rkInf = rkn(1,i_rct) * (Tn(i_z) ** rkn(2,i_rct)) &
          * exp(rkn(3,i_rct) / Tn(i_z))
        rk0 = rkn(4,i_rct) * (Tn(i_z) ** rkn(5,i_rct)) &
          * exp(rkn(6,i_rct) / Tn(i_z))
        if(rkn(10,i_rct) == zero) then
          rk(i_rct,i_z) = rk0 * rkInf * den(i_z,0) / (rkInf + rk0 * den(i_z,0))
        else
          FC = rkn(10,i_rct)
          PFL = log10(rk0 * den(i_z,0) / rkInf)
          NF = 0.75_wp - 1.27_wp * log10(FC)
          CF = -0.4_wp - 0.67_wp * log10(FC)
          FF =10._wp**(log10(FC)/(one+((PFL+CF)/(NF-0.14_wp*(PFL+CF)))**2))
          rk(i_rct,i_z) = FF*(rk0*rkInf/(rk0*den(i_z,0)+rkInf))
        end if
      end do

    else if (chem_type(i_rct) == 4) then
    ! association & radiative association reactions
      do i_z = 1, n_z
        rkInf = rkn(1,i_rct) * (Tn(i_z) ** rkn(2,i_rct)) &
          * exp(rkn(3,i_rct) / Tn(i_z))
        rk0 = rkn(4,i_rct) * (Tn(i_z) ** rkn(5,i_rct)) &
          * exp(rkn(6,i_rct) / Tn(i_z))
        rkRad = rkn(7,i_rct) * (Tn(i_z) ** rkn(8,i_rct)) &
          * exp(rkn(9,i_rct) / Tn(i_z))
        if(rkn(10,i_rct) == zero) then
          rk(i_rct,i_z) = rkRad + rk0 * rkInf * den(i_z,0) &
            / (rkInf + rk0 * den(i_z,0))
        else
          FC = rkn(10,i_rct)
          PFL = log10(rk0*den(i_z,0)/rkInf)
          FCL = log10(FC)
          NF = 0.75_wp - 1.27_wp*FCL
          CF = -0.4_wp-0.67_wp*FCL
          FF =10._wp**(FCL/(one+((PFL+CF)/(NF-0.14_wp*(PFL+CF)))**2))
          XF = FF/(one-FF)
          rk(i_rct,i_z) = min(rkInf, rkRad + FF*(rk0*rkInf*den(i_z,0)/(rk0*den(i_z,0)+rkInf)))
        end if
      end do

    else if (chem_type(i_rct) == 5) then
    ! association reactions (Sander's formula)
    ! k0 * [M] * kinf / (k0 * [M] + kinf)
    !   * 0.6 ^ (1 / (1 + (log10(k0 * [M] / kinf))^2)
      do i_z = 1, n_z
        rkInf = rkn(1,i_rct) * (Tn(i_z) ** rkn(2,i_rct)) &
          * exp(rkn(3,i_rct) / Tn(i_z))
        rk0 = rkn(4,i_rct) * (Tn(i_z) ** rkn(5,i_rct)) &
          * exp(rkn(6,i_rct) / Tn(i_z))
        rexp = one / (one + (log10(rk0 * den(i_z,0) / rkInf))**two)
        rk(i_rct,i_z) = rk0 * rkInf * den(i_z,0) / (rkInf + rk0 * den(i_z,0)) &
          * 0.6_wp ** rexp
      end do

    ! Tabulated reaction rates (deactivated for now)
    ! else if(chem_type(i_rct) == 6) then
      ! nt = ntab(i_rct)
      ! do i_z = 1, n_z
        ! if(Tn(i_z) < tmp_rct(1,nt)) then
          ! nt1 = 1
          ! nt2 = 2
          ! rtmp = zero
        ! else if(Tn(i_z) > tmp_rct(ntmp_rct(nt),nt)) then
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
        ! rk(i_rct,i_z) = (one-rtmp)*(one-rprs)*rct_tab(np1,nt1,nt) &
                 ! + rtmp*(one-rprs)*rct_tab(np2,nt1,nt) &
                 ! + (one-rtmp)*rprs*rct_tab(np1,nt2,nt) &
                 ! + rtmp*rprs*rct_tab(np2,nt2,nt) 
        ! rk(i_rct,i_z) = exp(rk(i_rct,i_z))
      ! end do

    ! --- Ion reactions -------------------------------------------------------
    else if (chem_type(i_rct) == -1) then
    ! unimolecular reaction
      do i_z = 1, n_z
        rk(i_rct,i_z) = rkn(1,i_rct) * (Tn(i_z) ** rkn(2,i_rct)) &
          * two * exp(rkn(3,i_rct) / Tn(i_z)) / (one + exp(rkn(4,i_rct) / Tn(i_z)))
      end do

    else if (chem_type(i_rct) == -2) then
    ! normal two-body reaction
      do i_z = 1, n_z
        rk(i_rct,i_z)=rkn(1,i_rct)*(Tn(i_z)**rkn(2,i_rct)) &
          *two*exp(rkn(3,i_rct)/Tn(i_z))/(one+exp(rkn(4,i_rct)/Tn(i_z)))
      end do

    else if (chem_type(i_rct) == -3) then
    ! 3-body reaction
      do i_z = 1, n_z
        rk(i_rct,i_z)=rkn(1,i_rct)*(Tn(i_z)**rkn(2,i_rct)) &
          *two*exp(rkn(3,i_rct)/Tn(i_z))/(one+exp(rkn(4,i_rct)/Tn(i_z))) 
      end do

    else if (chem_type(i_rct) == -4) then
    ! electron recombination
      do i_z = 1, n_z
        rk(i_rct,i_z)=rkn(1,i_rct)*(te(i_z)**rkn(2,i_rct)) &
          *two*exp(rkn(3,i_rct)/te(i_z))/(one+exp(rkn(4,i_rct)/te(i_z)))
      end do
    
    else if (chem_type(i_rct) .ne. 0) then
      write(*,'("Error! Invalid setting for reaction ", I0, ": ", A)') &
        i_rct, ctitle(i_rct)
      write(*,'("Exiting...")')
      stop
    end if

  end do

end subroutine read_reactions

end module
