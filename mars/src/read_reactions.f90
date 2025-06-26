! This subroutine reads in the list of reactions and their rate constants in 
! nreactions.csv and ireactions.csv. 
subroutine read_reactions

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use global_variables
  use subs, only : find_name

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  implicit none
  ! file unit numbers for neutral and ion reactions
  integer :: fid_nrct, fid_irct
  ! number of neutral and ion reactions
  integer :: n_nrct, n_irct
  
  ! loop variables
  integer :: i_rct, i, nr
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
     ! DO nr = 1, ntab_rct
        ! READ(40,"(A12,3X,A12,3X,A12,3X,A12,3X,A12)") (fn(i),i=1,5)
        ! DO i = 1, 5
           ! itab_rct(i,nr) = FIND_NAME(fn(i),name)
        ! END DO
        ! READ(40,*) nprs_rct(nr),ntmp_rct(nr)
        ! READ(40,*) (tmp_rct(nt,nr),nt=1,ntmp_rct(nr))
        ! DO np = 1, nprs_rct(nr)
           ! READ(40,*) rdum,(rct_tab(np,nt,nr),nt=1,ntmp_rct(nr))
           ! plog_rct(np,nr) = LOG(rdum)
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
  nr = 0
  do i_rct = 1, n_nrct
    read(fid_nrct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) nr = nr + 1
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
    if (tm_chem_type /= 0) nr = nr + 1
  end do
  n_rct = nr ! assign to global variable

  rewind(unit=fid_irct)

  ! now allocate global variables and actually populate them
  allocate(ctitle(n_rct), chem_type(n_rct), rk(10,n_rct), irct(5,n_rct))
  ! allocate(ntab(n_rct))

  ! neutral reactions
  read (fid_nrct,'(A)') ts_line
  read (fid_nrct,'(A)') ts_line
  read (fid_nrct,'(A)') ts_line  
  nr = 0
  do i_rct = 1, n_nrct
    read(fid_nrct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) then
      nr = nr + 1
      chem_type(nr) = tm_chem_type
      ctitle(nr) = ts_species(1)//' + '//ts_species(2)//' = '// & 
        ts_species(3)//' + '//ts_species(4)//' + '//ts_species(5)
      ! determine indices for products
      do concurrent (i = 3:5)
        irct(i,nr) = find_name(ts_species(i), sp_list)
      end do
      
      ! determine indices for reactants
      ! unimolecular reaction
      if (abs(chem_type(nr)) == 1) then
        rk(:,nr) = t_rk
        i = 1
        irct(i,nr) = find_name(ts_species(i), sp_list)
        if (irct(i,nr) <= 0) then
          write(*,"(I6,':',2X,A)") nr, ctitle(nr)
          write(*,"('Error: Reactant ',I2,' not found: ',A12)") & 
            i, ts_species(i)
          stop
        end if

      ! bimolecular & trimolecular reaction
      else if ((abs(chem_type(nr)) > 1) .and. (abs(chem_type(nr)) < 6)) then
        rk(:,nr) = t_rk
        do i = 1, 2
          irct(i,nr) = find_name(ts_species(i), sp_list)
          if (irct(i,nr) <= 0) then
            write(*,"(I6,':',2X,A)") nr, ctitle(nr)
            write(*,"('Error: Reactant ',I2,' not found: ',A12)") & 
              i, ts_species(i)
            stop
          end if
        end do
         
      ! Tabulated reaction rates (deactivated for now)
      ! else if(abs(chem_type(nr)) == 6) then
        ! do i = 1, 2
          ! irct(i,nr) = find_name(ts_species(i), sp_list)
          ! if (irct(i,nr) <= 0) then
            ! write(*,"(I6,':',2X,A)") nr, ctitle(nr)
            ! write(*,"('Error: Reactant ',I2,' not found: ',A12)") & 
              ! i, ts_species(i)
            ! stop
          ! end if
        ! end do
        ! do nt = 1, ntab_rct
          ! ltab = .true.
          ! do i = 1, 5
            ! if(irct(i,nr) /= itab_rct(i,nt)) ltab = .false.
        ! end do
        ! if(ltab) then
          ! ntab(nr) = nt
          ! EXIT
        ! end if
      ! end do
      ! if(.not.ltab) then
        ! write(*,"(I6,':',2X,A)") nr,ctitle(nr)
        ! write(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") nr
        ! STOP
      ! end if

      else 
        write(*,"('Error: Invalid setting for neutral reaction ', I0, &
          '! Exiting...')") i_rct
        stop
      end if
    end if
  end do
  close(unit=fid_nrct)

  ! ion reactions
  read (fid_irct,'(A)') ts_line
  read (fid_irct,'(A)') ts_line
  read (fid_irct,'(A)') ts_line  
  do i_rct = 1, n_irct
    read(fid_irct,'(A)') ts_line
    read(ts_line,*) tn1, (ts_species(i),i=1,5), tm_chem_type, (t_rk(i),i=1,10)
    if (tm_chem_type /= 0) then
      nr = nr + 1
      chem_type(nr) = tm_chem_type
      ctitle(nr) = ts_species(1)//' + '//ts_species(2)//' = '// & 
        ts_species(3)//' + '//ts_species(4)//' + '//ts_species(5)
      ! determine indices for products
      do concurrent (i = 3:5)
        irct(i,nr) = find_name(ts_species(i), sp_list)
      end do
      
      ! determine indices for reactants
      ! unimolecular reaction
      if (abs(chem_type(nr)) == 1) then
        rk(:,nr) = t_rk
        i = 1
        irct(i,nr) = find_name(ts_species(i), sp_list)
        if (irct(i,nr) <= 0) then
          write(*,"(I6,':',2X,A)") nr, ctitle(nr)
          write(*,"('Error: Reactant ',I0,' not found: ', A12)") &
            i, ts_species(i)
          stop
        end if

      ! bimolecular & trimolecular reaction
      else if ((abs(chem_type(nr)) > 1) .and. (abs(chem_type(nr)) < 6)) then
        rk(:,nr) = t_rk
        do i = 1, 2
          irct(i,nr) = find_name(ts_species(i), sp_list)
          if (irct(i,nr) <= 0) then
            write(*,"(I6,':',2X,A)") nr, ctitle(nr)
            write(*,"('Error: Reactant ',I0,' not found: ', A12)") & 
              i, ts_species(i)
            stop
          end if
        end do

      ! Tabulated reaction rates (deactivated for now)
      ! else if(ABS(chem_type(nr)) == 6) then
        ! do i = 1, 2
          ! irct(i,nr) = FIND_NAME(ts_species(i),sp_list)
            ! if (irct(i,nr) <= 0) then
              ! write(*,"(I6,':',2X,A)") nr,ctitle(nr)
              ! write(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)") nr,i,ts_species(i)
              ! STOP
            ! end if
          ! end do
           
          ! do nt = 1, ntab_rct
            ! ltab = .true.
            ! do i = 1, 5
              ! if(irct(i,nr) /= itab_rct(i,nt)) ltab = .false.
            ! end do
            ! if(ltab) then
              ! ntab(nr) = nt
              ! EXIT
            ! end if
          ! end do
          ! if(.not.ltab) then
            ! write(*,"(I6,':',2X,A)") nr,ctitle(nr)
            ! write(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") nr
            ! STOP
          ! end if

      else 
        write(*,"('Error: Invalid setting for ion reaction ', I0, &
          '! Exiting...')") i_rct
        stop
      end if
    end if
  end do
  close(unit=fid_irct)

  ! deallocate(itab_rct)

end subroutine read_reactions
