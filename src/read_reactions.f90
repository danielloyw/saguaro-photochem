SUBROUTINE READ_REACTIONS

  USE PRECISION
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_NAME
  IMPLICIT NONE

  !  .. Local Variables

  CHARACTER(len=256) :: cline, sdum, header
  CHARACTER(len=12), DIMENSION(5) :: fn 
  LOGICAL :: ltab
  INTEGER :: ntab_rct, nrctmax_neut, nrctmax_ions
  INTEGER :: nr, nt, np, nx, i, ndum, idum
  REAL(RP), DIMENSION(10) :: rckd
  REAL(RP) :: rdum
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itab_rct
  INTEGER, PARAMETER :: ntmp_tab_max = 5
  INTEGER, PARAMETER :: nprs_tab_max = 30

  INTEGER :: ios

  !
  !  .. Read Tabulated Reaction Rate Coefficients
  !

  OPEN(unit=40,file='nreactions.tab',status='old')
     READ(40,*) ntab_rct
     ALLOCATE(itab_rct(5,ntab_rct),nprs_rct(ntab_rct),ntmp_rct(ntab_rct),                           &
          plog_rct(nprs_tab_max,ntab_rct),tmp_rct(ntmp_tab_max,ntab_rct),                           &
          rct_tab(nprs_tab_max,ntmp_tab_max,ntab_rct))
     DO nr = 1, ntab_rct
        READ(40,"(A12,3X,A12,3X,A12,3X,A12,3X,A12)") (fn(i),i=1,5)
        DO i = 1, 5
           itab_rct(i,nr) = FIND_NAME(fn(i),name)
        END DO
        READ(40,*) nprs_rct(nr),ntmp_rct(nr)
        READ(40,*) (tmp_rct(nt,nr),nt=1,ntmp_rct(nr))
        DO np = 1, nprs_rct(nr)
           READ(40,*) rdum,(rct_tab(np,nt,nr),nt=1,ntmp_rct(nr))
           plog_rct(np,nr) = LOG(rdum)
        END DO
     END DO
  CLOSE(unit=40)
  rct_tab = LOG(rct_tab)

  !
  !  .. Reactions from formulae
  !
   
  OPEN(unit=63,file='nreactions.csv',status='old',action='read')
  OPEN(unit=64,file='ireactions.csv',status='old',action='read')

  READ (63,'(A)') cline
  READ (cline,*) nrctmax_neut, sdum
  READ (63,'(A)') header
  READ (63,'(A)') header

  nr = 0
  DO nx = 1, nrctmax_neut
     READ(63,'(A)') cline
     READ(cline,*) ndum, (fn(i),i=1,5), idum, (rckd(i),i=1,10)
     IF (idum /= 0) nr = nr + 1
  END DO

  REWIND(63)

  READ (64,'(A)') cline
  READ (cline,*) nrctmax_ions, sdum
  READ (64,'(A)') header
  READ (64,'(A)') header

  DO nx = 1, nrctmax_ions
     READ(64,'(A)') cline
     READ(cline,*) ndum, (fn(i),i=1,5), idum, (rckd(i),i=1,10)
     IF (idum /= 0) nr = nr + 1
  END DO
  nrct = nr   

  REWIND(64)

  ALLOCATE(itype(nrct),rck(10,nrct),irct(5,nrct),ctitle(nrct),ntab(nrct))

  READ (63,'(A)') header
  READ (63,'(A)') header
  READ (63,'(A)') header  
  nr  = 0
  DO nx = 1, nrctmax_neut
     READ(63,'(A)') cline
     READ(cline,*) ndum, (fn(i),i=1,5), idum, (rckd(i),i=1,10)
     IF ( idum /= 0 ) THEN

        nr = nr + 1
        itype(nr) = idum
        ctitle(nr) = fn(1)//' + '//fn(2)//' = '//fn(3)//' + '//fn(4)//' + '//fn(5)
        DO i = 3, 5
           irct(i,nr) = FIND_NAME(fn(i),name)
        END DO

        IF (ABS(itype(nr)) == 1) THEN                                             ! Unimolecular Reaction
           rck(1:10,nr) = rckd(1:10)
           i = 1
           irct(i,nr) = FIND_NAME(fn(i),name)
           IF (irct(i,nr) <= 0) THEN
              WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
              WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
              STOP
           END IF

        ELSE IF ((ABS(itype(nr)) > 1) .and. (ABS(itype(nr)) < 6)) THEN            ! Bimolecular & Trimolecular Reactions
           rck(1:10,nr) = rckd(1:10)
           DO i = 1, 2
              irct(i,nr) = FIND_NAME(fn(i),name)
              IF (irct(i,nr) <= 0) THEN
                 WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
                 WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
                 STOP
              END IF
           END DO
           
        ELSE IF(ABS(itype(nr)) == 6) THEN                                         ! Tabulated reaction rates
           DO i = 1, 2
              irct(i,nr) = FIND_NAME(fn(i),name)
              IF (irct(i,nr) <= 0) THEN
                 WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
                 WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
                 STOP
              END IF
           END DO
           DO nt = 1, ntab_rct
              ltab = .true.
              DO i = 1, 5
                 IF(irct(i,nr) /= itab_rct(i,nt)) ltab = .false.
              END DO
              IF(ltab) THEN
                 ntab(nr) = nt
                 EXIT
              END IF
           END DO
           IF(.not.ltab) THEN
              WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
              WRITE(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") nr
              STOP
           END IF

        ELSE IF(ABS(itype(nr)) == 7) THEN                                         ! Heterogenous reaction 
           rck(1:10,nr) = rckd(1:10)
           i = 1
           irct(i,nr) = FIND_NAME(fn(i),name)
           IF (irct(i,nr) <= 0) THEN
              WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
              WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
              STOP
           END IF
        END IF
     END IF
  END DO
  CLOSE(unit=63)

  READ (64,'(A)') header
  READ (64,'(A)') header
  READ (64,'(A)') header  
  DO nx = 1, nrctmax_ions
     READ(64,'(A)') cline
     READ(cline,*) ndum, (fn(i),i=1,5), idum, (rckd(i),i=1,10)
     IF ( idum /= 0 ) THEN

        nr = nr + 1
        itype(nr) = idum
        ctitle(nr) = fn(1)//' + '//fn(2)//' = '//fn(3)//' + '//fn(4)//' + '//fn(5)
        DO i = 3, 5
           irct(i,nr) = FIND_NAME(fn(i),name)
        END DO

        IF (ABS(itype(nr)) == 1) THEN                                             ! Unimolecular Reaction
           rck(1:10,nr) = rckd(1:10)
           i = 1
           irct(i,nr) = FIND_NAME(fn(i),name)
           IF (irct(i,nr) <= 0) THEN
              WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
              WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
              STOP
           END IF

        ELSE IF ((ABS(itype(nr)) > 1) .and. (ABS(itype(nr)) < 6)) THEN            ! Bimolecular & Trimolecular Reactions
           rck(1:10,nr) = rckd(1:10)
           DO i = 1, 2
              irct(i,nr) = FIND_NAME(fn(i),name)
              IF (irct(i,nr) <= 0) THEN
                 WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
                 WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
                 STOP
              END IF
           END DO

        ELSE IF(ABS(itype(nr)) == 6) THEN                                         ! Tabulated reaction rates
           DO i = 1, 2
              irct(i,nr) = FIND_NAME(fn(i),name)
              IF (irct(i,nr) <= 0) THEN
                 WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
                 WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: ',A12)"),nr,i,fn(i)
                 STOP
              END IF
           END DO
           
           DO nt = 1, ntab_rct
              ltab = .true.
              DO i = 1, 5
                 IF(irct(i,nr) /= itab_rct(i,nt)) ltab = .false.
              END DO
              IF(ltab) THEN
                 ntab(nr) = nt
                 EXIT
              END IF
           END DO
           IF(.not.ltab) THEN
              WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
              WRITE(*,"(' TABULATED DATA NOT FOUND FOR REACTION ',I4,', ABORTING ....')") nr
              STOP
           END IF

        ELSE IF(ABS(itype(nr)) == 7) THEN                                         ! Heterogenous reaction 
           rck(1:10,nr) = rckd(1:10)
           i = 1
           irct(i,nr) = FIND_NAME(fn(i),name)
           IF (irct(i,nr) <= 0) THEN
              WRITE(*,"(I6,':',2X,A)") nr,ctitle(nr)
              WRITE(*,"(' REACTION',I6,' REACTANT',I2,' NOT FOUND: , ABORTING ... ',A12)"),nr,i,fn(i)
              STOP
           END IF
        END IF
     END IF
  END DO
  CLOSE(unit=64)


963 FORMAT(I4,1X,A12,'+',A12,'-->',A12,'+',A12,'+',A12,'+',I1,3(',',ES9.2,',',F7.3,',',F8.2),',',F6.3,3(',',A24))
940 FORMAT(A12,3X,A12,3X,A12,3X,A12,3X,A12)

  DEALLOCATE(itab_rct)

  RETURN
END SUBROUTINE READ_REACTIONS
