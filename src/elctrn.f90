SUBROUTINE ELCTRN 

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_NAME, FIND_BIN, DIFCS_SEC

  !
  !  .. Declarations
  !

  IMPLICIT NONE

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nprdcts_el
  CHARACTER(len=87), ALLOCATABLE, DIMENSION(:,:) :: elrct_title
  INTEGER, ALLOCATABLE, DIMENSION(:) :: loab_el
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lobr_el
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lopr_el
  INTEGER :: nbrn_max, nprmax_el, nabs_excite, nabs_dissoc, nabs_ioniz, nabs_el
  CHARACTER(len=12) :: yname
  CHARACTER(len=12) :: xname
  CHARACTER(len=256) :: header, cline
  CHARACTER(len=12), DIMENSION(6) :: fm
  INTEGER :: na, n,nm, nb, nx, np, ndum, i, ik, ni, ns, ne, j, npe, iflag, nstart, nend
  REAL(RP) :: yenrg, rdum, E1, E2, E4, Eion, sumf_ion, sumf_sec, enge
  CHARACTER(len=1) :: iret

  iflag = 0

!  IF (iflag == 1) THEN

  !
  !  .. Main Section
  !

  OPEN(unit=60, file='../data/electrons/eimpact.dat', status='old', action='read')

  READ (60,*) nelb, nabs_el_thk, nabs_el_thn, nbrn_max, nprmax_el ! number of bins, number of species (thick), number of species (thin) , ,

  nabs_el = nabs_el_thk + nabs_el_thn ! number of species

  ALLOCATE(elctreV(nelb),elctDeV(nelb),eCS_tot_elast(nelb,nabs_el_thk),                             &
       eCS_tot_inelast(nelb,nabs_el_thk),eCS(nelb,nabs_el_thk,nbrn_max),                            &
       state(nabs_el_thk,nbrn_max),enrgE(nabs_el_thk,nbrn_max),ipath(nabs_el_thk,3),                &
       eCS_exc(nelb,nabs_el_thk),eCS_ion(nelb,nabs_el_thk),                          &
       brat_el(nelb,nbrn_max,nabs_el),                                    &
       nbrnch_el(nabs_el),crs_tot_inel(nelb,nabs_el))

  ALLOCATE(nprdcts_el(nbrn_max,nabs_el),loab_el(nabs_el),lobr_el(nprmax_el,nabs_el),     &
       lopr_el(nprmax_el,nbrn_max,nabs_el),elrct_title(nbrn_max,nabs_el))

  READ(60,"(A)") header
  READ(60,*) (elctreV(ne), ne=1,nelb)
  READ(60,"(A)") header
  READ(60,*) (elctDeV(ne), ne=1,nelb)

  nert = 0; na = 0
  DO n = 1, nabs_el_thk     
     READ(60,"(A12,3I4)") xname, nabs_excite, nabs_dissoc, nabs_ioniz ! number of excited, dissociation and ionization states
     nm = FIND_NAME(xname,name)
     IF(nm < 0) THEN
        WRITE(*,"(' ERROR ')")
        STOP
     END IF
     na = na + 1
     loab_el(na) = nm 
     ipath(na,1) = nabs_excite
     ipath(na,2) = nabs_dissoc
     ipath(na,3) = nabs_ioniz
     nbrnch_el(na) = nabs_dissoc + nabs_ioniz
     READ(60,"(A)") header
     READ(60,*) (eCS_tot_elast(i,na),i=1,nelb)
     READ(60,"(A)") header
     READ(60,*) (eCS_tot_inelast(i,na),i=1,nelb)
     ! .. read excited states
     DO nx = 1, nabs_excite
        READ (60,"(A12,F9.3)") yname, yenrg
        state(na,nx) = yname
        enrgE(na,nx) = yenrg
        READ(60,*) (eCS(ne,na,nx),ne=1,nelb)
     END DO
     ! .. read dissociation & ionization states 
     DO nb = 1, nbrnch_el(na)
        READ (60,"(A12,F9.3)") yname, yenrg
        state(na,nabs_excite+nb) = yname
        enrgE(na,nabs_excite+nb) = yenrg
        READ (60,"(A)") cline
        ! .. read products
        elrct_title(nb,na) = cline ! equation
        READ (cline,"(5(A12,3X),A12)") (fm(j),j=1,6)
        np = 0
        DO j = 3, 6
           nm = FIND_NAME(fm(j),name)
           IF ((nm > 0) .and. (nm <= nsp)) THEN
              np = np + 1
              lopr_el(np,nb,na) = nm
              lobr_el(np,na) = nb
           END IF
        END DO
        nprdcts_el(nb,na) = np
        READ(60,*) (eCS(ne,na,nabs_excite+nb),ne=1,nelb)
     END DO
  END DO
  
  DO n = 1, nabs_el_thn                     
     READ(60,'(A12, I4)') xname, ndum
     nm = FIND_NAME(xname,name)
     IF(nm <= 0) THEN
        WRITE(*,"(' ERROR')")
        STOP
     END IF
     na = na + 1
     loab_el(na) = nm
     nbrnch_el(na) = ndum
     READ(60,'(A)') header
     READ(60,'(10ES11.3)') (crs_tot_inel(ne,na),ne=1,nelb)
     DO nb = 1, nbrnch_el(na)
        READ (60,'(A)') cline
        elrct_title(nb,na) = cline
        READ (cline,"((5(A12,3X),A12,F6.1))") (fm(j),j=1,6), rdum
        np = 0
        DO j = 3, 6
           nm = FIND_NAME(fm(j),name)
           IF ((nm > 0) .and. (nm <= nsp)) THEN
              np = np + 1
              lopr_el(np,nb,na) = nm
              lobr_el(np,na) = nb
           END IF
        END DO
        nprdcts_el(nb,na) = np
        READ(60,'(10ES11.3)') (brat_el(ne,nb,na),ne=1,nelb)
     END DO
  END DO

  CLOSE(60)
  
  ALLOCATE(sum_cs_rees_ion(nelb,nabs_el,100),sum_cs_rees_sec(nelb,nabs_el,100))
  
  !  .. Total cross sections for e-impact excitation and ionization of CO2 and N2

  DO na = 1, nabs_el_thk
     DO ne = 1, nelb
!        eCS_exc(ne,na) = SUM(eCS(ne,na,1:(ipath(1,1)+ipath(1,2))))       
!        eCS_ion(ne,na) = SUM(eCS(ne,na,1:ipath(1,3)))
        nstart = 1
        nend = ipath(na,1) + ipath(na,2)
        eCS_exc(ne,na) = SUM(eCS(ne,na,nstart:nend))       
        nstart = nend + 1
        nend = nend + ipath(na,3)
        eCS_ion(ne,na) = SUM(eCS(ne,na,nstart:nend))
     END DO
  END DO
  
  !  .. Calculation of normalization factor for differential cross section
  
  DO na = 1, nabs_el_thk
     DO ik = 1, ipath(na,3) 
        DO i = 1, nelb
           enge = enrgE(na,ipath(na,1)+ipath(na,2)+ik)
           if(elctreV(i) > enge) then
              E1 = elctreV(i)-enge  
              ni = FIND_BIN(elctreV,elctDeV,E1)
              if(elctreV(ni) > E1) then ! how is remaining energy higher?
                 E4 = elctreV(ni)-half*elctDeV(ni)
                 ni = ni-1
              else if(elctreV(ni) <= E1) then
                 E4 = elctreV(ni)+half*elctDeV(ni)
              endif
              Eion = elctreV(i) - E4
              nm = ni 
              ns = FIND_BIN(elctreV,elctDeV,enge)   !why is this here? Does nothing
              sumf_ion = zero
              sumf_sec = zero
              if(nm == 0) then
                 sum_cs_rees_ion(i,na,ik) = one 
                 sum_cs_rees_sec(i,na,ik) = one
              else
                 DO j = 1, nm
                    E2 = elctreV(i)-elctreV(j) 
                    sumf_ion = sumf_ion + DIFCS_SEC(E2-Eion,elctreV(i),Eion)*elctDeV(j)
                    sumf_sec = sumf_sec + DIFCS_SEC(elctreV(j),elctreV(i),Eion)*elctDeV(j)
                 ENDDO
                 sum_cs_rees_ion(i,na,ik) = sumf_ion
                 sum_cs_rees_sec(i,na,ik) = sumf_sec
                 IF(ABS(sumf_sec) < 1.E-30_RP) THEN
                    WRITE(*,"(' sumf_sec = ',ES11.3,' nm = ',I6)") sumf_sec,nm
                    DO j = 1, nm
                       WRITE(*,"(I6,2ES11.3)") j, elctDeV(j),  DIFCS_SEC(elctreV(j),elctreV(i),Eion)
                    END DO
                    READ(*,"(A)") iret
                 END IF
              endif
           endif
        ENDDO
     ENDDO
  ENDDO

  !
  !  .. Move things to more convenient arrays
  !

  nert = 0
  DO na = 1, nabs_el
     DO nb = 1, nbrnch_el(na)
        nert = nert + 1
    END DO
  END DO

  ALLOCATE(etitle(nert),iert(5,nert))
  npe = 0; iert(:,:) = 0
  DO na = 1, nabs_el
     DO nb = 1, nbrnch_el(na)
        npe = npe + 1
        etitle(npe) = elrct_title(nb,na)
        iert(1,npe) = loab_el(na)
        DO np = 1, nprdcts_el(nb,na)
           iert(np+1,npe) = lopr_el(np,nb,na)
        END DO
    END DO
  END DO

  !  .. Deallocate local arrays

  DEALLOCATE(nprdcts_el,loab_el,lobr_el,lopr_el,elrct_title)

  !  .. Allocate some external arrays

!  END IF

  ALLOCATE(esrc(nphrt,nlev,nelb),Sel(nlev,nelb),eflux(nlev,nelb),pS(nlev),eph(nert,nlev),rpe(nert,nlev))
  esrc = zero
  Sel = zero
  eflux(:,:) = zero
  pS = zero
  eph = zero
  rpe = zero

  RETURN
END SUBROUTINE ELCTRN
