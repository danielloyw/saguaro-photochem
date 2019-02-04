SUBROUTINE COMPOUT

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : INDEXX

  !
  !  .. Declarations
  !

  IMPLICIT none
  
  REAL(RP) :: sm,sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,rfac
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: pcolrate
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ecolrate
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: colrate
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indxp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indxr
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: cprd_ext, cprd_ph, clss_ph, cprd_pe, clss_pe, cprd_chem,   &
       clss_chem, cprd_net, clss_net, crc, ctp, cbt, cbl, csf, cdens, tbal
  REAL(RP) :: cprd_ext_oxy,cprd_ph_oxy,clss_ph_oxy,cprd_pe_oxy,clss_pe_oxy,cprd_chem_oxy,      &
       clss_chem_oxy,cprd_net_oxy,clss_net_oxy,crc_oxy,ctp_oxy,cbt_oxy,cbl_oxy
  REAL(RP) :: cprd_ext_hyd,cprd_ph_hyd,clss_ph_hyd,cprd_pe_hyd,clss_pe_hyd,cprd_chem_hyd,      &
       clss_chem_hyd,cprd_net_hyd,clss_net_hyd,crc_hyd,ctp_hyd,cbt_hyd,cbl_hyd
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indx
  INTEGER :: nz, nm, np, nx, nr, nl, ne, nw, nw_high, na, nb

  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wcrs
  
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: solar_flux, COphotorate
  
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: den_oxy, mol_oxy, flx_oxy, prext_oxy, pr_ph_oxy, ls_ph_oxy, &
       pr_pe_oxy, ls_pe_oxy, pr_chem_oxy, ls_chem_oxy, pr_oxy, ls_oxy, rcdn_oxy, div_flx_oxy

  REAL(RP), ALLOCATABLE, DIMENSION(:) :: den_hyd, mol_hyd, flx_hyd, prext_hyd, pr_ph_hyd, ls_ph_hyd, &
       pr_pe_hyd, ls_pe_hyd, pr_chem_hyd, ls_chem_hyd, pr_hyd, ls_hyd, rcdn_hyd, div_flx_hyd

  CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: heads, dheads

  ALLOCATE(pcolrate(nphrt,3),ecolrate(nert),colrate(nrct,3),indxr(nrct),indxp(nphrt), cprd_ext(nsp),    &
       cprd_ph(nsp), clss_ph(nsp), cprd_pe(nsp), clss_pe(nsp), cprd_chem(nsp), clss_chem(nsp),      &
       cprd_net(nsp), clss_net(nsp), crc(nsp), ctp(nsp), cbt(nsp), cbl(nsp), csf(nsp), indx(nsp),   &
       cdens(nsp), tbal(nsp))

  ALLOCATE(den_oxy(nlev), mol_oxy(nlev), flx_oxy(nlev), prext_oxy(nlev), pr_ph_oxy(nlev),           &
       ls_ph_oxy(nlev), pr_pe_oxy(nlev), ls_pe_oxy(nlev), pr_chem_oxy(nlev), ls_chem_oxy(nlev),     &
       pr_oxy(nlev), ls_oxy(nlev), rcdn_oxy(nlev), div_flx_oxy(nlev))

  ALLOCATE(den_hyd(nlev), mol_hyd(nlev), flx_hyd(nlev), prext_hyd(nlev), pr_ph_hyd(nlev),           &
       ls_ph_hyd(nlev), pr_pe_hyd(nlev), ls_pe_hyd(nlev), pr_chem_hyd(nlev), ls_chem_hyd(nlev),     &
       pr_hyd(nlev), ls_hyd(nlev), rcdn_hyd(nlev), div_flx_hyd(nlev))

  ALLOCATE(wcrs(ncrsA+ncrsB_low+ncrsC), solar_flux(ncrsA+ncrsB_low+ncrsC-2, nlev))
  ALLOCATE(COphotorate(ncrsA+ncrsB_low, nlev))
  
  ALLOCATE(heads(nsp+10),dheads(2*nsp+3))
  nlev = 241
  OPEN (unit=61,file='../runs/'//trim(runID)//'/output/atm1D.out',status='unknown',action='write')
  
     WRITE (61,"(2I6)") nlev, nsp
     WRITE (61,"(' MOLECULES')")
     WRITE (61,"(10(2X,A12,1X))") (name(nm),nm=1,nsp)
     WRITE (61,"(' ALTITUDE (cm)')") 
     WRITE (61,"(10ES15.7)") (z(nz), nz=1, nlev)
     WRITE (61,"(' RADIUS (cm)')") 
     WRITE (61,"(10ES15.7)") (rz(nz), nz=1, nlev)
     WRITE (61,"(' GRAVITY (cm s-2)')")
     WRITE (61,"(10ES15.7)") (grv(nz), nz=1, nlev)
     WRITE (61,"(' NEUTRAL TEMPERATURE (Kelvins)')") 
     WRITE (61,"(10ES15.7)") (tn(nz), nz=1, nlev)
     WRITE (61,"(' ELECTRON TEMPERATURE (Kelvins)')") 
     WRITE (61,"(10ES15.7)") (te(nz), nz=1, nlev)
     WRITE (61,"(' PRESSURE (dyne/cm^2)')") 
     WRITE (61,"(10ES15.7)") (prs(nz), nz=1, nlev)
     WRITE (61,"(' MASS DENSITY (g cm^-3)')") 
     WRITE (61,"(10ES15.7)") (rho(nz), nz=1, nlev)
     WRITE (61,"(' MEAN MOLECULAR WEIGHT (amu)')") 
     WRITE (61,"(10ES15.7)") (mass(nz), nz=1, nlev)
     WRITE (61,"(' EDDY COEFFICIENT (cm^2s^-1)')") 
     WRITE (61,"(10ES15.7)") (ek(nz), nz=1, nlev)
     WRITE (61,"(' TOTAL DENISTY (cm^-3)')")
     WRITE (61,"(10ES15.7)") (den(nz,0), nz=1, nlev) 
     DO nm = 1, nsp
        WRITE (61,"(A12,F5.2)") name(nm)
        WRITE (61,"(10ES15.7)") (den(nz,nm), nz=1, nlev)
     END DO
   
  CLOSE(unit=61)

  OPEN (unit=61,file='atm1D.out',status='unknown',action='write')
  
     WRITE (61,"(2I6)") nlev, nsp
     WRITE (61,"(' MOLECULES')")
     WRITE (61,"(10(2X,A12,1X))") (name(nm),nm=1,nsp)
     WRITE (61,"(' ALTITUDE (cm)')") 
     WRITE (61,"(10ES15.7)") (z(nz), nz=1, nlev)
     WRITE (61,"(' RADIUS (cm)')") 
     WRITE (61,"(10ES15.7)") (rz(nz), nz=1, nlev)
     WRITE (61,"(' GRAVITY (cm s-2)')")
     WRITE (61,"(10ES15.7)") (grv(nz), nz=1, nlev)
     WRITE (61,"(' NEUTRAL TEMPERATURE (Kelvins)')") 
     WRITE (61,"(10ES15.7)") (tn(nz), nz=1, nlev)
     WRITE (61,"(' ELECTRON TEMPERATURE (Kelvins)')") 
     WRITE (61,"(10ES15.7)") (te(nz), nz=1, nlev)
     WRITE (61,"(' PRESSURE (dyne/cm^2)')") 
     WRITE (61,"(10ES15.7)") (prs(nz), nz=1, nlev)
     WRITE (61,"(' MASS DENSITY (g cm^-3)')") 
     WRITE (61,"(10ES15.7)") (rho(nz), nz=1, nlev)
     WRITE (61,"(' MEAN MOLECULAR WEIGHT (amu)')") 
     WRITE (61,"(10ES15.7)") (mass(nz), nz=1, nlev)
     WRITE (61,"(' EDDY COEFFICIENT (cm^2s^-1)')") 
     WRITE (61,"(10ES15.7)") (ek(nz), nz=1, nlev)
     WRITE (61,"(' TOTAL DENSITY (cm^-3)')")
     WRITE (61,"(10ES15.7)") (den(nz,0), nz=1, nlev) 
     DO nm = 1, nsp
        WRITE (61,"(A12,F5.2)") name(nm)
        WRITE (61,"(10ES15.7)") (den(nz,nm), nz=1, nlev)
     END DO
   
  CLOSE(unit=61)

  heads(1) = 'ALT (km)'
  heads(2) = 'RAD (km)'
  heads(3) = 'GRV (cm/s2)'
  heads(4) = 'Tn (K)'
  heads(5) = 'Te (K)'
  heads(6) = 'PRS (ubar)'
  heads(7) = 'RHO (gm/cm3)'
  heads(8) = 'MMW (amu)'
  heads(9) = 'Kzz (cm2/s)'
  heads(10) = 'Ntot (cm-3)'
  DO nm = 1, nsp
     heads(10+nm) = name(nm)
  END DO

  OPEN (unit=61,file='../runs/'//trim(runID)//'/output/atm1D.csv',status='unknown',action='write')
     WRITE (61,"(A12,500(',',A12))") (heads(nm),nm=1,10+nsp)
     DO nz = 1, nlev
        WRITE (61,"(ES15.7,500(',',ES15.7))") 1.E-5*z(nz),1.E-5*rz(nz),grv(nz),tn(nz),te(nz),prs(nz),rho(nz),mass(nz),ek(nz),den(nz,0),(xmol(nz,nm),nm=1,nsp)
     END DO
  CLOSE(unit=61)

  !
  !  .. Diffusion Coefficients
  !

  dheads(1) = 'ALT (km)'
  dheads(2) = 'Kzz (cm^2s^-1)'
  DO nm = 1, nsp
     dheads(2+nm) = name(nm)
  END DO
  dheads(3+nsp) = 'HT (km)'
  DO nm = 1, nsp
     dheads(nsp+3+nm) = name(nm)
  END DO

  OPEN (unit=61,file='../runs/'//trim(runID)//'/output/diff.csv',status='unknown',action='write')
     WRITE (61,"(A12,500(',',A12))") (dheads(nm),nm=1,2*nsp+3)
     DO nz = 1, nlev
        WRITE (61,"(ES15.7,500(',',ES15.7))") 1.E-5*z(nz),ek(nz),(df(nz,nm),nm=1,nsp),(1.E-5_RP*ht(nz,nm),nm=0,nsp)
     END DO
  CLOSE(unit=61)


  !
  !  .. Solar Fluxes
  !

  OPEN(unit=62,file='../runs/'//TRIM(runID)//'/output/solar_flux.out',status='unknown')
     WRITE (62,"(2I6)") ncrsA+ncrsB_low-2+ncrsC, nlev
     WRITE (62,"('Wavelength (Angstroms)')") 
     wcrs(1:ncrsA) = wcrsA
     wcrs(ncrsA+1:ncrsA+ncrsB_low-2) = wcrsB_low(2:ncrsB_low-2)
     wcrs(ncrsA+ncrsB_low-1:ncrsA+ncrsB_low-2+ncrsC) = wcrsC
     WRITE (62,"(10ES11.3)") (wcrs(nw), nw=1, ncrsA+ncrsB_low-2+ncrsC)
     DO nz = 1, nlev
        WRITE(62,"(F11.3)") 1.E-5_RP*z(nz)
        DO nw = 1, ncrsA
           solar_flux(nw,nz) = fsolA(nw)*trnA(nw,nz)
        END DO
        DO nw = 1, ncrsB_low-2
           sm = zero
           DO nw_high = (nw-1)*50000+1, nw*50000
              sm = sm + fsolB(nw_high)*trnB(nw_high,nz)
           END DO
           solar_flux(ncrsA+nw,nz) = sm
        END DO
        DO nw = 1, ncrsC
           solar_flux(ncrsA+ncrsB_low-2+nw,nz) = fsolC(nw)*trnC(nw,nz)
        END DO
        WRITE(62,"(10ES11.3)") (solar_flux(nw,nz), nw=1,ncrsA+ncrsB_low-2+ncrsC)
     END DO
  CLOSE(unit=62)
  
  
  OPEN(unit=70,file='../runs/'//TRIM(runID)//'/output/COphotorates.out',status='unknown')
    WRITE (70,"(2I6)") ncrsA+ncrsB_low-2, nlev
    WRITE (70,"('Wavelength (Angstroms)')")
    WRITE (70,"(10ES11.3)") (wcrs(nw), nw=1, ncrsA+ncrsB_low-2)
    DO nz = 1, nlev
       WRITE(70,"(F11.3)") 1.E-5_RP*z(nz)
       na = 4
       nb = 4
       DO nw = 1, ncrsA
          COphotorate(nw,nz) = xcrsA(nw,na)*bratA(nw, nb, na)*fsolA(nw)*trnA(nw,nz)
       END DO
       na = 22
       nb = 1
       DO nw = 1, ncrsB_low-2
          sm = zero
          DO nw_high = (nw-1)*50000+1, nw*50000
             sm = sm + xcrsB(nw_high,na)*bratB(nw_high, nb, na)*fsolB(nw_high)*trnB(nw_high,nz)
          END DO
          COphotorate(ncrsA+nw,nz) = sm
       END DO
        WRITE(70,"(10ES11.3)") (COphotorate(nw,nz), nw=1,ncrsA+ncrsB_low-2)
     END DO
  CLOSE(unit=70)

  !
  !  .. Photolysis Rates
  !
 
  WHERE (rpt < 1.0E-99_RP) rpt = zero
  OPEN(unit=62,file='../runs/'//TRIM(runID)//'/output/photorates.out',status='unknown')
     WRITE (62,"(2I6)") nphrt, nlev
     WRITE (62,"(' ALTITUDE (km)')") 
     WRITE (62,"(10ES11.3)") (1.E-5*z(nz), nz=1, nlev)
     DO np = 1, nphrt
        WRITE(62,"(A)") ptitle(np)
        WRITE(62,"(10ES11.3)") (rpt(np,nz),nz=1,nlev)
     END DO
  CLOSE(unit=62)

  !
  !  .. Electron Fluxes
  !

  OPEN(unit=62,file='../runs/'//TRIM(runID)//'/output/eflux.out',status='unknown')
     WRITE (62,"(2I6)") nelb, nlev
     WRITE (62,"(' ENERGY GRID (eV)')") 
     WRITE (62,"(10ES11.3)") (elctreV(ne), ne=1, nelb)
     DO nz = 1, nlev
        WRITE(62,"(F11.3)") 1.E-5_RP*z(nz)
        WRITE(62,"(10ES11.3)") (eFLUX(nz,ne), ne=1,nelb)
     END DO
  CLOSE(unit=62)

  !
  !  .. Suprathermal Electron Rates
  !

  WHERE (rpe < 1.0E-99_RP) rpe = zero
  OPEN(unit=62,file='../runs/'//TRIM(runID)//'/output/elerates.out',status='unknown')
     WRITE (62,"(3I6)") nert, nlev
     WRITE (62,"(' ALTITUDE (km)')") 
     WRITE (62,"(10ES11.3)") (1.E-5*z(nz), nz=1, nlev)
     DO np = 1, nert
        WRITE(62,"(A)") etitle(np)
        WRITE(62,"(10ES11.3)") (rpe(np,nz),nz=1,nlev)
     END DO
  CLOSE(unit=62)

  !
  !  .. Chemical Rates
  !

!  nrmax = 0
!  DO nr = 1, nrct
!     IF(itype(nr) /= 1) nrmax = nrmax + 1
!  END DO

  WHERE (rt < 1.0E-99_RP) rt = zero
  OPEN(unit=63,file='../runs/'//TRIM(runID)//'/output/ratecoeff.out',status='unknown')
     WRITE (63,"(2I6)") nrct, nlev
     WRITE (63,"(' ALTITUDE (km)')") 
     WRITE (63,"(10ES11.3)") (1.E-5*z(nz), nz=1, nlev)
     DO nr = 1, nrct
!        IF(itype(nr) /= 1) THEN
        WRITE(63,"(A)") ctitle(nr)
        WRITE(63,"(10ES11.3)") (rt(nr,nz),nz=1,nlev)
!        END IF
     END DO
  CLOSE(unit=63)

  WHERE (rct < 1.0E-99_RP) rct = zero
  OPEN(unit=63,file='../runs/'//TRIM(runID)//'/output/chemrates.out',status='unknown')
     WRITE (63,"(2I6)") nrct, nlev
     WRITE (63,"(' ALTITUDE (km)')") 
     WRITE (63,"(10ES11.3)") (1.E-5*z(nz), nz=1, nlev)
     DO nr = 1, nrct
        WRITE(63,"(A)") ctitle(nr)
        WRITE(63,"(10ES11.3)") (rct(nr,nz),nz=1,nlev)
     END DO
  CLOSE(unit=63)

  !
  !  .. Column-integrated Rates for each reaction
  !


  !
  !   .. Photons
  !

  DO np = 1, nphrt
     sum0 = zero
     DO nz = nibot, nlev-1
        rfac=((half*(rz(nz+1)+rz(nz)))**2/RPLANET**2)*half*(rz(nz+1)-rz(nz))
        sum0=sum0 + rfac*(rpt(np,nz)+rpt(np,nz+1))
     END DO
     pcolrate(np,1) = sum0
     sum0 = zero
     DO nz = 1, nibot-1
        rfac=((half*(rz(nz+1)+rz(nz)))**2/RPLANET**2)*half*(rz(nz+1)-rz(nz))
        sum0=sum0 + rfac*(rpt(np,nz)+rpt(np,nz+1))
     END DO
     pcolrate(np,2) = sum0
     pcolrate(np,3) = pcolrate(np,1) + pcolrate(np,2)
  END DO

  OPEN(unit=74,file='../runs/'//TRIM(runID)//'/output/pcolrates.out',status='unknown',action='write')
     CALL INDEXX(pcolrate(1:nphrt,3),indxp)
     DO np = nphrt, 1, -1
        WRITE(74,"(I4,3(2X,ES11.3),2X,A87)") indxp(np),pcolrate(indxp(np),1), pcolrate(indxp(np),2),pcolrate(indxp(np),3),ptitle(indxp(np))
     END DO
  CLOSE(unit=74)

  !
  !  .. Suprathermal Electrons
  !

  DO np = 1, nert
     sum0 = zero
     DO nz = 1, nlev-1
        rfac=((half*(rz(nz+1)+rz(nz)))**2/RPLANET**2)*half*(rz(nz+1)-rz(nz))
        sum0=sum0 + rfac*(rpe(np,nz)+rpe(np,nz+1))
     END DO
     ecolrate(np) = sum0
  END DO

  OPEN(unit=74,file='../runs/'//TRIM(runID)//'/output/ecolrates.out',status='unknown',action='write')
     CALL INDEXX(ecolrate,indxp)
     DO np = nert, 1, -1
        WRITE(74,"(I4,2X,ES11.3,2X,A87)") indxp(np),ecolrate(indxp(np)), etitle(indxp(np))
     END DO
  CLOSE(unit=74)

  !
  !  .. Chemical Reactions
  !

  DO nr = 1, nrct
     sum0 = zero
     DO nz = nibot, nlev-1
        rfac=((half*(rz(nz+1)+rz(nz)))**2/RPLANET**2)*half*(rz(nz+1)-rz(nz))
        sum0=sum0 + rfac*(rct(nr,nz)+rct(nr,nz+1))
     END DO
     colrate(nr,1) = sum0
     sum0 = zero
     DO nz = 1, nibot-1
        rfac=((half*(rz(nz+1)+rz(nz)))**2/RPLANET**2)*half*(rz(nz+1)-rz(nz))
        sum0=sum0 + rfac*(rct(nr,nz)+rct(nr,nz+1))
     END DO
     colrate(nr,2) = sum0
     colrate(nr,3) = colrate(nr,1) + colrate(nr,2)
  END DO

  OPEN(unit=75,file='../runs/'//TRIM(runID)//'/output/colrates.out',status='unknown',action='write')
     CALL INDEXX(colrate(1:nrct,3),indxr)
     DO nr = nrct, 1, -1
        WRITE(75,"(I4,3(2X,ES11.3),2X,A73)") indxr(nr),colrate(indxr(nr),1),colrate(indxr(nr),2),colrate(indxr(nr),3),ctitle(indxr(nr)) 
     END DO
  CLOSE(unit=75)

   !
   !  .. Column-integrated density each species
   !

   DO nm = 1, nsp
      sum0=zero
      DO nl = 1, nlev-1
         rfac=((half*(rz(nl+1)+rz(nl)))**2/RPLANET**2)*half*(rz(nl+1)-rz(nl))
         sum0=sum0 + rfac*(den(nl,nm)+den(nl+1,nm))
      END DO
      cdens(nm) = sum0 
   END DO

   !
   !  .. Column-integrated rates for each species
   !

   DO nm = 1, nsp
      sum0=zero; sum1=zero; sum2=zero; sum3=zero; sum4=zero; sum5=zero; sum6=zero; sum7=zero
      DO nl = 1, nlev-1
         rfac=((half*(rz(nl+1)+rz(nl)))**2/RPLANET**2)*half*(rz(nl+1)-rz(nl))
         sum0=sum0 + rfac*(prext(nl,nm)+prext(nl+1,nm))
         sum1=sum1 + rfac*(pr_ph(nl,nm)+pr_ph(nl+1,nm))
         sum2=sum2 + rfac*(ls_ph(nl,nm)+ls_ph(nl+1,nm))
         sum3=sum3 + rfac*(pr_pe(nl,nm)+pr_pe(nl+1,nm))
         sum4=sum4 + rfac*(ls_pe(nl,nm)+ls_pe(nl+1,nm))
         sum5=sum5 + rfac*(pr_chem(nl,nm)+pr_chem(nl+1,nm))
         sum6=sum6 + rfac*(ls_chem(nl,nm)+ls_chem(nl+1,nm))
         sum7=sum7 + rfac*(rcdn(nl,nm)+rcdn(nl+1,nm))
      END DO
      cprd_ext(nm) = sum0
      cprd_ph(nm) = sum1
      clss_ph(nm) = -sum2
      cprd_pe(nm) = sum3
      clss_pe(nm) = -sum4
      cprd_chem(nm) = sum5
      clss_chem(nm) = -sum6
      cprd_net(nm) = sum0 + sum1 + sum3 + sum5
      clss_net(nm) = -sum2 - sum4 - sum6 
      crc(nm) = -sum7
      ctp(nm) = -flx(nlev,nm)*rz(nlev)**2/RPLANET**2
      cbt(nm) = flx(1,nm)
      cbl(nm) = cprd_net(nm)+clss_net(nm)+crc(nm)+ctp(nm)+cbt(nm)
      csf(nm) = -crc(nm) - cbt(nm)
      tbal(nm) = cdens(nm)/abs(cbl(nm))
   END DO

   !
   !  .. Summary file with column integrated balances for each molecule
   !

   CALL INDEXX(abs(cbl),indx)
   OPEN(unit=64,file='../runs/'//TRIM(runID)//'/output/balance.out',status='unknown')
      WRITE(64,"('    Name  Status  Density    Ext/GCR     Bot Flux    Top Flux      Prod       Loss       Condense     Balance  Time Const')")   
      DO nx = nsp, 1, -1
         nm = indx(nx)
         IF(istat(nm) > 0) THEN
         WRITE(64,"(A12,I3,10(2X,ES10.3))") name(nm),istat(nm),cdens(nm),cprd_ext(nm),cbt(nm),ctp(nm),cprd_net(nm),clss_net(nm), &
             crc(nm),cbl(nm),tbal(nm)
         END IF
      END DO
      DO nx = nsp, 1, -1
         nm = indx(nx)
         IF(istat(nm) <= 0) THEN
         WRITE(64,"(A12,I3,10(2X,ES10.3))") name(nm),istat(nm),cdens(nm),cprd_ext(nm),cbt(nm),ctp(nm),cprd_net(nm),clss_net(nm), &
             crc(nm),cbl(nm),tbal(nm)
         END IF
      END DO
   CLOSE(unit=64)


   !
   !  .. Sort
   !

   CALL INDEXX(csf,indx)
   OPEN(unit=63,file='../runs/'//TRIM(runID)//'/output/mcolrates.out',status='unknown',action='write')
      WRITE(63,"('                 Flux Bot     Condense     Flux Top        Prod         Loss      External      Balance')")
      DO nx = nsp, 1, -1
         nm = indx(nx)
         WRITE(63,"(A12,7(2X,ES11.3))") name(nm),cbt(nm),crc(nm),ctp(nm),cprd_net(nm),clss_net(nm),cprd_ext(nm),cbl(nm)
      END DO
   CLOSE(unit=63)

   !
   !  .. Species Files
   !

   WHERE (prext < 1.E-99_RP) prext = zero
   WHERE (pr_ph < 1.E-99_RP) pr_ph = zero
   WHERE (ls_ph < 1.E-99_RP) ls_ph = zero
   WHERE (pr_pe < 1.E-99_RP) pr_pe = zero
   WHERE (ls_pe < 1.E-99_RP) ls_pe = zero
   WHERE (pr_chem < 1.E-99_RP) pr_chem = zero
   WHERE (ls_chem < 1.E-99_RP) ls_chem = zero
   WHERE (pr < 1.E-99_RP) pr = zero
   WHERE (ls < 1.E-99_RP) ls = zero
   WHERE (rcdn < 1.E-99_RP) rcdn = zero

   DO nx = 1, SIZE(locp)           
      nm = locp(nx)
      OPEN(unit=64,file='../runs/'//TRIM(runID)//'/output/molecules/'//name(nm)(1:LEN_TRIM(name(nm)))//'.OUT',status='unknown')
         WRITE(64,"(ES11.3,' ;Column External Production  ')") cprd_ext(nm)
         WRITE(64,"(ES11.3,' ;Column Photon Production    ')") cprd_ph(nm)
         WRITE(64,"(ES11.3,' ;Column Photon Loss          ')") clss_ph(nm)
         WRITE(64,"(ES11.3,' ;Column Electron Production  ')") cprd_pe(nm)
         WRITE(64,"(ES11.3,' ;Column Electron Loss        ')") clss_pe(nm)
         WRITE(64,"(ES11.3,' ;Column chem Production      ')") cprd_chem(nm)
         WRITE(64,"(ES11.3,' ;Column chem Loss            ')") clss_chem(nm)
         WRITE(64,"(ES11.3,' ;Column Net Production       ')") cprd_net(nm)
         WRITE(64,"(ES11.3,' ;Column Net Loss             ')") clss_net(nm)
         WRITE(64,"(ES11.3,' ;Column Net Condensation     ')") crc(nm)
         WRITE(64,"(ES11.3,' ;Escape Flux                 ')") ctp(nm)
         WRITE(64,"(ES11.3,' ;Surface Flux                ')") cbt(nm)
         WRITE(64,"(ES11.3,' ;Net Balance                 ')") cbl(nm)
         WRITE(64,932)
         DO nl = 1, nlev
            WRITE(64,933) cm_to_km*z(nl), den(nl,nm), den(nl,nm)/den(nl,0), flx(nl,nm),             &
                prext(nl,nm), pr_ph(nl,nm), ls_ph(nl,nm), pr_pe(nl,nm), ls_pe(nl,nm),               &
                pr_chem(nl,nm), ls_chem(nl,nm), pr(nl,nm), ls(nl,nm), rcdn(nl,nm), div_flx(nl,nm),  &
                pr(nl,nm)-ls(nl,nm)-rcdn(nl,nm)-div_flx(nl,nm)
         END DO
     CLOSE(unit=64)
  END DO

   DO nx = 1, SIZE(ldcp)           
      nm = ldcp(nx)
      OPEN(unit=65,file='../runs/'//TRIM(runID)//'/output/molecules/'//name(nm)(1:LEN_TRIM(name(nm)))//'.OUT',status='unknown')
         WRITE(65,"(ES11.3,' ;Column External Production  ')") cprd_ext(nm)
         WRITE(65,"(ES11.3,' ;Column Photon Production    ')") cprd_ph(nm)
         WRITE(65,"(ES11.3,' ;Column Photon Loss          ')") clss_ph(nm)
         WRITE(65,"(ES11.3,' ;Column Electron Production  ')") cprd_pe(nm)
         WRITE(65,"(ES11.3,' ;Column Electron Loss        ')") clss_pe(nm)
         WRITE(65,"(ES11.3,' ;Column chem Production      ')") cprd_chem(nm)
         WRITE(65,"(ES11.3,' ;Column chem Loss            ')") clss_chem(nm)
         WRITE(65,"(ES11.3,' ;Column Net Production       ')") cprd_net(nm)
         WRITE(65,"(ES11.3,' ;Column Net Loss             ')") clss_net(nm)
         WRITE(65,"(ES11.3,' ;Column Net Condensation     ')") crc(nm)
         WRITE(65,"(ES11.3,' ;Escape Flux                 ')") ctp(nm)
         WRITE(65,"(ES11.3,' ;Surface Flux                ')") cbt(nm)
         WRITE(65,"(ES11.3,' ;Net Balance                 ')") cbl(nm)
         WRITE(65,932)
         DO nl = 1, nlev
            WRITE(65,933) cm_to_km*z(nl), den(nl,nm), den(nl,nm)/den(nl,0), flx(nl,nm),             &
                prext(nl,nm), pr_ph(nl,nm), ls_ph(nl,nm), pr_pe(nl,nm), ls_pe(nl,nm),               &
                pr_chem(nl,nm), ls_chem(nl,nm), pr(nl,nm), ls(nl,nm), rcdn(nl,nm), div_flx(nl,nm),  &
                pr(nl,nm)-ls(nl,nm)-rcdn(nl,nm)-div_flx(nl,nm)
         END DO
     CLOSE(unit=65)
  END DO

  !
  !  Summary files for atoms
  !

  !  .. Oxygen

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*den(nl,nm)
     END DO
     den_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     mol_oxy(nl) = den_oxy(nl)/den(nl,0)
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*flx(nl,nm)
     END DO
     flx_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*prext(nl,nm)
     END DO
     prext_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*pr_ph(nl,nm)
     END DO
     pr_ph_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*ls_ph(nl,nm)
     END DO
     ls_ph_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*pr_pe(nl,nm)
     END DO
     pr_pe_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*ls_pe(nl,nm)
     END DO
     ls_pe_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*pr_chem(nl,nm)
     END DO
     pr_chem_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*ls_chem(nl,nm)
     END DO
     ls_chem_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*pr(nl,nm)
     END DO
     pr_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*ls(nl,nm)
     END DO
     ls_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*rcdn(nl,nm)
     END DO
     rcdn_oxy(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + noxy(nm)*div_flx(nl,nm)
     END DO
     div_flx_oxy(nl) = sm
  END DO

  !

  !
  !  .. Column-integrated rates for O
  !

  sum0=zero; sum1=zero; sum2=zero; sum3=zero; sum4=zero; sum5=zero; sum6=zero; sum7=zero
  DO nl = 1, nlev-1
     rfac=((half*(rz(nl+1)+rz(nl)))**2/RPLANET**2)*half*(rz(nl+1)-rz(nl))
     sum0=sum0 + rfac*(prext_oxy(nl)+prext_oxy(nl+1))
     sum1=sum1 + rfac*(pr_ph_oxy(nl)+pr_ph_oxy(nl+1))
     sum2=sum2 + rfac*(ls_ph_oxy(nl)+ls_ph_oxy(nl+1))
     sum3=sum3 + rfac*(pr_pe_oxy(nl)+pr_pe_oxy(nl+1))
     sum4=sum4 + rfac*(ls_pe_oxy(nl)+ls_pe_oxy(nl+1))
     sum5=sum5 + rfac*(pr_chem_oxy(nl)+pr_chem_oxy(nl+1))
     sum6=sum6 + rfac*(ls_chem_oxy(nl)+ls_chem_oxy(nl+1))
     sum7=sum7 + rfac*(rcdn_oxy(nl)+rcdn_oxy(nl+1))
  END DO
  cprd_ext_oxy = sum0
  cprd_ph_oxy = sum1
  clss_ph_oxy = -sum2
  cprd_pe_oxy = sum3
  clss_pe_oxy = -sum4
  cprd_chem_oxy = sum5
  clss_chem_oxy = -sum6
  cprd_net_oxy = sum0 + sum1 + sum3 + sum5
  clss_net_oxy = -sum2 - sum4 - sum6 
  crc_oxy = -sum7
  ctp_oxy = -flx_oxy(nlev)*rz(nlev)**2/RPLANET**2
  cbt_oxy = flx_oxy(1)
  cbl_oxy = cprd_net_oxy+clss_net_oxy+crc_oxy+ctp_oxy+cbt_oxy
  
  OPEN(unit=65,file='../runs/'//TRIM(runID)//'/output/molecules/O_ATOMS.OUT',status='unknown')
     WRITE(65,"(ES11.3,' ;Column External Production  ')") cprd_ext_oxy
     WRITE(65,"(ES11.3,' ;Column Photon Production    ')") cprd_ph_oxy
     WRITE(65,"(ES11.3,' ;Column Photon Loss          ')") clss_ph_oxy
     WRITE(65,"(ES11.3,' ;Column Electron Production  ')") cprd_pe_oxy
     WRITE(65,"(ES11.3,' ;Column Electron Loss        ')") clss_pe_oxy
     WRITE(65,"(ES11.3,' ;Column chem Production      ')") cprd_chem_oxy
     WRITE(65,"(ES11.3,' ;Column chem Loss            ')") clss_chem_oxy
     WRITE(65,"(ES11.3,' ;Column Net Production       ')") cprd_net_oxy
     WRITE(65,"(ES11.3,' ;Column Net Loss             ')") clss_net_oxy
     WRITE(65,"(ES11.3,' ;Column Net Condensation     ')") crc_oxy
     WRITE(65,"(ES11.3,' ;Escape Flux                 ')") ctp_oxy
     WRITE(65,"(ES11.3,' ;Surface Flux                ')") cbt_oxy
     WRITE(65,"(ES11.3,' ;Net Balance                 ')") cbl_oxy
     WRITE(65,962)
     DO nl = 1, nlev
        WRITE(65,933) cm_to_km*z(nl), den_oxy(nl), mol_oxy(nl), flx_oxy(nl),                     &
             prext_oxy(nl), pr_ph_oxy(nl), ls_ph_oxy(nl), pr_pe_oxy(nl), ls_pe_oxy(nl),          &
             pr_chem_oxy(nl), ls_chem_oxy(nl), pr_oxy(nl), ls_oxy(nl), rcdn_oxy(nl),             &
             pr_oxy(nl)-ls_oxy(nl), div_flx_oxy(nl), pr_oxy(nl)-ls_oxy(nl)-rcdn_oxy(nl)-div_flx_oxy(nl)
     END DO
  CLOSE(unit=65)

  !
  !  Hydrogen
  !

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*den(nl,nm)
     END DO
     den_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     mol_hyd(nl) = den_hyd(nl)/den(nl,0)
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*flx(nl,nm)
     END DO
     flx_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*prext(nl,nm)
     END DO
     prext_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*pr_ph(nl,nm)
     END DO
     pr_ph_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*ls_ph(nl,nm)
     END DO
     ls_ph_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*pr_pe(nl,nm)
     END DO
     pr_pe_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*ls_pe(nl,nm)
     END DO
     ls_pe_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*pr_chem(nl,nm)
     END DO
     pr_chem_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*ls_chem(nl,nm)
     END DO
     ls_chem_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*pr(nl,nm)
     END DO
     pr_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*ls(nl,nm)
     END DO
     ls_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*rcdn(nl,nm)
     END DO
     rcdn_hyd(nl) = sm
  END DO

  DO nl = 1, nlev
     sm = zero
     DO nm = 1, nsp-1
        sm = sm + nhyd(nm)*div_flx(nl,nm)
     END DO
     div_flx_hyd(nl) = sm
  END DO

  !

  !
  !  .. Column-integrated rates for O
  !

  sum0=zero; sum1=zero; sum2=zero; sum3=zero; sum4=zero; sum5=zero; sum6=zero; sum7=zero
  DO nl = 1, nlev-1
     rfac=((half*(rz(nl+1)+rz(nl)))**2/RPLANET**2)*half*(rz(nl+1)-rz(nl))
     sum0=sum0 + rfac*(prext_hyd(nl)+prext_hyd(nl+1))
     sum1=sum1 + rfac*(pr_ph_hyd(nl)+pr_ph_hyd(nl+1))
     sum2=sum2 + rfac*(ls_ph_hyd(nl)+ls_ph_hyd(nl+1))
     sum3=sum3 + rfac*(pr_pe_hyd(nl)+pr_pe_hyd(nl+1))
     sum4=sum4 + rfac*(ls_pe_hyd(nl)+ls_pe_hyd(nl+1))
     sum5=sum5 + rfac*(pr_chem_hyd(nl)+pr_chem_hyd(nl+1))
     sum6=sum6 + rfac*(ls_chem_hyd(nl)+ls_chem_hyd(nl+1))
     sum7=sum7 + rfac*(rcdn_hyd(nl)+rcdn_hyd(nl+1))
  END DO
  cprd_ext_hyd = sum0
  cprd_ph_hyd = sum1
  clss_ph_hyd = -sum2
  cprd_pe_hyd = sum3
  clss_pe_hyd = -sum4
  cprd_chem_hyd = sum5
  clss_chem_hyd = -sum6
  cprd_net_hyd = sum0 + sum1 + sum3 + sum5
  clss_net_hyd = -sum2 - sum4 - sum6 
  crc_hyd = -sum7
  ctp_hyd = -flx_hyd(nlev)*rz(nlev)**2/RPLANET**2
  cbt_hyd = flx_hyd(1)
  cbl_hyd = cprd_net_hyd+clss_net_hyd+crc_hyd+ctp_hyd+cbt_hyd
  
  OPEN(unit=65,file='../runs/'//TRIM(runID)//'/output/molecules/H_ATOMS.OUT',status='unknown')
     WRITE(65,"(ES11.3,' ;Column External Production  ')") cprd_ext_hyd
     WRITE(65,"(ES11.3,' ;Column Photon Production    ')") cprd_ph_hyd
     WRITE(65,"(ES11.3,' ;Column Photon Loss          ')") clss_ph_hyd
     WRITE(65,"(ES11.3,' ;Column Electron Production  ')") cprd_pe_hyd
     WRITE(65,"(ES11.3,' ;Column Electron Loss        ')") clss_pe_hyd
     WRITE(65,"(ES11.3,' ;Column chem Production      ')") cprd_chem_hyd
     WRITE(65,"(ES11.3,' ;Column chem Loss            ')") clss_chem_hyd
     WRITE(65,"(ES11.3,' ;Column Net Production       ')") cprd_net_hyd
     WRITE(65,"(ES11.3,' ;Column Net Loss             ')") clss_net_hyd
     WRITE(65,"(ES11.3,' ;Column Net Condensation     ')") crc_hyd
     WRITE(65,"(ES11.3,' ;Escape Flux                 ')") ctp_hyd
     WRITE(65,"(ES11.3,' ;Surface Flux                ')") cbt_hyd
     WRITE(65,"(ES11.3,' ;Net Balance                 ')") cbl_hyd
     WRITE(65,962)
     DO nl = 1, nlev
        WRITE(65,933) cm_to_km*z(nl), den_hyd(nl), mol_hyd(nl), flx_hyd(nl),                     &
             prext_hyd(nl), pr_ph_hyd(nl), ls_ph_hyd(nl), pr_pe_hyd(nl), ls_pe_hyd(nl),          &
             pr_chem_hyd(nl), ls_chem_hyd(nl), pr_hyd(nl), ls_hyd(nl), rcdn_hyd(nl),             &
             pr_hyd(nl)-ls_hyd(nl), div_flx_hyd(nl), pr_hyd(nl)-ls_hyd(nl)-rcdn_hyd(nl)-div_flx_hyd(nl)
     END DO
  CLOSE(unit=65)


931 FORMAT(' Ext Prod = ',ES10.3,' Prod. = ',ES10.3,' Loss = ',ES10.3,                     &
         ' Condens. = ',ES10.3,' Escape Flux = ',ES10.3,' Surface Flux ',ES10.3,' Balance = ',ES10.3)
932 FORMAT('      alt     den        mole       flux       ext prd    pr_ph      ls_ph      ',      &
         'pr_pe      ls_pe     pr_chem   ls_chem    net prod   net loss    condense   -div_flx   balance')
933 FORMAT(F10.1,20(1X,ES10.3))

951 FORMAT(' Ext Prod = ',ES10.3,' Prod. = ',ES10.3,' Loss = ',ES10.3,                     &
         ' Condens. = ',ES10.3,' Escape Flux = ',ES10.3,' Surface Flux ',ES10.3,' Balance = ',ES10.3)
952 FORMAT('      alt     den        mole       flux       ext prd    prod       loss   ', &
 '   -div_flx    balance')
953 FORMAT(F10.1,10(1X,ES10.3))

962 FORMAT('      alt     den        mole       flux       ext prd    pr_ph      ls_ph      ',      &
         'pr_pe      ls_pe     pr_chem   ls_chem    net prod   net loss    condense    pr-ls    -div_flx   balance')

  ! ##########################################################################
  ! #                                                                        #
  ! #                           Solar Flux Stuff                             #
  ! #                                                                        #
  ! ##########################################################################

!  OPEN(unit=81,file='../runs/'//TRIM(runID)//'/output/opdepthC.out',status='unknown')
!     WRITE(81,'(2I6,F10.3)') nlev, nsol,cos_sza
!     WRITE(81,*) ' ALTITUDE (km)'
!     WRITE(81,'(10ES11.3)') (cm_to_km*z(nz),nz=1,nlev)
!     WRITE(81,*) ' PRESSURE (dynes/cm^2)'
!     WRITE(81,'(10ES11.3)') (prs(nz),nz=1,nlev)
!     WRITE(81,*) ' WAVELENGTH (Angstroms)'
!     WRITE(81,'(10ES11.3)') wsolC
!     WRITE(81,*) ' FLUX (photons cm-2 s-1)'
!     WRITE(81,'(10ES11.3)') fsolC
!     DO nw = 1, ncrsC
!        WRITE(81,'(F10.3)') wsolC(nw)
!        WRITE(81,*) ' TOTAL OPTICAL DEPTH'
!        WRITE(81,'(10ES11.3)') (trn(nw,nz),nz=1,nlev)
!        WRITE(81,*) ' GAS ABSORBER OPTICAL DEPTH'
!        WRITE(81,'(10ES11.3)') (tau_gas(nz,nw),nz=1,nlev)
!        WRITE(81,*) ' AEROSOl OPTICAL DEPTH'
!        WRITE(81,'(10ES11.3)') (tau_aer(nz,nw)/cos_sza,nz=1,nlev)
!        WRITE(81,*) ' RAYLEIGH OPTICAL DEPTH'
!        WRITE(81,'(10ES11.3)') (tau_ray(nz,nw)/cos_sza,nz=1,nlev)
!     END DO
!  CLOSE(unit=81)

  !
  !  .. Write output file with solar flux for plotting
  !

!  DO nz = 1, nlev
!     fb0(nz,1:5) = zero
!     DO nw = 1, nsol
!        IF(wsol(nw) < 1000.) THEN
!           fb0(nz,1) = fb0(nz,1) + fsol(nw)*EXP(-tau(nz,nw))
!        ELSE IF((wsol(nw)>=1000.).and.(wsol(nw)<1450.)) THEN
!           fb0(nz,2) = fb0(nz,2) + fsol(nw)*EXP(-tau(nz,nw))
!        ELSE IF((wsol(nw)>=1450.).and.(wsol(nw)<1680.)) THEN
!           fb0(nz,3) = fb0(nz,3) + fsol(nw)*EXP(-tau(nz,nw))
!        ELSE IF((wsol(nw)>=1680.).and.(wsol(nw)<2000.)) THEN
!           fb0(nz,4) = fb0(nz,4) + fsol(nw)*EXP(-tau(nz,nw))
!        ELSE IF(wsol(nw)>2000.) THEN
!           fb0(nz,5) = fb0(nz,5) + fsol(nw)*EXP(-tau(nz,nw))
!        END IF
!     END DO
!  END DO

!  OPEN(unit=82,file='../runs/'//TRIM(runID)//'/output/solbnd.out',status='unknown')
!     DO nz = 1, nlev
!        WRITE(82,'(F9.1,11(2X,ES11.3))') cm_to_km*z(nz),prs(nz),(fb0(nz,nw),nw=1,5)
!     END DO
!  CLOSE(unit=82)


!     OPEN(unit=40,file='eph.out',status='unknown')
!        DO npe = 1, nert
!           WRITE(40,"(I6,':  e + ',A12,' --> ',3(2X,A12))") npe,name(iele(1,npe)),(name(iele(np,npe)),np=2,4)
!           WRITE(40,"(10ES11.3)") (eph(npe,nl),nl=nibot,nlev)
!        END DO
!     CLOSE(unit=40)


!           OPEN(unit=30, file='electron_flux.out',status='unknown')
!              write(30,"(2I6)") nelb,nlev
!              write(30,"(' Energy Grid (eV)')")
!              write(30,"(10ES11.3)") (elctreV(ne),ne=1,nelb)
!              DO nl = nibot, nlev
!                 WRITE(30,"(' alt = ',ES11.3)") 1.E-5*z(nl)
!                 write(30,"(10ES11.3)") (eFLUX(nl,ne),ne=1,nelb)
!              END DO
!           CLOSE(30)

           
  DEALLOCATE(pcolrate,ecolrate,colrate,indxr,indxp,cprd_ext,cprd_ph,clss_ph,cprd_pe,clss_pe,        &
       cprd_chem,clss_chem,cprd_net,clss_net,crc,ctp,cbt,cbl,csf,wcrs,solar_flux,COphotorate)
  RETURN
END SUBROUTINE COMPOUT
