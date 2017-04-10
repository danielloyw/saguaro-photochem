!###################################################################################################
!#                                                                                                 #
!# Calculation of electron mean intensity based on the photoelectron flux provided by soldep1.f90. #
!# The subroutine assumes no vertical transport and solves backwards for the mean intensity at     #
!# each energy bin assuming zero input for the highest energy bin.                                 #
!#                                                                                                 #
!###################################################################################################

SUBROUTINE el_dep
!(nsp,Sel,den,te,elctreV,elctDeV,eCS_tot_inelast,eCS_ion,eCS_excite,eCS_ionize, & 
!     eCS_dissoc,state_excite,state_dissoc,state_ionize,enrg_excite,enrg_dissoc,enrg_ionize,ipath,eFLUX,pS) 

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL
  USE SUBS


  !
  !  .. Declarations 
  !
 
  IMPLICIT NONE

  !
  !  .. Local Variables
  !

  !INTEGER(IP) :: nlev
  !INTEGER(IP) :: nelb, nabs_el, nstep, nstep_max
  INTEGER(IP) :: i, is, j, ik, k, m, n, ne, ni, nl, nm , ns, nmin, nmax, kk, nstep,nstep_max, ism
  REAL(RP), DIMENSION(SIZE(elctreV)) :: Selz
  REAL(RP), DIMENSION(SIZE(elctreV)+1,SIZE(elctreV)+1) :: Ael, Ael_ext, Ael_sec, Ael_ion
  REAL(RP), DIMENSION(SIZE(elctreV)) :: SUMN
  REAL(RP), DIMENSION(SIZE(elctreV)) :: SUM_SEC
  REAL(RP), DIMENSION(SIZE(elctreV),nabs_el,5,20) :: sum_cs_rees_ion, sum_cs_rees_sec
  REAL(RP) :: SUMf, sum1,sum2, sum_tot, sm, rn, rnp1, Ann, An1n1, E1, E1_min,  &
       E1_max, E2, E2_min, E2_max, step, Ei,Eion, sm_ion, sm_sec,E3,E4,Emin,Emax
  REAL(RP) :: sumf_ion, sumf_sec
  REAL(RP), DIMENSION(SIZE(ipath,1)) :: ipl
  REAL(RP), DIMENSION(SIZE(elctreV)) :: src_plasma
  REAL(RP) :: eCStot,alpha, E_kk, E_min, E_mid
  REAL(RP) :: sum_enrg, sum_src,sm_n
  REAL(RP) :: sum_s_c, sum_t_e, sum_exc,sum_isj
  REAL(RP) :: error, sum_error, sum_er_tot,error2,sum_src_N
  REAL(RP) :: sigma_excite_new, sigma_dissoc_new, sigma_dif_ion, sigma_dif_sec,source_ph
  REAL(RP), DIMENSION(SIZE(elctreV),SIZE(ipath,1)) :: sum_cs_dis, sum_cs_ext, sum_cs_ion, cs_tot
  REAL(RP), DIMENSION(size(elctreV)+1) :: enrg
  REAL(RP), DIMENSION(size(elctreV)+1) :: flx,f_ion,f_sec,sm_f,csE,cs_new
  REAL(RP) :: E_pls, E_mns, Thml_Loss_pls, Thml_Loss_mns, ngas_loss_pls, ngas_loss_mns, F_pls, F_mns, gas_source,dene
  REAL(RP), DIMENSION(SIZE(elctreV)+1,SIZE(ipath,1),ipath(1,2)) :: dis_cs
  REAL(RP), DIMENSION(SIZE(elctreV)+1,SIZE(ipath,1),ipath(1,1)) :: ext_cs 


  !
  !  .. Parameters for gap-filling
  ! 
 
  REAL(RP) :: xmin, xn, x18, sigma, gauss, gauss18, delta, ratio

  !
  !  .. Parameters for secondary electron production
  
  REAL(RP) :: de, f1, f2, pSz

  !
  !  .. Initialize
  ! 

  !nlev = SIZE(te)
  !nelb = SIZE(elctreV)
  !nabs_el = SIZE(ipath,1)
  nstep_max = 1
  write(*,*) 'nelb=',nelb
  write(*,*) 'nabs_el=',nabs_el
  write(*,*) 'n_ion_states=',SIZE(ipath,2)
  write(*,*) ipath(1,1),ipath(1,2)
  ALLOCATE(eFLUX(nlev,nelb))
  enrg(:) = elctreV(:)

!------------------------------------------------------------
! IONIZATION POTENTIALS FOR SECONDARY ELECTRON CALCULATIONS |
!------------------------------------------------------------
  ipl(1) = 15.4_RP
  ipl(2) = 13.595_RP
  ipl(3) = 24.581_RP
  ipl(4) = 15.6_RP
  ipl(5) = 12.6_RP


!------------------------------------------------------------
!  .. Calculation of normalization factors for differential |
!  .. cross section of e-impact. This is only for the       |
!  .. ionization processes where a new electron emerges.    |
!------------------------------------------------------------
  OPEN(53,file='Normalization')
  DO is = 1, nabs_el
     ism  = loab_el(is)      ! This is the species loop
     write(53,'("Species: ", A10,1x,"Ionization states: ", i3)') name(ism), ipath(is,3)
     Do ik = 1, ipath(is,3)  ! This is the ionization processes loop for this species.
        write(53,'("State: ",A10,1x,"IP: ",f8.3)') state_ionize(is,ik), enrg_ionize(is,ik)
        DO i = 1, nelb       ! This is the electron energy loop
           nstep=nstep_max   ! This is the number of sub-bins in which each energy bin will be divided  
           if(enrg_ionize(is,ik) > elctDeV(i)) nstep=nstep_max
           dE = elctDeV(i)/nstep
           E_min = elctreV(i)-0.5_RP*elctDeV(i)
           kk = 1
           do kk = 1, nstep  ! This is the energy loop for each sub-bin
              E_kk = E_min + dE*kk - 0.5_RP*dE
              if(nstep == 1) E_kk = enrg(i)
              E1 = E_kk-enrg_ionize(is,ik) !E1 is the energy of degraded electron if we start from sub-bin kk
              IF(E1 > zero) then 
                 ni = FIND_BIN(elctreV,elctDeV,nelb,E1)
                 IF(ni < i) then 
!                 if(E1 < elctreV(ni) .and. ni /= 1 ) ni = ni - 1
                 ELSEIF(ni == i) then 
                    ni = ni-1
                 ENDIF
                 Eion = E_kk -  enrg(ni)
                 f_ion(:) = zero
                 f_sec(:) = zero
                 DO j = 1,ni
                    E3 = enrg(j) + Eion
                    E2 = E_kk - enrg(j)
                    IF(E3 > E_kk) E3 = E_kk
                    f_ion(j) = DIFCS_SEC(E_kk-enrg(j)-Eion,E_kk,Eion)
                    f_sec(j) = DIFCS_SEC(enrg(j),E_kk,Eion)
!!$                    f_ion(j) = DIFCS_REES(E2,E_kk)
!!$                    f_sec(j) = DIFCS_REES(E3,E_kk)
                 ENDDO
                 sumf_ion = zero
                 sumf_sec = zero
                 Do j = 1, ni
                    sumf_ion = sumf_ion + f_ion(j)*elctDeV(j) 
                    sumf_sec = sumf_sec + f_sec(j)*elctDeV(j) 
                 Enddo
                 if(sumf_ion == zero) sumf_ion = one
                 if(sumf_sec == zero) sumf_sec = one
                 sum_cs_rees_ion(i,is,ik,kk) = sumf_ion
                 sum_cs_rees_sec(i,is,ik,kk) = sumf_sec
                 write(53,'(3(i3,1x),3(F8.2,1x),2(E12.5,1x)))') & 
                      i,kk,ni,enrg(i),enrg(ni),eion, &
                      sum_cs_rees_ion(i,is,ik,kk), &
                      sum_cs_rees_sec(i,is,ik,kk)
              ENDIF               
           enddo
        ENDDO
     ENDDO
  ENDDO
  CLOSE(53)
!!$ write(*,*) 'Pass Normalization'

!-------------------------------------------------------
!  .. Calculate electron flux at each energy and level |
!-------------------------------------------------------
  open(96,file='energy-balance.txt')
  ALTI: DO i = 1, nlev
     Ael(:,:) = zero
     Ael_ext(:,:) = zero
     Ael_sec(:,:) = zero
     Ael_ion(:,:) = zero
     Selz(:) = Sel(i,:)
     !Selz(nelb) = 1._RP
     dene = den_el(i)
     write(*,'(i3,6(Es10.3,1x))') i, dene,den(i,1),den(i,2),den(i,3),den(i,4),den(i,5)
     sum_cs_ext(:,:) = zero
     sum_cs_dis(:,:) = zero
     sum_cs_ion(:,:) = zero

     DO n = nelb,2,-1
        DO is = 1, nabs_el
           ism = loab_el(is)

!--------------------------------------------------------
! .. Ionization of neutrals & production of secondaries |
!--------------------------------------------------------
           Do j = 1, ipath(is,3)
              nstep=nstep_max 
              if(enrg_dissoc(is,j) > elctDeV(n)) nstep=nstep_max
              dE = elctDeV(n)/nstep
              E_min = elctreV(n)-0.5_RP*elctDeV(n)
              do kk = 1, nstep
                 E_kk = E_min + dE*kk - 0.5_RP*dE
                 if(nstep == 1) E_kk = enrg(n)

                 E1 = E_kk-enrg_ionize(is,j)
                 IF(E1 > zero) then 
                    E2 = one
                    ni = FIND_BIN(elctreV,elctDeV,nelb,E1)
                    IF(ni < n) then 
                       m = ni
!                    if(E1 < elctreV(ni) .and. ni /= 1) m = ni - 1
                    ELSEIF(ni == n) then 
                       m = ni-1
                    ENDIF
                    Eion = E_kk - enrg(m)
                    cs_new(:) = zero
                    DO k = 1, m   !<--- loop for energy of secondary
                       E3 = enrg(k) + Eion   !E3=W=Esecondary + Ionization energy
                       IF(E3 <= E_kk) then 
                          if(ni < n-1) E2 = ((enrg_ionize(is,j)+enrg(k))/(E_kk-enrg(m-k+1)))
                          if(ni == n-1) E2 = ((enrg_ionize(is,j)+enrg(k))/(E_kk-enrg(m-k+1)))
!!$    sigma = DIFCS_REES(E_kk-enrg(k),E_kk) * eCS_dissoc(n,is,j) / (nstep*SUM_CS_REES_ion(n,is,j,kk))
                          sigma = DIFCS_SEC(E_kk-enrg(k)-Eion,E_kk,Eion) * eCS_ionize(n,is,j) /&
                               (nstep*SUM_CS_REES_ion(n,is,j,kk))
                          cs_new(k) = sigma * E2
                          Ael_ion(k,n) = Ael_ion(k,n) + sigma * den(i,ism) * E2 
!!$    sigma = DIFCS_REES(E3,E_kk) * eCS_dissoc(n,is,j) / (nstep*SUM_CS_REES_sec(n,is,j,kk))
                          sigma = DIFCS_SEC(enrg(k),E_kk,Eion) * eCS_ionize(n,is,j) / &
                               (nstep*SUM_CS_REES_sec(n,is,j,kk))
                          Ael_sec(k,n) = Ael_sec(k,n) + sigma * den(i,ism) 
                       ENDIF
!!$                       IF(E3 > E_kk) then  !if the primary is totally degraded only a secondary emerges
!!$                          E3 = E_kk
!!$                          sigma = zero
!!$                          if(SUM_CS_REES_sec(n,is,j,kk) /= zero) &
!!$!                               sigma = DIFCS_REES(E3,E_kk) * eCS_dissoc(n,is,j) / (nstep*SUM_CS_REES_sec(n,is,j,kk))
!!$                               sigma = DIFCS_SEC(enrg(k),E_kk,Eion) * eCS_dissoc(n,is,j) / (nstep*SUM_CS_REES_sec(n,is,j,kk))
!!$                          if(SUM_CS_REES_sec(n,is,j,kk) == zero) &
!!$                               sigma = eCS_dissoc(n,is,j)
!!$                          Ael_sec(k,n) = Ael_sec(k,n) + sigma * den(i,is)
!!$                       ENDIF
                    ENDDO
                    ! Calculate the effective ionization cross section due to grid formulation
                    sm=zero
                    DO k = 1,m
                       sm = sm + cs_new(k)*elctDeV(k)
                    ENDDO
                    sum_cs_ion(n,is) = sum_cs_ion(n,is) + sm
                 ENDIF
              enddo
           EndDo

!------------------------------------------------------------
! .. Production of energy E electrons from the degradation  |
! .. of higher energy electrons by neutral gas excitation   |
!------------------------------------------------------------
           Do j = 1, ipath(is,1)
              nstep=nstep_max 
              if(enrg_excite(is,j) > elctDeV(n)) nstep=nstep_max
              dE = elctDeV(n)/nstep
              E_min = elctreV(n)-0.5_RP*elctDeV(n)
              do kk = 1, nstep
                 E_kk = E_min + dE*kk - 0.5_RP*dE
                 if(nstep == 1) E_kk = enrg(n)
                 E1 = E_kk - enrg_excite(is,j)
                 IF(E1 > zero) then 
                    ni =  FIND_BIN(elctreV,elctDeV,nelb,E1)
                    IF(ni < n) then 
                       m = ni
!                    if(E1 < elctreV(ni) .and. ni /= 1) m = ni - 1
                       E2 = (enrg_excite(is,j))/(E_kk-enrg(m))
                       sigma_excite_new = eCS_excite(n,is,j)*E2/nstep
                       Ael_ext(m,n) = Ael_ext(m,n) + den(i,ism)*sigma_excite_new 
                    ELSEIF(ni == n) then 
                       E2 = ((enrg_excite(is,j))/(E_kk-enrg(ni-1)))
                       sigma_excite_new = eCS_excite(n,is,j)*E2/nstep
                       Ael_ext(ni-1,n) = Ael_ext(ni-1,n) + den(i,ism)*sigma_excite_new 
                    ENDIF
                    sum_cs_ext(n,is) = sum_cs_ext(n,is) + sigma_excite_new
                 ENDIF
              enddo
           EndDo

!------------------------------------------------------------
! .. Production of energy E electrons from the degradation  |
! .. of higher energy electrons by neutral gas dissociation |
!------------------------------------------------------------
           Do j = 1, ipath(is,2)
              nstep=nstep_max 
              if(enrg_dissoc(is,j) > elctDeV(n)) nstep=nstep_max
              dE = elctDeV(n)/nstep
              E_min = elctreV(n)-0.5_RP*elctDeV(n)
              do kk = 1, nstep
                 E_kk = E_min + dE*kk - 0.5_RP*dE
                 if(nstep == 1) E_kk = enrg(n)
                 E1 = E_kk - enrg_dissoc(is,j)
                 IF(E1 > zero) then 
                    ni =  FIND_BIN(elctreV,elctDeV,nelb,E1)
                    IF(ni < n) then 
                       m = ni
!                    if(E1 < elctreV(ni) .and. ni /= 1) m = ni - 1
                       E2 = (enrg_dissoc(is,j))/(E_kk-enrg(m))
                       sigma_dissoc_new = eCS_dissoc(n,is,j)*E2/nstep
                       Ael_ext(m,n) = Ael_ext(m,n) + den(i,ism)*sigma_dissoc_new 
                    ELSEIF(ni == n) then 
                       E2 = ((enrg_dissoc(is,j))/(E_kk-enrg(ni-1)))
                       sigma_dissoc_new = eCS_dissoc(n,is,j)*E2/nstep
                       Ael_ext(ni-1,n) = Ael_ext(ni-1,n) + den(i,ism)*sigma_dissoc_new 
                    ENDIF
                    sum_cs_dis(n,is) = sum_cs_dis(n,is) + sigma_dissoc_new
                 ENDIF
              enddo
           EndDo
!--------------FINISHED CALCULATING SOURCE TERMS FOR EACH ENERGY-------------------------

        ENDDO
     ENDDO

     Ael = Ael_ion + Ael_sec

     DO is = 1,nabs_el
        DO n = 1, nelb
           CS_tot(n,is) = sum_cs_dis(n,is)+sum_cs_ext(n,is)+sum_cs_ion(n,is)
        ENDDO
     ENDDO

!------------------------------------------------------------
! CALCULATE FLUXES AT EACH ENERGY BIN FOR SPECIFIC ALTITUDE |
!------------------------------------------------------------
     flx(:) = zero
     DO n = nelb-1, 1, -1
        E_pls = elctreV(n+1) + 0.5_RP*elctDeV(n+1)
        E_mns = elctreV(n) - 0.5_RP*elctDeV(n)
        if(n==1) E_mns = 0.05_RP
! .. Thermal e loss
        dE = E_pls-E_mns
        E_mid = elctreV(n+1) - 0.5_RP*elctDeV(n+1)
        Thml_Loss_pls = dene*Loss(enrg(n+1),dene,te(i))/(enrg(n+1)-enrg(n)) 
!Loss(E_pls,dene,te(i))/dE  !  (Loss(E_pls,dene,te(i)) - Loss(E_mid,dene,te(i)))/dE   !
        E_mid = elctreV(n) + 0.5_RP*elctDeV(n)
        Thml_Loss_mns = dene*Loss(enrg(n),dene,te(i))/(enrg(n+1)-enrg(n))  
!Loss(E_mns,dene,te(i))/dE  !(Loss(E_mid,dene,te(i)) - Loss(E_mns,dene,te(i)))/dE   !
! .. Neutral gas loss
        ngas_loss_pls = elctDev(n+1)*SUM(den(i,loab_el(1:nabs_el))*CS_tot(n+1,1:nabs_el))/dE
        ngas_loss_mns = elctDeV(n)*SUM(den(i,loab_el(1:nabs_el))*CS_tot(n,1:nabs_el))/dE
! .. Production by collisions with neutrals leading to ionization & secondaries
        F_pls = zero
        F_mns = zero
        do k = n+2,nelb
           F_pls = F_pls + Ael(n+1,k)*flx(k)*elctDeV(k)
        enddo
        do k = n+1,nelb
           F_mns = F_mns + Ael(n,k)*flx(k)*elctDeV(k)
        enddo
        gas_source = (elctDeV(n+1)*F_pls + elctDeV(n)*F_mns)/dE
! .. Production by collisions with neutrals leading to excitation
        F_pls = zero
        F_mns = zero
        do k = n+2,nelb
           F_pls = F_pls + Ael_ext(n+1,k)*flx(k)
        enddo
        do k = n+1,nelb
           F_mns = F_mns + Ael_ext(n,k)*flx(k)
        enddo
        gas_source = gas_source + (elctDeV(n+1)*F_pls + elctDeV(n)*F_mns)/dE
! .. Source term from photoelectrons
        source_ph = (elctDeV(n)*selz(n) + elctDeV(n+1)*selz(n+1))/dE
! .. Electron flux at the boundaries of each bin
        flx(n) = ( source_ph + gas_source + Thml_loss_pls*flx(n+1) - ngas_loss_pls*flx(n+1) ) & 
             / ( Thml_loss_mns + ngas_loss_mns ) 
        if(flx(n) < zero) flx(n) = zero
!        write(*,'(8(E12.5,1x))') enrg(n),selz(n),gas_source,Thml_loss_pls,Thml_loss_mns, &
!             ngas_loss_pls,ngas_loss_mns,flx(n)
     ENDDO
! .. Electron flux at the midpoint of each bin
     eFLUX(i,1:nelb) = flx(1:nelb)

!------------------------------
! Energy balance calculations |
!------------------------------

!!$     sum_error = zero
!!$     DO n = 1,1 ! nelb
!!$
!!$        sm_f(:) = zero
!!$        do m = 1, n
!!$           do k = n+1,nelb
!!$             sm_f(m) = sm_f(m) + Ael(m,k)*flx(k)*elctDeV(k) + &
!!$                                 Ael_ext(m,k)*flx(k)
!!$           enddo
!!$        enddo
!!$        sum_s_c = zero
!!$        do m = 1, n
!!$           sum_s_c = sum_s_c + sm_f(m)*enrg(m)*elctDeV(m)
!!$        enddo
!!$
!!$        sum_t_e = zero
!!$        do m = n+1,nelb
!!$           sum_t_e = sum_t_e + flx(m)*Loss(enrg(m),dene,te(i))*dene*elctDeV(m)
!!$        enddo
!!$
!!$        csE(:) = zero
!!$        do m = n+1,nelb
!!$           sum_isj = zero
!!$           do is=1,nabs_el 
!!$              do j =1,ipath(is,1) 
!!$                 sum_isj = sum_isj + den(i,is)*ecs_excite(m,is,j)*enrg_excite(is,j)
!!$              enddo
!!$              do j =1,ipath(is,2) 
!!$                 sum_isj = sum_isj + den(i,is)*ecs_dissoc(m,is,j)*enrg_dissoc(is,j)
!!$              enddo
!!$           enddo
!!$           csE(m) = sum_isj
!!$        enddo
!!$
!!$        sum_exc = zero
!!$        do m = n+1, nelb
!!$           sum_exc = sum_exc + csE(m)*flx(m)*elctDeV(m)
!!$        enddo
!!$        
!!$        sum_enrg = sum_s_c + sum_t_e + sum_exc 
!!$
!!$        sum_src = zero
!!$        SUM_SRC_N = ZERO
!!$        do k=n+1,nelb                       
!!$           sum_src = sum_src + Selz(k)*enrg(k)*elctDeV(k)
!!$           sum_src_N = sum_src_N + Selz(k)*elctDeV(k)
!!$        enddo
!!$        error = 0.
!!$        if(sum_src /= 0.) error = 100.*(sum_enrg-sum_src)/sum_src
!!$
!!$        sm_n = 0.
!!$        do j = 1,n
!!$           do k = n+1,nelb
!!$              sm_n = sm_n + (Ael_ion(j,k) + Ael_ext(j,k) )*flx(k)*elctDeV(j)
!!$           enddo
!!$        enddo
!!$        error2 = 100.*(sum_src_N-sm_n)/sum_src_N
!!$        if(n.eq.1) write(97,'(i3,1x,f12.5,1x,7(es11.3,1x),f8.3)') i, sum_s_c, sum_exc,  & 
!!$         sum_t_e, sum_enrg, sum_src, Sum(Selz(n:nelb)), eFLUX(i,n), error
!!$       if(i.eq.1.and.n.eq.1) write(96,'(a3,1x,a12,1x,7(a11,1x),2x,a8,1x,2(A11,1x),a8)')      &
!!$            ' ','Energy','PhotoEl','FLUX','Sec+Ion','Excit','Thrm.El','Total','Primary','Error','N_prim','N_c','error'
!!$       if(i.eq.1) write(96,'(i3,1x,f12.5,1x,7(es11.3,1x),2x,f8.3,1x,2(es11.3,1x),f8.3)')       & 
!!$             n, enrg(n),selz(n),flx(n),sum_s_c, sum_exc, sum_t_e, sum_enrg, sum_src, error, sum_src_N,sm_n, error2
!!$     sum_error = sum_error + abs(error)/nelb
!!$     ENDDO
!!$!     write(96,'(i3,1x,f12.5,1x,f12.5)') i, z(i)*1.E-5, sum_error  !/nelb
!!$     sum_er_tot = sum_er_tot + (sum_error)/nlev
   
  END DO ALTI
!!$  close(50)
!!$  write(96,*) 'Total error = ',sum_er_tot
!!$  close(96)
  write(*,*) 'pass calculations'

!**************************************
! Secondary electron production rates *
!**************************************
  ALLOCATE(pS(nlev))
  DO k = 1, nlev
     pSz = 0.
     do i = 1, nabs_el
        do n = 2, nelb  !-1
           pSz = pSz + den(k,i)*eCS_ion(n,i)*eflux(k,n)*elctDeV(n) !0.5*(f1+f2)*de
        enddo
     enddo
     pS(k) = pSz
  ENDDO
  write(*,*) 'PASS secondary'

!-----------------
! OUTPUT RESULTS |
!-----------------
  OPEN(20,file='eFLUXES.out')
  write(20,'(2(i3,1x))') nlev,nelb
  write(20,'(200(ES11.3,1x))') enrg
  DO k = 1, nlev 
     write(20,'(f10.3,1x,201(ES11.3,1x))') z(k),ps(k),(eFLUX(k,n),n=1,nelb)
  ENDDO
  CLOSE(20)


!!$open(30,file='eflux_1010.out')
!!$write(30,'(A15,1x,A15)') 'Energy','Flux'
!!$write(30,'(A15,1x,A15)') '(eV)','(elec.cm-2s-1)'
!!$do ne=1,nelb+1
!!$   write(30,'(2(ES15.5,1x,ES15.5))') enrg(ne),flx(ne),eflux(1,ne)
!!$enddo
!!$close(30)

  RETURN
END SUBROUTINE el_dep
