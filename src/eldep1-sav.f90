!###################################################################################################
!#                                                                                                 #
!# Calculation of electron mean intensity based on the photoelectron flux provided by soldep1.f90. #
!# The subroutine assumes no vertical transport and solves backwards for the mean intensity at     #
!# each energy bin assuming zero input for the highest energy bin.                                 #
!#                                                                                                 #
!###################################################################################################

SUBROUTINE ELDEP1

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_BIN, DIFCS_SEC

  !
  !  .. Declarations 
  !
 
  IMPLICIT NONE

  INTEGER :: nabs_el
  INTEGER :: i, na, j, k, m, n, ni, nm, npe, nb, nl, ne, nj
  INTEGER, DIMENSION(2) :: isp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: Selz, src_plasma
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: Ael, Ael_ext, Ael_sec, Ael_ion, sum_cs_ion, sum_cs_ext
  REAL(RP) :: sum1,sum2, sum_tot, sm, rn, rnp1, Ann, An1n1, E1, E2, E2_min, E2_max, step, Ei,Eion,  &
       sm_ion, sm_sec, E3, E4, tei, eden, sigma_excite_new, sigma_dif_ion, sigma_dif_sec, pSz,  &
       cs_ion, cs_sec
  CHARACTER(len=1) :: iret

  !
  !  .. Initialize
  ! 

  nabs_el = SIZE(ipath,1)
  eFLUX(:,:) = zero
  isp(1) = iCO2
  isp(2) = iN2

  ALLOCATE(Selz(nelb),src_plasma(nelb),Ael(nelb,nelb), Ael_ext(nelb,nelb), Ael_sec(nelb,nelb),      &
       Ael_ion(nelb,nelb),sum_cs_ion(nelb,nabs_el),sum_cs_ext(nelb,nabs_el))

  !
  !  .. Calculate electron flux at each energy and level
  !

  ALTI: DO i = nlev, nibot, -1

     Ael(:,:) = zero
     Ael_ext(:,:) = zero
     Ael_sec(:,:) = zero
     Ael_ion(:,:) = zero
     Selz(:)=Sel(i,:)
     tei = te(i)
     eden = den(i,iELE)


     IF( SUM(Selz(1:nelb)) > zero ) THEN

        sum_cs_ext(:,:) = zero 
        sum_cs_ion(:,:) = zero

        ENERGY: DO m = nelb,2,-1    

           SPECIES: DO na = 1, nabs_el

              ! Secondary electron production

              DO j = 1, ipath(na,3)               ! Loop over ionization channels

                 nj = ipath(na,1)+ipath(na,2)+j
                 E1 = elctreV(m)-enrgE(na,nj)
                 IF(E1 > zero) THEN 
                    ni = FIND_BIN(elctreV,elctDeV,E1)
                    IF((elctreV(ni) > E1) .and. (ni /= 1)) THEN
                       E4 = elctreV(ni)-half*elctDeV(ni)
                       ni =ni-1
                    ELSE
                       E4 = elctreV(ni)+half*elctDeV(ni)
                    END IF
                    Eion = elctreV(m) - E4
                    cs_ion = sum_cs_rees_ion(m,na,j)
                    cs_sec = sum_cs_rees_sec(m,na,j)
                    IF(ni > 0) THEN  
                       DO nm = ni, 1, -1  
                          IF(nm == ni) THEN
                             E3 = enrgE(na,nj)/Eion
                          ELSE
                             E3 = one
                          END IF
                          E2_min = elctreV(nm)-half*elctDeV(nm)
                          E2_max = elctreV(nm)+half*elctDeV(nm)
                          step = (E2_max-E2_min)/one
                          sm_ion = zero
                          sm_sec = zero
                          DO k=1,1
                             Ei = E2_min + (k-1)*step + half*step
                             sm_ion = sm_ion + step*DIFCS_SEC(elctreV(m)-Ei-Eion,elctreV(m),Eion)*eCS(m,na,nj)/cs_ion
                             sm_sec = sm_sec + step*DIFCS_SEC(Ei,elctreV(m),Eion)*eCS(m,na,nj)/cs_sec
                          END DO
                          sigma_dif_ion = sm_ion/elctDeV(nm)
                          sigma_dif_sec = sm_sec/elctDeV(nm)
                          Ael_ion(nm,m) = Ael_ion(nm,m) + den(i,isp(na))*elctDeV(m)*(sigma_dif_ion)*E3 
                          Ael_sec(nm,m) = Ael_sec(nm,m) + den(i,isp(na))*elctDeV(m)*(sigma_dif_sec)     
                       END DO
                    END IF
                    sum_cs_ion(m,na) = sum_cs_ion(m,na) + eCS(m,na,nj) 
                 END IF

              END DO

              ! Excitation to discrete states

              DO j = 1, ipath(na,1)+ipath(na,2)
                 sigma_excite_new = zero
                 E1 = elctreV(m)-enrgE(na,j)
                 nm = FIND_BIN(elctreV,elctDeV,E1)
                 sigma_excite_new = eCS(m,na,j)
                 IF((E1>zero) .and. (nm<m)) then
                    E1 = enrgE(na,j)/(elctreV(m)-elctreV(nm))
                    E3 = one
!                    sigma_excite_new = eCS(m,na,j) 
                    Ael_ext(nm,m) = Ael_ext(nm,m) + den(i,isp(na))*sigma_excite_new*E1
                    !              sigma_excite_new = eCS(m,na,j)*E1
                    !              Ael_ext(nm,m) = Ael_ext(nm,m) + den(i,isp(na))*sigma_excite_new
                 ENDIF
                 IF((E1>zero) .and. (nm==m)) then
                    E3 = enrgE(na,j)/(elctreV(nm)-elctreV(nm-1))
!                    sigma_excite_new = eCS(m,na,j) 

!  ### 5/13/15 RVY
!  ### Took a guess that /E3 was incorrect and should really be *E3 (based on *E1 above)
!  ### have a querry in to PL.  In the meantime, I tried a run with the change.  Results look better.

                    Ael_ext(nm-1,m) = Ael_ext(nm-1,m) + den(i,isp(na))*sigma_excite_new*E3
!                    Ael_ext(nm-1,m) = Ael_ext(nm-1,m) + den(i,isp(na))*sigma_excite_new/E3
                    !              sigma_excite_new = eCS(m,na,j)/E3 
                    !              Ael_ext(nm-1,m) = Ael_ext(nm-1,m) + den(i,isp(na))*sigma_excite_new
                 ENDIF
                 sum_cs_ext(m,na) = sum_cs_ext(m,na) + sigma_excite_new 
                 !           sum_cs_ext(m,na) = sum_cs_ext(m,na) + sigma_excite_new*E3
              END DO


           END DO SPECIES
        END DO ENERGY

!        OPEN(unit=61,file='Ael_ext.out',status='unknown')
!        WRITE(61,"(10ES11.3)") (elctreV(m),m=1,nelb)
!        DO m = 1, nelb
!           WRITE(61,"('E = ',ES11.3)") elctreV(m) 
!           WRITE(61,"(10ES11.3)") (Ael_ext(ne,m),ne=1,nelb)
!        END DO
!        CLOSE(unit=61)

!        OPEN(unit=61,file='Ael_sec.out',status='unknown')
!        WRITE(61,"(10ES11.3)") (elctreV(m),m=1,nelb)
!        DO m = 1, nelb
!           WRITE(61,"('E = ',ES11.3)") elctreV(m) 
!           WRITE(61,"(10ES11.3)") (Ael_sec(ne,m),ne=1,nelb)
!        END DO
!        CLOSE(unit=61)

!        OPEN(unit=61,file='Ael_ion.out',status='unknown')
!        WRITE(61,"(10ES11.3)") (elctreV(m),m=1,nelb)
!        DO m = 1, nelb
!           WRITE(61,"('E = ',ES11.3)") elctreV(m) 
!           WRITE(61,"(10ES11.3)") (Ael_ion(ne,m),ne=1,nelb)
!        END DO
!        CLOSE(unit=61)

!        OPEN(unit=61,file='esum_cs.out',status='unknown')
!        WRITE(61,"(' CO2 ext')")
!        WRITE(61,"(10ES11.3)") (sum_cs_ext(ne,1),ne=1,nelb)
!        WRITE(61,"(' CO2 ion')")
!        WRITE(61,"(10ES11.3)") (sum_cs_ion(ne,1),ne=1,nelb)
!        WRITE(61,"(' N2 ext')")
!        WRITE(61,"(10ES11.3)") (sum_cs_ext(ne,2),ne=1,nelb)
!        WRITE(61,"(' N2 ion')")
!        WRITE(61,"(10ES11.3)") (sum_cs_ion(ne,2),ne=1,nelb)
!        CLOSE(unit=61)

        Ael = Ael_ext + Ael_ion + Ael_sec

        ! Electron flux at each energy for level i

        eFLUX(i,nelb) = Selz(nelb)/(den(i,iCO2)*(sum_cs_ext(nelb,1)+sum_cs_ion(nelb,1)) &
             +den(i,iN2)*(sum_cs_ext(nelb,2)+sum_cs_ion(nelb,2)))

        DO n = nelb-1, 1, -1
           sum1 = DOT_PRODUCT(Ael(n,n+1:nelb),eFLUX(i,n+1:nelb))
           sum2 = DOT_PRODUCT(Ael(n+1,n+2:nelb),eFLUX(i,n+2:nelb))
           sum_tot = half * (sum1 + sum2)
           rn = eden*LOSS(elctreV(n),eden,tei)
           rnp1 = eden*LOSS(elctreV(n+1),eden,tei)
           Ann = (rn / (elctreV(n+1)-elctreV(n))) + half*( den(i,iCO2)* (sum_cs_ext(n,1) + sum_cs_ion(n,1) ) &
                + den(i,iN2)* (sum_cs_ion(n,2)   + sum_cs_ext(n,2)) )
           An1n1 = (rnp1 / (elctreV(n+1)-elctreV(n))) - half*( den(i,iCO2)* (sum_cs_ext(n+1,1) + sum_cs_ion(n+1,1)) & 
                + den(i,iN2)* (sum_cs_ion(n+1,2) + sum_cs_ext(n+1,2)) )
           IF (ABS(Ann) < 1.E-20_RP) THEN
              eFLUX(i,n) = zero
           ELSE
              eFLUX(i,n) = ( half*(Selz(n)+Selz(n+1)) + An1n1*eFLUX(i,n+1) + sum_tot) / Ann
              IF(eFLUX(i,n) > 1.E10) THEN
                 WRITE(*,*) i,n
                 WRITE(*,"(' sum_tot = ',ES12.4)") sum_tot
                 WRITE(*,"(' elctreV(n), elctreV(n+1) = ',2ES12.4)") elctreV(n), elctreV(n+1)
                 WRITE(*,"(' eFLUX(i,n), , eFLUX(i,n+1) = ',2ES12.4)") eFLUX(i,n), eFLUX(i,n+1)
                 WRITE(*,"(' Selz(n), Selz(n+1) = ',2ES12.4)") Selz(n),Selz(n+1)
                 WRITE(*,"(' Ann, An1n1 = ',2ES12.4)") Ann, An1n1
                 WRITE(*,"(' rn, rnp1 = ',2ES12.4)") rn, rnp1
                 WRITE(*,"(' first term = ',ES12.4)") half*(den(i,iCO2)*(sum_cs_ext(n,1) + sum_cs_ion(n,1) ) + den(i,iN2)* (sum_cs_ion(n,2)   + sum_cs_ext(n,2)) )
                 WRITE(*,"(' second term = ',ES12.4)") half*(den(i,iCO2)*(sum_cs_ext(n+1,1) + sum_cs_ion(n+1,1)) + den(i,iN2)* (sum_cs_ion(n+1,2) + sum_cs_ext(n+1,2)) )
                 WRITE(*,"(' sum_cs_ext(n,1), sum_cs_ion(n,1) =',2ES12.4)") sum_cs_ext(n,1), sum_cs_ion(n,1)
                 WRITE(*,"(' sum_cs_ext(n+1,1), sum_cs_ion(n+1,1) =',2ES12.4)") sum_cs_ext(n+1,1), sum_cs_ion(n+1,1)
                 WRITE(*,"(' sum_cs_ext(n,2), sum_cs_ion(n,2) = ',2ES12.4)") sum_cs_ext(n,2), sum_cs_ion(n,2)
                 WRITE(*,"(' sum_cs_ext(n+1,2), sum_cs_ion(n+1,2) =',2ES12.4)") sum_cs_ext(n+1,2), sum_cs_ion(n+1,2) 
                 WRITE(*,"(' den(i,iCO2), den(i,iN2) = ',2ES12.4)") den(i,iCO2), den(i,iN2)
                 STOP
              ENDIF
           ENDIF
           IF(eFLUX(i,n) <= zero ) eFLUX(i,n) = zero

        END DO

!        IF(i == nlev) THEN
!           OPEN(unit=61,file='Ael.out',status='unknown')
!           DO m = 1, nelb
!              WRITE(61,"(' ne = ',I6)") m
!              WRITE(61,"(10ES11.3)") (Ael(m,n),n=1,nelb)
!           END DO
!           CLOSE(unit=61)
!        ENDIF

     END IF
  
  END DO ALTI

!  OPEN(unit=60,file='eflux.out',status='unknown')
!     DO i = nibot, nlev
!        WRITE(60,"(' z = ',F11.3)") 1.E-5_RP*z(i)
!        WRITE(60,"(10ES11.3)") (eFLUX(i,n),n=1,nelb)
!     END DO
!  CLOSE(unit=60)


   
  ! Secondary production rates

  DO k = nibot, nlev
     pSz = zero
     do i = 1, nabs_el
        do n = 1, nelb  !-1
           pSz = pSz + den(k,isp(i))*eCS_ion(n,i)*eflux(k,n)*elctDeV(n) !0.5*(f1+f2)*de
        enddo
     enddo
     pS(k) = pSz
  ENDDO

  
  npe = 0
  DO na = 1, nabs_el_thk
     DO nb = 1, nbrnch_el(na)
        npe = npe + 1
        DO nl = nibot, nlev
           sm = zero
           DO ne = 1, nelb
              sm = sm + eflux(nl,ne)*eCS(ne,na,ipath(na,1)+nb)*elctDeV(ne) 
           END DO
           eph(npe,nl) = sm
        END DO
     END DO
  END DO
  
  DO n = 1, nabs_el_thn ! why is this different?
     na = n + nabs_el_thk
     DO nb = 1, nbrnch_el(na)
        npe = npe + 1
        DO nl = nibot, nlev
           sm = zero
           DO ne = 1, nelb
              sm = sm + brat_el(ne,nb,na)*eFLUX(nl,ne)*crs_tot_inel(ne,na)*elctDeV(ne) 
           END DO
           eph(npe,nl) = sm
        END DO
     END DO
  END DO

  DEALLOCATE(Selz,src_plasma,Ael,Ael_ext,Ael_sec,Ael_ion,sum_cs_ion,sum_cs_ext)

  RETURN

  CONTAINS

  FUNCTION LOSS(E,NE,TE)
    USE PRECISION
    USE CONSTANTS
    IMPLICIT NONE
    REAL(RP) :: LOSS, E, NE, TE, Ee
!    IF(E <= 0._RP) THEN
!	WRITE(*,"(' LOSS: E = ',ES15.7)") E
!        STOP
!    END IF
!    IF(Ne <= 0._RP) THEN
!	WRITE(*,"(' LOSS: Ne = ',ES15.7)") Ne 
!        STOP
!    END IF
    Ee = 8.618E-5_RP*TE
    IF (Ee > E) THEN
	LOSS = zero
    ELSE 
	LOSS = (3.37e-12_RP/((E**0.94_RP)*(NE**0.03_RP)))*((E-Ee)/(E-0.53_RP*Ee))**2.36_RP
    END IF
!    IF (Ee > E) LOSS = zero !(3.37e-12/((E**0.94)*(NE**0.03)))*((Ee-E)/(Ee-0.53*E))**2.36  !zero
!   IF (Ee > E) write(*,*) 'LOSS = ',LOSS
  END FUNCTION LOSS

END SUBROUTINE eldep1
