!###################################################################################################
!#                                                                                                 #
!# Calculation of electron mean intensity based on the photoelectron flux provided by soldep1.f90. #
!# The subroutine assumes no vertical transport and solves backwards for the mean intensity at     #
!# each energy bin assuming zero input for the highest energy bin.                                 #
!#                                                                                                 #
!###################################################################################################


! #
! #  Remember to check extrapolation of e + CO2 --> CO2+(A,B,C) to higher energies.
! #

! #
! #  Is DIFCS_SEC specific to N2 or does it work for CO2 as well?
! #

SUBROUTINE ELDEP1

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_BIN, DIFCS_SEC

  !
  !  .. Declarations 
  !
 
  IMPLICIT NONE

  INTEGER :: i, na, j, n, ni, nm, npe, nb, nl, ne, nj
  INTEGER, DIMENSION(2) :: isp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: Selz, src_plasma
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: Ael, Ael_ext, Ael_sec, Ael_ion, sum_cs_ion, sum_cs_ext
  REAL(RP) :: sum1,sum2, sum_tot, sm, rn, rnp1, Ann, An1n1, E1, E2, Eion,  &
       sm_ion, sm_sec, E3, E4, eden, sigma_scaled, sigma_dif_ion, sigma_dif_sec, pSz,  &
       cs_ion, cs_sec
  CHARACTER(len=1) :: iret

  !
  !  .. Initialize
  ! 

  eFLUX(:,:) = zero
  isp(1) = iCO2
  isp(2) = iN2

  ALLOCATE(Selz(nelb),src_plasma(nelb),Ael(nelb,nelb), Ael_ext(nelb,nelb), Ael_sec(nelb,nelb),      &
       Ael_ion(nelb,nelb),sum_cs_ion(nelb,nabs_el_thk),sum_cs_ext(nelb,nabs_el_thk))

  !
  !  .. Calculate electron flux at each energy and level
  !

  DO nl = nibot, nlev

     Ael(:,:) = zero             ! total attenuation
     Ael_ext(:,:) = zero         ! attenuation by excitation of neutrals
     Ael_sec(:,:) = zero         ! attenuation by secondary electron
     Ael_ion(:,:) = zero         ! attenuation by ionization
     Selz(:)=Sel(nl,:)
     eden = den(nl,iELE)

     IF( SUM(Selz(1:nelb)) > zero ) THEN
        sum_cs_ext(:,:) = zero ! total excitation cross-section
        sum_cs_ion(:,:) = zero ! total ionization cross-section

        DO ne = 2, nelb
           DO na = 1, nabs_el_thk
              
              ! Excitation to discrete states

              DO j = 1, ipath(na,1)+ipath(na,2)
                 E1 = elctreV(ne)-enrgE(na,j) ! exiting energy of electron
                 IF(E1 > zero) THEN
                    nm = FIND_BIN(elctreV,elctDeV,E1)
                    IF(nm < ne) THEN
                       sigma_scaled = eCS(ne,na,j)*(enrgE(na,j)/(elctreV(ne)-elctreV(nm)))
                       Ael_ext(nm,ne) = Ael_ext(nm,ne) + den(nl,isp(na))*sigma_scaled
                    ELSE IF(nm == ne) THEN
                       sigma_scaled = eCS(ne,na,j)*(enrgE(na,j)/(elctreV(nm)-elctreV(nm-1)))
                       Ael_ext(nm-1,ne) = Ael_ext(nm-1,ne) + den(nl,isp(na))*sigma_scaled
                    ENDIF
                    sum_cs_ext(ne,na) = sum_cs_ext(ne,na) + sigma_scaled
                 END IF
              END DO
			  
			  ! Ionization channels
		
			  DO j = 1, ipath(na,3)
                 nj = ipath(na,1)+ipath(na,2)+j
                 E1 = elctreV(ne)-enrgE(na,nj) ! remaining energy after ionization
                 IF(E1 > zero) THEN ! ionization occurs
                    ni = FIND_BIN(elctreV,elctDeV,E1)
                    IF((elctreV(ni) > E1) .and. (ni /= 1)) THEN
                       E4 = elctreV(ni)-half*elctDeV(ni)
                       ni = ni-1
                    ELSE
                       E4 = elctreV(ni)+half*elctDeV(ni)
                    END IF
                    Eion = elctreV(ne) - E4 ! "discretized" ionization energy
                    cs_ion = sum_cs_rees_ion(ne,na,j)
                    cs_sec = sum_cs_rees_sec(ne,na,j)
                    IF(ni > 0) THEN  ! will this criterion ever fail?
                       DO nm = 1, ni
                          IF(nm == ni) THEN
                             E3 = enrgE(na,nj)/Eion
                          ELSE
                             E3 = one
                          END IF
                          sigma_dif_ion = DIFCS_SEC(elctreV(ne)-elctreV(nm)-Eion,elctreV(ne),Eion)*eCS(ne,na,nj)/cs_ion ! divide to ensure normalized
                          sigma_dif_sec = DIFCS_SEC(elctreV(nm),elctreV(ne),Eion)*eCS(ne,na,nj)/cs_sec
                          Ael_ion(nm,ne) = Ael_ion(nm,ne) + den(nl,isp(na))*elctDeV(ne)*(sigma_dif_ion)*E3 
                          Ael_sec(nm,ne) = Ael_sec(nm,ne) + den(nl,isp(na))*elctDeV(ne)*(sigma_dif_sec)     
                       END DO
                    END IF
                    sum_cs_ion(ne,na) = sum_cs_ion(ne,na) + eCS(ne,na,nj) 
                 END IF
              END DO

           END DO
        END DO

        Ael = Ael_ext + Ael_ion + Ael_sec

        !
        !   .. Calculate Electron flux 
        !

        eFLUX(nl,nelb) = Selz(nelb)/(den(nl,iCO2)*(sum_cs_ext(nelb,1)+sum_cs_ion(nelb,1)) &
             +den(nl,iN2)*(sum_cs_ext(nelb,2)+sum_cs_ion(nelb,2)))
        DO n = nelb-1, 1, -1 ! eqn 20 in photoelectron notes (with dE(j) =dE(j+1))
           sum1 = DOT_PRODUCT(Ael(n,n+1:nelb),eFLUX(nl,n+1:nelb))
           sum2 = DOT_PRODUCT(Ael(n+1,n+2:nelb),eFLUX(nl,n+2:nelb))
           sum_tot = half * (sum1 + sum2)
           rn = eden*LOSS(elctreV(n),eden,te(nl))
           rnp1 = eden*LOSS(elctreV(n+1),eden,te(nl))
           Ann = (rn / (elctreV(n+1)-elctreV(n))) + half*( den(nl,iCO2)* (sum_cs_ext(n,1) + sum_cs_ion(n,1) ) &
                + den(nl,iN2)* (sum_cs_ion(n,2)   + sum_cs_ext(n,2)) )
           An1n1 = (rnp1 / (elctreV(n+1)-elctreV(n))) - half*(den(nl,iCO2)*(sum_cs_ext(n+1,1) + sum_cs_ion(n+1,1)) & 
                + den(nl,iN2)*(sum_cs_ion(n+1,2) + sum_cs_ext(n+1,2)))
           IF (Ann < 1.0E-20) THEN ! if Ann = 0
              eFLUX(nl,n) = 0.57 * eden
           ELSE
              eFLUX(nl,n) = ( half*(Selz(n)+Selz(n+1)) + An1n1*eFLUX(nl,n+1) + sum_tot) / Ann
           END IF
        END DO

     END IF
  
  END DO 
   
  ! Secondary production rates

  DO nl = nibot, nlev
     pSz = zero
     do i = 1, nabs_el_thk
        do n = 1, nelb  !-1
           pSz = pSz + den(nl,isp(i))*eCS_ion(n,i)*eFLUX(nl,n)*elctDeV(n) 
        enddo
     enddo
     pS(nl) = pSz
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
  
  DO n = 1, nabs_el_thn
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

  FUNCTION LOSS(E,NE,TE) ! loss of suprathermal electrons with energy E into thermal with density NE and temperature TE
    USE PRECISION
    USE CONSTANTS
    IMPLICIT NONE
    REAL(RP) :: LOSS, E, NE, TE, Ee

    Ee = 8.618E-5_RP*TE
    IF (Ee > E) THEN
	LOSS = zero
    ELSE 
	LOSS = (3.37e-12_RP/((E**0.94_RP)*(NE**0.03_RP)))*((E-Ee)/(E-0.53_RP*Ee))**2.36_RP
    END IF

	END FUNCTION LOSS

END SUBROUTINE eldep1
