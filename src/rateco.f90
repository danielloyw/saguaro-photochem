 SUBROUTINE RATECO

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : LOCATE

  IMPLICIT NONE

  !  .. Internal Variables

  REAL(RP) :: rk0, rk1, rk2
  REAL(RP) :: FC, PFL, FCL, NF, CF, FF, XF
  REAL(RP) :: rtmp, rprs
  INTEGER :: nr, nl, nm, nt1, nt2, np1, np2, nt


  ALLOCATE(rt(nrct,nlev),rct(nrct,nlev)) 

  
  !  .. Calculate neutral reaction rates at temperature
  
  DO nr = 1, nrct
     IF ((itype(nr) < -4) .or. (itype(nr) > 7)) THEN

        WRITE(*,"(5X,' ERROR, INCORRECT REACTION TYPE, nr = ',I4,' itype = ',I3,' ABORTING ...')") nr, itype(nr)
        STOP

     ELSE IF (itype(nr) == 1) THEN                                          ! Unimolecular

        irct(2,nr) = 0
        DO nl = 1, nlev
           rt(nr,nl) = rck(1,nr)
        END DO

     ELSE IF(itype(nr) == 2) THEN                                           ! Bimolecular

        DO nl = 1, nlev
           rt(nr,nl) = rck(1,nr)*(tn(nl)**rck(2,nr))*EXP(rck(3,nr)/tn(nl))
        END DO

     ELSE IF (itype(nr) == 3) THEN                                          ! Association

        DO nl = 1, nlev
           rk1 = rck(1,nr)*(tn(nl)**rck(2,nr))*EXP(rck(3,nr)/tn(nl))
           rk0 = rck(4,nr)*(tn(nl)**rck(5,nr))*EXP(rck(6,nr)/tn(nl))
           IF(rck(10,nr)== zero) THEN
              rt(nr,nl) = rk0*rk1*den(nl,0)/(rk0*den(nl,0)+rk1)
           ELSE
              FC = rck(10,nr)
              PFL = LOG10(rk0*den(nl,0)/rk1)
              FCL = LOG10(FC)
              NF = 0.75_RP - 1.27_RP*FCL
              CF = -0.4_RP-0.67_RP*FCL
              FF =10._RP**(FCL/(one+((PFL+CF)/(NF-0.14_RP*(PFL+CF)))**2))
!              rt(nr,nl) = FF*(rk0*rk1*den(nl,0)/(rk0*den(nl,0)+rk1))
              rt(nr,nl) = FF*(rk0*rk1/(rk0*den(nl,0)+rk1))
           END IF
        END DO

     ELSE IF (itype(nr) == 4) THEN                                           ! Association & Radiative Association

        DO nl = 1, nlev
           rk1 = rck(1,nr)*(tn(nl)**rck(2,nr))*EXP(rck(3,nr)/tn(nl))
           rk0 = rck(4,nr)*(tn(nl)**rck(5,nr))*EXP(rck(6,nr)/tn(nl))
           rk2 = rck(7,nr)*(tn(nl)**rck(8,nr))*EXP(rck(9,nr)/tn(nl))
           IF(rck(10,nr)== zero) THEN
              rt(nr,nl) = rk2 + rk0*rk1*den(nl,0)/(rk0*den(nl,0)+rk1)
           ELSE
              FC = rck(10,nr)
              PFL = LOG10(rk0*den(nl,0)/rk1)
              FCL = LOG10(FC)
              NF = 0.75_RP - 1.27_RP*FCL
              CF = -0.4_RP-0.67_RP*FCL
              FF =10._RP**(FCL/(one+((PFL+CF)/(NF-0.14_RP*(PFL+CF)))**2))
              XF = FF/(one-FF)
!              rt(nr,nl) = (rk0*den(nl,0)*XF+rk2)*rk1/(rk0*den(nl,0)*XF+rk1)
              rt(nr,nl) = MIN(rk1, rk2 + FF*(rk0*rk1*den(nl,0)/(rk0*den(nl,0)+rk1)))

!              rk1p = rk1 - rk2
!              rt(nr,nl) = rk2 + FF*(rk0*rk1p*den(nl,0)/(rk0*den(nl,0)+rk1p)))


           END IF
        END DO

     ELSE IF (itype(nr) == 5) THEN                                           ! special for OH + CO 

        DO nl = 1, nlev
           rk1 = rck(1,nr)*(tn(nl)**rck(2,nr))*EXP(rck(3,nr)/tn(nl))
           rk0 = rck(4,nr)*(tn(nl)**rck(5,nr))*EXP(rck(6,nr)/tn(nl))
           IF(rck(10,nr)== zero) THEN
              rt(nr,nl) = rk0*rk1/(rk0*den(nl,0)+rk1)
           ELSE
              FC = rck(10,nr)
              PFL = LOG10(rk0*den(nl,0)/rk1)
              FCL = LOG10(FC)
              NF = 0.75_RP - 1.27_RP*FCL
              CF = -0.4_RP-0.67_RP*FCL
              FF =10._RP**(FCL/(one+((PFL+CF)/(NF-0.14_RP*(PFL+CF)))**2))
              rt(nr,nl) = FF*(rk0*rk1/(rk0*den(nl,0)+rk1))
           END IF
        END DO

     ELSE IF(itype(nr) == 6) THEN                                                 !   Tabulated reactions

        nt = ntab(nr)
        DO nl = 1, nlev
           IF(tn(nl) < tmp_rct(1,nt)) THEN
              nt1 = 1
              nt2 = 2
              rtmp = zero
           ELSE IF(tn(nl) > tmp_rct(ntmp_rct(nt),nt)) THEN
              nt1 = ntmp_rct(nt)-1
              nt2 = nt1 + 1
              rtmp = one
           ELSE
              nt1 = LOCATE(tmp_rct(1:ntmp_rct(nt),nt),tn(nl))
              nt2 = nt1 + 1
              rtmp = (tn(nl)-tmp_rct(nt1,nt))/(tmp_rct(nt2,nt)-tmp_rct(nt1,nt))
           END IF
           IF (LOG(prs(nl)) < plog_rct(1,nt)) THEN
              rprs = zero
              np1 = 1
              np2 = 2
           ELSE IF (LOG(prs(nl)) > plog_rct(nprs_rct(nt),nt)) THEN
              rprs = one
              np1 = nprs_rct(nt)-1
              np2 = np1 + 1
           ELSE
              np1 = LOCATE(plog_rct(1:nprs_rct(nt),nt),LOG(prs(nl)))
              np2 = np1 + 1
              rprs = (LOG(prs(nl))-plog_rct(np1,nt))/(plog_rct(np2,nt)-plog_rct(np1,nt))
           END IF
           rt(nr,nl) = (one-rtmp)*(one-rprs)*rct_tab(np1,nt1,nt)                      &
                     + rtmp*(one-rprs)*rct_tab(np2,nt1,nt)                            &
                     + (one-rtmp)*rprs*rct_tab(np1,nt2,nt)                            &
                     + rtmp*rprs*rct_tab(np2,nt2,nt) 
           rt(nr,nl) = EXP(rt(nr,nl))
        END DO

     ELSE IF (itype(nr) == 7) THEN                                                ! Heterogenous

        nm = irct(1,nr)
        DO nl = 1, nlev
           rt(nr,nl) = rck(1,nr)*SQRT(rgas*tn(nl)/two/pi/mmw(nm))*surfarea(nl)/den(nl,0)
        END DO

     ELSE IF (itype(nr) == -1) THEN                  ! Unimolecular Reaction

        DO nl = 1, nlev
           rt(nr,nl)=rck(1,nr)*(300._RP/tn(nl))**rck(2,nr)                                          &
             *two*EXP(rck(3,nr)/tn(nl))/(one+EXP(rck(4,nr)/tn(nl)))
        END DO

     ELSE IF (itype(nr) == -2) THEN                  ! normal two-body reaction

        DO nl = 1, nlev
           rt(nr,nl)=rck(1,nr)*(300._RP/tn(nl))**rck(2,nr)                                          &
                *two*EXP(rck(3,nr)/tn(nl))/(one+EXP(rck(4,nr)/tn(nl)))
        END DO

     ELSE IF (itype(nr) == -3) THEN                  ! 3-body reaction

        DO nl = 1, nlev
           rt(nr,nl)=rck(1,nr)*(300._RP/tn(nl))**rck(2,nr)                                          &
                *two*EXP(rck(3,nr)/tn(nl))/(one+EXP(rck(4,nr)/tn(nl))) 
        END DO

     ELSE IF (itype(nr) == -4) THEN                  ! electron recombination

        DO nl = 1, nlev
           rt(nr,nl)=rck(1,nr)*(300._RP/te(nl))**rck(2,nr)                                          &
                *two*EXP(rck(3,nr)/te(nl))/(one+EXP(rck(4,nr)/te(nl)))
        END DO

     END IF

  END DO

  RETURN
END SUBROUTINE RATECO
