              DO j = 1, ipath(na,1)+ipath(na,2)
                 E1 = elctreV(m)-enrgE(na,j)
                 nm = FIND_BIN(elctreV,elctDeV,E1)
                 IF((E1>zero) .and. (nm<m)) THEN
                    Ael_ext(nm,m) = Ael_ext(nm,m) + den(i,isp(na))*eCS(m,na,j)*(enrgE(na,j)/(elctreV(m)-elctreV(nm)))
                 ELSE IF((E1>zero) .and. (nm==m)) THEN
                    Ael_ext(nm-1,m) = Ael_ext(nm-1,m) + den(i,isp(na))*eCS(m,na,j)*(enrgE(na,j)/(elctreV(nm)-elctreV(nm-1)))
                 ENDIF
                 sum_cs_ext(m,na) = sum_cs_ext(m,na) + eCS(m,na,j) 
              END DO
