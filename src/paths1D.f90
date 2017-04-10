SUBROUTINE PATHS1D

  USE PRECISION
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : LOCATE
  IMPLICIT none
  
  !  .. Local Variables
  
  REAL(RP) :: sin_sza
  REAL(RP) :: rtan
  REAL(RP) :: rtn2
  INTEGER :: nz, mz, ntn

  !  .. Parameters
  
  REAL(RP), PARAMETER :: zero=0., one=1., two=2.
  
  !  .. Initialization

  ALLOCATE(ntan(nlev),ds(nlev,nlev))                         ! move this to inside PATHS1
  ds = 1.0E+33_RP

     sin_sza = SQRT(1.-cos_sza*cos_sza)
     IF (cos_sza >= zero) THEN                                ! Day
        illum = 1
        DO nz = 1, nlev-1
           rtan = rz(nz) * sin_sza ! tangent radius of light ray
           rtn2 = rtan * rtan
           ntan(nz) = nz
           DO mz = nz, nlev-1
              ds(mz,nz)=SQRT(rz(mz+1)*rz(mz+1)-rtn2)-SQRT(rz(mz)*rz(mz)-rtn2)
           END DO
        END DO
     ELSE IF((cos_sza < zero) .and. (rz(nlev)*sin_sza > RSHADOW)) THEN   ! Twilight
        illum = 0
        nbot = LOCATE(rz*sin_sza,RSHADOW)+1
!        WRITE(*,*) ' PATHS: nt = ',nt,' nbot = ',nbot(nt)
        DO nz = nbot, nlev-1              
           rtan = rz(nz) * sin_sza
           ntn = LOCATE(rz, rtan) + 1
           ntan(nz) = ntn 
           rtn2 = rtan * rtan	
           ds(ntn,nz) = SQRT(rz(ntn+1)*rz(ntn+1) - rtn2) 
           DO mz = ntn+1, nlev-1   
              ds(mz,nz)=SQRT(rz(mz+1)*rz(mz+1)-rtn2)-SQRT(rz(mz)*rz(mz)-rtn2)
           END DO
        END DO
        rtan = rz(nlev) * sin_sza
        ntan(nlev) = LOCATE(rz, rtan) + 1
     ELSE IF((cos_sza < zero).and.(rz(nlev)*sin_sza <= RSHADOW)) THEN  ! Night
        illum = -1
!        WRITE(*,*) ' IN SHADOW'
        DO nz = 1, nlev
           DO mz = nz, nlev
              ds(mz,nz) = 1.0E+33_RP
           END DO
        END DO
     END IF
  
  RETURN
END SUBROUTINE PATHS1D
