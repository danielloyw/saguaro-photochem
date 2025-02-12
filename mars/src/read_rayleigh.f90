SUBROUTINE READ_RAYLEIGH

  !  .. Modules

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : INTRP

  !  .. Declarations

  INTEGER :: n, nw, nl, nrtab
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: wrtab, crtab

  !  .. Rayleigh Scattering Optical Depth

  OPEN(unit=64,file='../data/photons/CO2_Rayleigh.dat',status='old',action='read')
     READ(64,*) nrtab ! number of entries
     ALLOCATE(wrtab(nrtab),crtab(nrtab))
     DO n = 1, nrtab
        READ(64,*) wrtab(n), crtab(n) ! wavelength, scattering cross section
     END DO
  CLOSE(unit=64)
     
  ALLOCATE(crs_ray(ncrsC),tau_ray(nlev,ncrsC))
  CALL INTRP(wrtab,crtab,wcrsC,crs_ray)

  DO nw = 1, ncrsC
     DO nl = 1, nlev
        IF(wcrsC(nw) > 500._RP) THEN
           tau_ray(nl,nw) = crs_ray(nw)*col(nl)
        ELSE
           tau_ray(nl,nw) = zero
        END IF
     END DO
  END DO

  DEALLOCATE(wrtab,crtab)

  RETURN
END SUBROUTINE READ_RAYLEIGH
