SUBROUTINE READ_AEROSOL

  !  .. Modules

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : LOCATE

  !  .. Declarations

!  INTEGER :: nw, na
!  REAL(RP), ALLOCATABLE, DIMENSION(:) :: zaer, wvaer
!  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: kaer_tab
!  REAL(RP) :: dum

  ALLOCATE(tau_aer(nlev,ncrsC))
  ALLOCATE(surfarea(nlev))

  tau_aer = zero
  surfarea = zero

  RETURN
END SUBROUTINE READ_AEROSOL
