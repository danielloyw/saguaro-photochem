MODULE CONSTANTS

  USE PRECISION

  !                                                                  
  ! .. ordinary numbers                                             
  !                                                                  

  REAL(RP), PARAMETER :: zero  = 0.0
  REAL(RP), PARAMETER :: half  = 0.5
  REAL(RP), PARAMETER :: one   = 1.0
  REAL(RP), PARAMETER :: two   = 2.0
  REAL(RP), PARAMETER :: three = 3.0
  REAL(RP), PARAMETER :: four  = 4.0
  REAL(RP), PARAMETER :: five  = 5.0
  REAL(RP), PARAMETER :: six   = 6.0
  REAL(RP), PARAMETER :: seven = 7.0
  REAL(RP), PARAMETER :: eight = 8.0
  REAL(RP), PARAMETER :: nine  = 9.0
  REAL(RP), PARAMETER :: ten   = 10.
  
  !                                                                  
  ! .. geometrical and physical constants in cgs units               
  !                                                                  
  
  REAL(RP), PARAMETER :: pi=3.1415926536       ! pi   
  REAL(RP), PARAMETER :: sqpi=1.7724538509     ! SQRT(pi)                       
  REAL(RP), PARAMETER :: fourpi=12.56637061    ! 4 * pi                        
  REAL(RP), PARAMETER :: amu=1.660539E-24      ! atomic mass unit            
  REAL(RP), PARAMETER :: rkb=1.38065E-16       ! Boltzmann's constant 
  REAL(RP), PARAMETER :: rgas = rkb/amu        ! gas constant        
  REAL(RP), PARAMETER :: crs0=2.654E-2         ! atomic cross section         
  REAL(RP), PARAMETER :: cs=2.997925E10        ! speed of light               
  REAL(RP), PARAMETER :: hplnk=6.626070E-27    ! Planck's constant            
  REAL(RP), PARAMETER :: Gcnst=6.6743E-8       ! Gravitational constant      
  REAL(RP), PARAMETER :: crad2=1.43878         ! hplnk*cs/rkb                  
  REAL(RP), PARAMETER :: hc=1.986446E-16       ! hplnk*cs                     
  REAL(RP), PARAMETER :: hck=1.43883           ! crad2                        
  REAL(RP), PARAMETER :: AU=1.495979E13        ! Astronomical unit            
  REAL(RP), PARAMETER :: nlosch=2.687E19       ! Loschmidt's number   
  
  !
  !  .. Conversion factors
  !
  
  REAL(RP), PARAMETER :: cm_to_km=1.E-5        
  REAL(RP), PARAMETER :: km_to_cm=1.E+5         
  REAL(RP), PARAMETER :: cm_to_aa=1.E+8
  REAL(RP), PARAMETER :: aa_to_cm=1.E-8
  REAL(RP), PARAMETER :: ev_to_egr=1.609E-12
  REAL(RP), PARAMETER :: erg_to_ev=6.215E+11
  REAL(RP), PARAMETER :: rad_to_deg=180./pi
  REAL(RP), PARAMETER :: deg_to_rad = pi/180.
  
END MODULE CONSTANTS




