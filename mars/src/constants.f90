! Defines various constants that are used in program
module constants
  use types, only: wp => dp

  ! Ordinary numbers, so that they have specified precision
  real(wp), parameter :: zero  = 0.0
  real(wp), parameter :: half  = 0.5
  real(wp), parameter :: one   = 1.0
  real(wp), parameter :: two   = 2.0
  real(wp), parameter :: three = 3.0
  real(wp), parameter :: four  = 4.0
  real(wp), parameter :: five  = 5.0
  real(wp), parameter :: six   = 6.0
  real(wp), parameter :: seven = 7.0
  real(wp), parameter :: eight = 8.0
  real(wp), parameter :: nine  = 9.0
  real(wp), parameter :: ten   = 10.0

  ! Geometrical and physical constants in cgs units
  real(wp), parameter :: pi      = 3.14159265359      ! pi
  real(wp), parameter :: sqrt_pi = 1.77245385091      ! SQRT(pi)
  real(wp), parameter :: fourpi  = 12.5663706144      ! 4 * pi
  real(wp), parameter :: amu     = 1.66053906892E-24  ! atomic mass unit
  real(wp), parameter :: kB      = 1.380649E-16       ! Boltzmann's constant
  real(wp), parameter :: Rgas    = kB/amu             ! gas constant
  real(wp), parameter :: crs0    = 2.654E-2           ! atomic cross section
  real(wp), parameter :: cs      = 2.99792458E10      ! speed of light
  real(wp), parameter :: Gcnst   = 6.67430E-8         ! gravitational constant
!  real(wp), parameter :: hplanck = 6.62607015E-27     ! Planck's constant
!  real(wp), parameter :: crad2   = 1.43878            ! hplanck*cs/kB
!  real(wp), parameter :: hc      = 1.986446E-16       ! hplanck*cs
  real(wp), parameter :: AU      = 1.495978707E13     ! Astronomical unit
  real(wp), parameter :: nlosch  = 2.686780111E19     ! Loschmidt's number

  ! Conversion factors
  real(wp), parameter :: cm_to_km   = 1.E-5
  real(wp), parameter :: km_to_cm   = 1.E+5
  real(wp), parameter :: cm_to_aa   = 1.E+8
  real(wp), parameter :: aa_to_cm   = 1.E-8
  real(wp), parameter :: ev_to_erg  = 1.602176634E-12
  real(wp), parameter :: erg_to_ev  = 6.24150907E+11
  real(wp), parameter :: rad_to_deg = 180.0_wp/pi
  real(wp), parameter :: deg_to_rad = pi/180.0_wp

end module constants