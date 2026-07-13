module diffusion

contains 

subroutine diffusion_setup
! This subroutine calculates diffusion coefficients for the continuity
! equation, and sets up the boundary conditions. Note that there is no ion
! diffusion. Please refer to the pchem_models.pdf notes. 
  
  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use types, only: wp => dp
  use constants
  use global_variables

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  ! polarizability of N2
  real(wp), parameter :: polN2 = 17.6E-25
  ! polarizability of CO2
  real(wp), parameter :: polCO2 = 2.63E-24
  ! half-step quantities
  real(wp), allocatable, dimension(:) :: Tn_mid, grv_mid, eK_mid, mass_mid 
  real(wp), allocatable, dimension(:,:) :: Ht_mid, Df_mid
  ! dTn/dr
  real(wp), allocatable, dimension(:) :: dTndr
  ! hard sphere collisional area
  real(wp) :: sigma
  ! reduced mass
  real(wp) :: mu

  ! loop variables
  integer :: i_sp1, i_sp2, i_z

  !----------------------------------------------------------------------------
  !  Calculate neutral molecular diffusion coefficients
  !----------------------------------------------------------------------------
  
  allocate(Df(0:n_z,n_sp))
  Df = zero
  
  do i_sp1 = 1, n_neu
    do i_z = 1, n_z
      ! hard sphere approximation (Chapman-Enskog)
      if (dtype(i_sp1) == 0) then                                                       
        sigma = 3.0e-15_wp ! pi*r^2 with r ~ 3 A
        mu = amu / (one/mmw(i_sp1) + one/mmw(iCO2))
        Df(i_z, i_sp1) = three / eight * sqrt_pi / sigma &
          * sqrt(kB * Tn(i_z) / mu / two) / den(i_z,0)
      
      ! Mason+Morrero [1970] formula 136
      else if (dtype(i_sp1) == 1) then                                                  
        Df(i_z, i_sp1) = Ad(i_sp1) * Tn(i_z)**sd1(i_sp1) & 
          * exp(-sd2(i_sp1) / Tn(i_z)) / prs(i_z)
      
      ! Mason+Morrero [1970] formula 135
      else if (dtype(i_sp1) == 2) then
        Df(i_z, i_sp1) = Ad(i_sp1) * Tn(i_z)**sd1(i_sp1) &
          * exp(-sd2(i_sp1) / Tn(i_z)) * exp(-sd3(i_sp1) / Tn(i_z)**two) &
          * log(phi(i_sp1) / kB / Tn(i_z))**two / prs(i_z)
      
      else
        write(*,'("Error: Invalid diffusion setting for ", A, & 
          "! Exiting ... ")') trim(adjustl(sp_list(i_sp1)))
        stop
      end if
    end do
  end do

  !----------------------------------------------------------------------------
  !  Calculate half-step variables
  !----------------------------------------------------------------------------
  
  ! extrapolate for ghost cell below
  do concurrent(i_sp1 = 1:n_neu)
     Df(0,i_sp1) = Df(1,i_sp1) * Df(1,i_sp1) / Df(2,i_sp1)
  end do
  
  allocate(Tn_mid(0:n_z), grv_mid(0:n_z), eK_mid(0:n_z), mass_mid(0:n_z))
  allocate(Ht_mid(0:n_z,0:n_sp), Df_mid(0:n_z,n_sp))
  allocate(dTndr(0:n_z))

  do i_z = 1, n_z-1
    Tn_mid(i_z) = half * (Tn(i_z+1) + Tn(i_z))
    dTndr(i_z) = (Tn(i_z+1) - Tn(i_z)) / dr_mid(i_z)
    grv_mid(i_z) = half * (grv(i_z+1) + grv(i_z))
    mass_mid = half * (mass(i_z+1) + mass(i_z))
    eK_mid(i_z) = half * (eK(i_z+1) + eK(i_z))
    Ht_mid(i_z,0) = kB * Tn_mid(i_z) / grv_mid(i_z) / mass_mid(i_z) / amu
    do i_sp1 = 1, n_diff
      i_sp2 = im_diff_all(i_sp1)
      Ht_mid(i_z, i_sp2) = kB * Tn_mid(i_z) / grv_mid(i_z) / mmw(i_sp2) / amu
      Df_mid(i_z, i_sp2) = half * (Df(i_z+1, i_sp2) + Df(i_z, i_sp2))
    end do
  end do
  
  ! extrapolate for ghost cell below
  i_z = 0
  Tn_mid(i_z) = Tn_mid(1) * Tn_mid(1) / Tn_mid(2)
  dTndr(i_z) = dTndr(1)
  grv_mid(i_z) = grv_mid(1) * grv_mid(1) / grv_mid(2)
  mass_mid(i_z) = mass_mid(1) * mass_mid(1) / mass_mid(2)
  eK_mid(i_z) = eK_mid(1) * eK_mid(1) / eK_mid(2)
  Ht_mid(i_z,0) = kB * Tn_mid(i_z) / grv_mid(i_z) / mass_mid(i_z) / amu
  do i_sp1 = 1, n_diff
    i_sp2 = im_diff_all(i_sp1)
    Ht_mid(i_z,i_sp2) = kB * Tn_mid(i_z) / grv_mid(i_z) / mmw(i_sp2) / amu
    Df_mid(i_z,i_sp2) = half * (Df(i_z+1,i_sp2) + Df(i_z,i_sp2))
  end do

  ! extrapolate for ghost cell above
  i_z = n_z
  Tn_mid(i_z) = Tn_mid(n_z-1) * Tn_mid(n_z-1) / Tn_mid(n_z-2)
  dTndr(i_z) = dTndr(n_z-1) * dTndr(n_z-1) / dTndr(n_z-2)
  grv_mid(i_z) = grv_mid(n_z-1) * grv_mid(n_z-1) / grv_mid(n_z-2)
  mass_mid(i_z) = mass_mid(n_z-1) * mass_mid(n_z-1) / mass_mid(n_z-2)
  eK_mid(i_z) = eK_mid(n_z-1) * eK_mid(n_z-1) / eK_mid(n_z-2)
  Ht_mid(i_z,0) = kB * Tn_mid(i_z) / grv_mid(i_z) / mass_mid(i_z) / amu
  do i_sp1 = 1, n_diff
    i_sp2 = im_diff_all(i_sp1)
    Ht_mid(i_z,i_sp2) = kB * Tn_mid(i_z) / grv_mid(i_z) / mmw(i_sp2) / amu
    Df_mid(i_z,i_sp2) = Df(i_z,i_sp2)
  end do
  
  !? why?
  ! do i_sp1 = 1, n_diff
    ! i_sp2 = im_diff_all(i_sp1)
    ! Lower Boundary
    ! if (ibnd(i_sp2,1) == 2) then
    ! specified velocity; negative upwards
    ! bval(i_sp2,1) = max(bval(i_sp2,1), &
        ! -eK(1) / (kB * Tn(1) / grv(1) / mmw(i_sp2) / amu))
    ! end if
  ! end do

  !----------------------------------------------------------------------------
  !  Implement fluxes
  !----------------------------------------------------------------------------

  ! Please refer to pchem_models.pdf notes. 

  allocate(alpha(0:n_z,n_diff), beta(0:n_z,n_diff))
  allocate(a(n_z,n_diff), b(n_z,n_diff), c(n_z,n_diff))
  alpha = zero
  beta = zero
  a = zero
  b = zero
  c = zero

  do i_sp1 = 1, n_diff
    i_sp2 = im_diff_all(i_sp1)
    do concurrent (i_z = 0:n_z)
      ! equations 7 from notes
      alpha(i_z,i_sp1) = (Df_mid(i_z,i_sp2) + eK_mid(i_z)) / dr_mid(i_z)
      beta(i_z,i_sp1) = half * (Df_mid(i_z,i_sp2) / Ht_mid(i_z,i_sp2) &
        + eK_mid(i_z) / Ht_mid(i_z,0) &
        + (Df_mid(i_z,i_sp2) + eK_mid(i_z)) * dTndr(i_z) / Tn_mid(i_z))
    end do
  end do

  do i_sp1 = 1, n_diff
    i_sp2 = im_diff_all(i_sp1)
     
    ! Interior
    ! equations 10 from notes
    do concurrent (i_z = 2:n_z-1)
      a(i_z,i_sp1) = -(alpha(i_z-1,i_sp1) - beta(i_z-1,i_sp1)) &
        * (rz_mid(i_z-1) / rz(i_z))**two / dr(i_z)
      b(i_z,i_sp1) = (alpha(i_z-1,i_sp1) + beta(i_z-1,i_sp1)) &
        * (rz_mid(i_z-1) / rz(i_z))**two / dr(i_z) &
        + ((alpha(i_z,i_sp1) - beta(i_z,i_sp1)) &
        * (rz_mid(i_z) / rz(i_z))**two / dr(i_z))
      c(i_z,i_sp1) = -(alpha(i_z,i_sp1) + beta(i_z,i_sp1)) &
        * (rz_mid(i_z) / rz(i_z))**two / dr(i_z)
    end do
    
    ! Lower Boundary ?fix with updated equations
    i_z = 1
    a(i_z,i_sp1) = zero    ! since there is nothing lower
    if (ibnd(i_sp2,1) == 1) then
    ! diffusive equilibrium (so div(flux) = 0)
      b(i_z,i_sp1) = zero
      c(i_z,i_sp1) = zero
    else if (ibnd(i_sp2,1) == 2) then
    ! fixed velocity; positive upwards
      b(i_z,i_sp1) = -two * bval(i_sp2,1) / dr(i_z) &
        + two * (alpha(i_z,i_sp1) - beta(i_z,i_sp1)) &
        * (rz_mid(i_z) / rz(i_z))**two / dr(i_z)
      c(i_z,i_sp1) = -two * (alpha(i_z,i_sp1) + beta(i_z,i_sp1)) &
        * (rz_mid(i_z) / rz(i_z))**two / dr(i_z)
    else if (ibnd(i_sp2,1) == 3) then
    ! fixed mole fraction (no contribution from above)
      b(i_z,i_sp1) = zero
      c(i_z,i_sp1) = zero
    else if (ibnd(i_sp2,1) == 4) then
    ! fixed flux
      b(i_z,i_sp1) = (alpha(i_z,i_sp1) - beta(i_z,i_sp1)) &
        * (rz_mid(i_z-1) / rz(i_z))**two / dr_mid(i_z)
      c(i_z,i_sp1) = -(alpha(i_z,i_sp1) + beta(i_z,i_sp1)) &
        * (rz_mid(i_z-1) / rz(i_z))**two / dr_mid(i_z)
    else
      write(*,'("Error: Invalid lower boundary condition for ", A, & 
        "! Exiting ... ")') trim(adjustl(sp_list(i_sp1)))
      stop
    end if

    ! Upper Boundary ?fix with updated equations
    i_z = n_z
    c(i_z,i_sp1) = zero    ! since there is nothing higher
    if (ibnd(i_sp2,2) == 1) then
    ! diffusive equilibrium (so div(flux) = 0)
      a(i_z,i_sp1) = zero
      b(i_z,i_sp1) = zero
    else if ((ibnd(i_sp2,2) == 2) .or. (ibnd(i_sp2,2) == 3)) then
    ! Jeans or fixed velocity; positive upwards
      if (ibnd(i_sp2,2) == 2) then
      ! Jeans velocity, bval(i_sp2,2) = 1: escaping at Jeans rate
        ! replace bval with actual velocity
        bval(i_sp2,2) = bval(i_sp2,2) * wJeans(mmw(i_sp2), rz(n_z), Tn(n_z))
      end if
      a(i_z,i_sp1) = -two * (rz_mid(i_z-1) / rz(i_z))**two &
        * (alpha(i_z-1,i_sp1) - beta(i_z-1,i_sp1)) / dr(i_z)
      b(i_z,i_sp1) = two * bval(i_sp2,2) / dr(i_z) &
        + two * (rz_mid(i_z-1) / rz(i_z))**two &
        * (alpha(i_z-1,i_sp1) + beta(i_z-1,i_sp1)) / dr(i_z)
    else if (ibnd(i_sp2,2) == 4) then
    ! fixed flux
      a(i_z,i_sp1) = -(alpha(i_z-1,i_sp1) - beta(i_z-1,i_sp1)) &
        * (rz_mid(i_z-2) / rz(i_z-1))**2 / dr_mid(i_z)
      b(i_z,i_sp1) = (alpha(i_z-1,i_sp1) + beta(i_z-1,i_sp1)) &
        * (rz_mid(i_z-2) / rz(i_z-1))**2 / dr_mid(i_z)
    else
      write(*,'("Error: Invalid upper boundary condition for ", A, &
        "! Exiting ... ")') trim(adjustl(sp_list(i_sp1)))
      stop
    end if
  end do

  deallocate(Tn_mid, grv_mid, eK_mid, mass_mid, Ht_mid, Df_mid, dTndr)

end subroutine diffusion_setup


real elemental function wJeans(m, r, T)
! This function calculates the Jeans escape speed, given the particle mass "m",
! radius from the center of the planet "r", and the temperature "T". 
! Lambda is the Jeans parameter. 
  use types, only: wp => dp
  use constants
  use global_variables
  implicit none
  real(wp), intent(in) :: m, r, T
  real(wp) :: lambda
  
  lambda = GM * m * amu / kB / T / r
  wJeans = sqrt(kB * T / two / pi / m / amu) * (one + lambda) * exp(-lambda)
end function wJeans

end module
