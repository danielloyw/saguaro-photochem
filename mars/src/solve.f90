module solve

contains

subroutine solve_chem
! This subroutine uses Newton's method to calculate the densities of all
! chemistry loop species given production and loss rates from chemical, photo
! and electron reactions. Note that the returned densities are not densities at
! equilibrium, and each iteration of Newton's method is not a time step. Please
! refer to pchem_models.pdf notes.

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use types, only: wp => dp
  use constants
  use global_variables
  use lapack95, only: gesv
  use photo, only: cal_photo_rates

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  implicit none
  ! number of species in chemistry loop + electrons
  integer :: n_chem_e

  ! --- Continuity equation ---------------------------------------------------
  ! Jacobian of rate terms: unit=s-1
  ! dim=(numerator species #, denominator species #)
  real(wp), allocatable, dimension(:,:) :: rjac
  ! full Jacobian: unit=s-1; dim=(numerator species #, denominator species #)
  real(wp), allocatable, dimension(:,:) :: jac
  ! continuity condition
  real(wp), allocatable, dimension(:) :: F_con
  ! change from one estimate of densities by Newton's method to next
  real(wp), allocatable, dimension(:) :: delN

  ! --- Convergence parameters ------------------------------------------------
  integer :: iter
  real(wp), parameter :: damping = 0.5_wp
  integer, parameter :: iter_max = 10
  real(wp), parameter :: tol = 1.0e-10_wp

  ! loop variables
  integer :: i_z, i_sp, i_sp2, i_r
  ! temporary / dummy variables
  integer :: ti_sp1, ti_sp2, ti_sp3, ti_sp4, ti_sp5
  integer :: ti_sp_c1, ti_sp_c2, ti_sp_c3, ti_sp_c4, ti_sp_c5
  real(wp) :: t_freq1, t_freq2

  !----------------------------------------------------------------------------
  !  Solve for densities at each altitude
  !----------------------------------------------------------------------------

  if (recal_photo_chem) then
    call cal_photo_rates
!    call cal_electrons
  end if

  if (have_ions) then
    n_chem_e = n_chem + 1
    ! recalculate electron density for charge conservation
    do concurrent (i_z = 1:n_z)
      den(i_z,n_sp) = dot_product(chrg(1:n_sp-1), den(i_z,1:n_sp-1))
    end do
  else
    n_chem_e = n_chem
  end if

  allocate(rjac(0:n_chem_e,0:n_chem_e)) ! 0th index to catch null reactants
  allocate(jac(n_chem_e,n_chem_e), F_con(n_chem_e), delN(n_chem_e))

  do i_z = 1, n_z
    iter = 0
    F_con = ten * tol ! initialize to fail convergence
    delN = ten * tol ! initialize to fail convergence
    NEWTON: do
      iter = iter + 1

      ! --- Set up terms in continuity equation -------------------------------
      ! note that there are no cross-altitude terms since there is no diffusion
      ! only Eq. 14b in notes
      rjac = zero

      ! Photo reactions
      pr_p(i_z,:) = zero; ls_p(i_z,:) = zero; rate_prct(:,i_z) = zero
      do i_r = 1, n_prct
        ti_sp1 = im_photo_all(1,i_r); ti_sp_c1 = im_all_chem(ti_sp1)
        ti_sp2 = im_photo_all(2,i_r); ti_sp_c2 = im_all_chem(ti_sp2)
        ti_sp3 = im_photo_all(3,i_r); ti_sp_c3 = im_all_chem(ti_sp3)
        ti_sp4 = im_photo_all(4,i_r); ti_sp_c4 = im_all_chem(ti_sp4)
        ti_sp5 = im_photo_all(5,i_r); ti_sp_c5 = im_all_chem(ti_sp5)

        rate_prct(i_r,i_z) = freq_prct(i_r,i_z) * den(i_z,ti_sp1)
        ls_p(i_z,ti_sp1) = ls_p(i_z,ti_sp1) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp2) = pr_p(i_z,ti_sp2) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp3) = pr_p(i_z,ti_sp3) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp4) = pr_p(i_z,ti_sp4) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp5) = pr_p(i_z,ti_sp5) + rate_prct(i_r,i_z)

        rjac(ti_sp_c1,ti_sp_c1) = rjac(ti_sp_c1,ti_sp_c1) - freq_prct(i_r,i_z)
        rjac(ti_sp_c2,ti_sp_c1) = rjac(ti_sp_c2,ti_sp_c1) + freq_prct(i_r,i_z)
        rjac(ti_sp_c3,ti_sp_c1) = rjac(ti_sp_c3,ti_sp_c1) + freq_prct(i_r,i_z)
        rjac(ti_sp_c4,ti_sp_c1) = rjac(ti_sp_c4,ti_sp_c1) + freq_prct(i_r,i_z)
        rjac(ti_sp_c5,ti_sp_c1) = rjac(ti_sp_c5,ti_sp_c1) + freq_prct(i_r,i_z)
      end do

      ! Electron reactions
      !? no Jacobian?
      pr_e(i_z,:) = zero; ls_e(i_z,:) = zero; rate_erct(:,i_z) = zero
      do i_r = 1, n_erct
        ti_sp1 = im_e_all(1,i_r)
        ti_sp2 = im_e_all(2,i_r)
        ti_sp3 = im_e_all(3,i_r)
        ti_sp4 = im_e_all(4,i_r)
        ti_sp5 = im_e_all(5,i_r)

        rate_erct(i_r,i_z) = freq_erct(i_r,i_z) * den(i_z,ti_sp1)
        ls_e(i_z,ti_sp1) = ls_e(i_z,ti_sp1) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp2) = pr_e(i_z,ti_sp2) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp3) = pr_e(i_z,ti_sp3) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp4) = pr_e(i_z,ti_sp4) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp5) = pr_e(i_z,ti_sp5) + rate_erct(i_r,i_z)
      end do

      ! Chemical reactions
      pr_c(i_z,:) = zero; ls_c(i_z,:) = zero; rate_rct(:,i_z) = zero
      do i_r = 1, n_rct
        if (chem_type(i_r) == 1) then    ! unimolecular
          ti_sp1 = ireactant(1,i_r); ti_sp_c1 = im_all_chem(ti_sp1)
          ti_sp3 = ireactant(3,i_r); ti_sp_c3 = im_all_chem(ti_sp3)
          ti_sp4 = ireactant(4,i_r); ti_sp_c4 = im_all_chem(ti_sp4)
          ti_sp5 = ireactant(5,i_r); ti_sp_c5 = im_all_chem(ti_sp5)

          rate_rct(i_r,i_z) = rk(i_r,i_z) * den(i_z,ti_sp1)
          ls_c(i_z,ti_sp1) = ls_c(i_z,ti_sp1) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp3) = pr_c(i_z,ti_sp3) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp4) = pr_c(i_z,ti_sp4) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp5) = pr_c(i_z,ti_sp5) + rate_rct(i_r,i_z)

          rjac(ti_sp_c1,ti_sp_c1) = rjac(ti_sp_c1,ti_sp_c1) - rk(i_r,i_z)
          rjac(ti_sp_c3,ti_sp_c1) = rjac(ti_sp_c3,ti_sp_c1) + rk(i_r,i_z)
          rjac(ti_sp_c4,ti_sp_c1) = rjac(ti_sp_c4,ti_sp_c1) + rk(i_r,i_z)
          rjac(ti_sp_c5,ti_sp_c1) = rjac(ti_sp_c5,ti_sp_c1) + rk(i_r,i_z)
        else
          ti_sp1 = ireactant(1,i_r); ti_sp_c1 = im_all_chem(ti_sp1)
          ti_sp2 = ireactant(2,i_r); ti_sp_c2 = im_all_chem(ti_sp2)
          ti_sp3 = ireactant(3,i_r); ti_sp_c3 = im_all_chem(ti_sp3)
          ti_sp4 = ireactant(4,i_r); ti_sp_c4 = im_all_chem(ti_sp4)
          ti_sp5 = ireactant(5,i_r); ti_sp_c5 = im_all_chem(ti_sp5)

          rate_rct(i_r,i_z) = rk(i_r,i_z) * den(i_z,ti_sp1) * den(i_z,ti_sp2)
          ls_c(i_z,ti_sp1) = ls_c(i_z,ti_sp1) + rate_rct(i_r,i_z)
          ls_c(i_z,ti_sp2) = ls_c(i_z,ti_sp2) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp3) = pr_c(i_z,ti_sp3) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp4) = pr_c(i_z,ti_sp4) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp5) = pr_c(i_z,ti_sp5) + rate_rct(i_r,i_z)

          t_freq1 = rk(i_r,i_z) * den(i_z,ti_sp2)
          t_freq2 = rk(i_r,i_z) * den(i_z,ti_sp1)
          rjac(ti_sp_c1,ti_sp_c1) = rjac(ti_sp_c1,ti_sp_c1) - t_freq1
          rjac(ti_sp_c1,ti_sp_c2) = rjac(ti_sp_c1,ti_sp_c2) - t_freq2
          rjac(ti_sp_c2,ti_sp_c1) = rjac(ti_sp_c2,ti_sp_c1) - t_freq1
          rjac(ti_sp_c2,ti_sp_c2) = rjac(ti_sp_c2,ti_sp_c2) - t_freq2
          rjac(ti_sp_c3,ti_sp_c1) = rjac(ti_sp_c3,ti_sp_c1) + t_freq1
          rjac(ti_sp_c3,ti_sp_c2) = rjac(ti_sp_c3,ti_sp_c2) + t_freq2
          rjac(ti_sp_c4,ti_sp_c1) = rjac(ti_sp_c4,ti_sp_c1) + t_freq1
          rjac(ti_sp_c4,ti_sp_c2) = rjac(ti_sp_c4,ti_sp_c2) + t_freq2
          rjac(ti_sp_c5,ti_sp_c1) = rjac(ti_sp_c5,ti_sp_c1) + t_freq1
          rjac(ti_sp_c5,ti_sp_c2) = rjac(ti_sp_c5,ti_sp_c2) + t_freq2
        end if
      end do

      do concurrent (i_sp = 1:n_sp)
        pr(i_z,i_sp) = pr_p(i_z,i_sp) + pr_e(i_z,i_sp) + pr_c(i_z,i_sp)
        ls(i_z,i_sp) = ls_p(i_z,i_sp) + ls_e(i_z,i_sp) + ls_c(i_z,i_sp)
        dNdt(i_z,i_sp) = (den(i_z,i_sp) - den_old(i_z,i_sp)) / tstep_chem
      end do

      ! if there are no chemical species, then just calculate production and
      ! loss rates without updating densities
      if (n_chem == 0) exit

      ! --- Test convergence --------------------------------------------------
      if ((iter > iter_max) .or. (maxval(abs(F_con)) < tol) &
        .or. (maxval(abs(delN)) < tol)) then
        !write(*,'(" CHEM NEWTON ended after iteration ", I0)') iter-1
        exit
      end if

      ! --- Solve continuity equation -----------------------------------------
      ! solve F(N) = 0 using Newton's method: N' = N - F(N)/(dF/dN)
      ! jac = -dF/dN; delN = N' - N; F_con = F(N)
      ! solve (-dF/dN)(N' - N) = F(N)
      jac = zero
      F_con = zero
      do i_sp = 1, n_chem
        F_con(i_sp) = pr(i_z,im_chem_all(i_sp)) - ls(i_z,im_chem_all(i_sp)) &
          - dNdt(i_z,im_chem_all(i_sp))
        jac(i_sp,i_sp) = one / tstep_chem
        do i_sp2 = 1, n_chem
          jac(i_sp,i_sp2) = jac(i_sp,i_sp2) - rjac(i_sp,i_sp2)
        end do
      end do

      ! Jacobian for electrons
      if (have_ions) then
        F_con(n_chem_e) = (dot_product(chrg(n_neu+1:n_sp-1), &
          den(i_z,n_neu+1:n_sp-1)) - den(i_z,n_sp)) / tstep_chem
        jac(n_chem_e,n_chem_e) = one / tstep_chem
        do i_sp = 1, n_chem
          jac(n_chem_e,i_sp) = - chrg(im_chem_all(i_sp)) / tstep_chem
        end do
      end if

      delN = F_con

      call gesv(jac, delN)       ! (-dF/dN)(N' - N) = F(N)
      ! The gesv(A,B) routine solves for X in the system of linear equations
      ! A*X = B. LU decomposition with partial pivoting and row interchanges is
      ! used to factor A as A = P*L*U, where P is a permutation matrix, L is
      ! unit lower triangular, and U is upper triangular. The factored form of
      ! A is then used to solve the system of equations A*X = B.
      ! A will be overwritten by the factors L and U from the factorization of
      ! A = P*L*U; the unit diagonal elements of L are not stored. B will be
      ! overwritten by the solution matrix X.

      ! --- Update densities -------------------------------------------------
      do concurrent (i_sp = 1:n_chem)
        ! limit maximum change to ensure smooth convergence
        if (delN(i_sp) < 0) then
          delN(i_sp) = max(delN(i_sp), -damping * den(i_z,im_chem_all(i_sp)))
        else
          delN(i_sp) = min(delN(i_sp), damping * den(i_z,im_chem_all(i_sp)))
        end if
        den(i_z,im_chem_all(i_sp)) = &
          max(den(i_z,im_chem_all(i_sp)) + delN(i_sp), eps)
      end do

      ! recalculate electron density for charge conservation
      if (have_ions) then
        den(i_z,n_sp) = dot_product(chrg(1:n_sp-1), den(i_z,1:n_sp-1))
      end if

      ! update total density
      den(i_z,0) = sum(den(i_z,1:n_sp-1))

    end do NEWTON
  end do

  deallocate(rjac, jac, F_con, delN)

end subroutine solve_chem


subroutine solve_diff
! This subroutine uses Newton's method to calculate the densities of all
! diffusion loop species given diffusion, and production and loss rates from
! chemical, photo and electron reactions. Note that the returned densities are
! not densities at equilibrium, and each iteration of Newton's method is not a
! time step. Please refer to pchem_models.pdf notes.

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use types, only: wp => dp
  use constants
  use global_variables
  use photo, only: cal_photo_rates

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  implicit none

  ! --- Continuity equation ---------------------------------------------------
  ! Jacobian of rate terms: unit=s-1
  ! dim=(numerator species #, denominator species #, altitude level)
  real(wp), allocatable, dimension(:,:,:) :: rjac
  ! full Jacobian: unit=s-1
  ! dim=(numerator species #, denominator species #, altitude level)
  real(wp), allocatable, dimension(:,:,:) :: jac_a, jac_b, jac_c
  ! continuity condition
  real(wp), allocatable, dimension(:,:) :: F_con
  ! change from one estimate of densities by Newton's method to next
  real(wp), allocatable, dimension(:,:) :: delN

  ! --- Convergence parameters ------------------------------------------------
  integer :: iter
  real(wp), parameter :: damping = 0.1_wp
  integer, parameter :: iter_max = 20
  real(wp), parameter :: tol_den = 1.0e-06_wp
  real(wp), parameter :: tol_bal = 1.0e-03_wp

  ! loop variables
  integer :: i_z, i_r, i_sp1, i_sp2, i_sp3
  ! temporary / dummy variables
  integer :: ti_sp1, ti_sp2, ti_sp3, ti_sp4, ti_sp5
  integer :: ti_sp_d1, ti_sp_d2, ti_sp_d3, ti_sp_d4, ti_sp_d5
  real(wp) :: t_freq1, t_freq2

  !----------------------------------------------------------------------------
  !  Solve for densities at next time step
  !----------------------------------------------------------------------------

  allocate(rjac(0:n_diff,0:n_diff,n_z))
  allocate(jac_a(n_diff,n_diff,n_z), jac_b(n_diff,n_diff,n_z))
  allocate(jac_c(n_diff,n_diff,n_z), F_con(n_diff,n_z), delN(n_diff,n_z))

  iter = 0
  F_con = ten * tol_den ! initialize to fail convergence
  delN = ten * tol_bal ! initialize to fail convergence
  
  NEWTON: do
    iter = iter + 1

    !? better here or outside loop like in CHEM?
    if (recal_photo_diff) then
       call cal_photo_rates
!       call cal_electrons
    end if

    ! --- Set up terms in continuity equation -------------------------------
    rjac = zero

    ! Photo reactions
    pr_p = zero; ls_p = zero; rate_prct = zero
    do i_r = 1, n_prct
      ti_sp1 = im_photo_all(1,i_r); ti_sp_d1 = im_all_diff(ti_sp1)
      ti_sp2 = im_photo_all(2,i_r); ti_sp_d2 = im_all_diff(ti_sp2)
      ti_sp3 = im_photo_all(3,i_r); ti_sp_d3 = im_all_diff(ti_sp3)
      ti_sp4 = im_photo_all(4,i_r); ti_sp_d4 = im_all_diff(ti_sp4)
      ti_sp5 = im_photo_all(5,i_r); ti_sp_d5 = im_all_diff(ti_sp5)
      do concurrent (i_z = 1:n_z)
        rate_prct(i_r,i_z) = freq_prct(i_r,i_z) * den(i_z,ti_sp1)
        ls_p(i_z,ti_sp1) = ls_p(i_z,ti_sp1) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp2) = pr_p(i_z,ti_sp2) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp3) = pr_p(i_z,ti_sp3) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp4) = pr_p(i_z,ti_sp4) + rate_prct(i_r,i_z)
        pr_p(i_z,ti_sp5) = pr_p(i_z,ti_sp5) + rate_prct(i_r,i_z)

        rjac(ti_sp_d1,ti_sp_d1,i_z) = rjac(ti_sp_d1,ti_sp_d1,i_z) &
          - freq_prct(i_r,i_z)
        rjac(ti_sp_d2,ti_sp_d1,i_z) = rjac(ti_sp_d2,ti_sp_d1,i_z) &
          + freq_prct(i_r,i_z)
        rjac(ti_sp_d3,ti_sp_d1,i_z) = rjac(ti_sp_d3,ti_sp_d1,i_z) &
          + freq_prct(i_r,i_z)
        rjac(ti_sp_d4,ti_sp_d1,i_z) = rjac(ti_sp_d4,ti_sp_d1,i_z) &
          + freq_prct(i_r,i_z)
        rjac(ti_sp_d5,ti_sp_d1,i_z) = rjac(ti_sp_d5,ti_sp_d1,i_z) &
          + freq_prct(i_r,i_z)
      end do
    end do

    ! Electron reactions
    !? no Jacobian?
    pr_e = zero; ls_e = zero; rate_erct = zero
    do i_r = 1, n_erct
      ti_sp1 = im_e_all(1,i_r)
      ti_sp2 = im_e_all(2,i_r)
      ti_sp3 = im_e_all(3,i_r)
      ti_sp4 = im_e_all(4,i_r)
      ti_sp5 = im_e_all(5,i_r)
      do concurrent (i_z = iz_bot:n_z)
        rate_erct(i_r,i_z) = freq_erct(i_r,i_z) * den(i_z,ti_sp1)
        ls_e(i_z,ti_sp1) = ls_e(i_z,ti_sp1) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp2) = pr_e(i_z,ti_sp2) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp3) = pr_e(i_z,ti_sp3) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp4) = pr_e(i_z,ti_sp4) + rate_erct(i_r,i_z)
        pr_e(i_z,ti_sp5) = pr_e(i_z,ti_sp5) + rate_erct(i_r,i_z)

!        rjac(ti_sp_d2,ti_sp_d1,i_z) = rjac(ti_sp_d2,ti_sp_d1,i_z) + freq_erct(i_r,i_z)
!        rjac(ti_sp_d3,ti_sp_d1,i_z) = rjac(ti_sp_d3,ti_sp_d1,i_z) + freq_erct(i_r,i_z)
!        rjac(ti_sp_d4,ti_sp_d1,i_z) = rjac(ti_sp_d4,ti_sp_d1,i_z) + freq_erct(i_r,i_z)
!        rjac(ti_sp_d5,ti_sp_d1,i_z) = rjac(ti_sp_d5,ti_sp_d1,i_z) + freq_erct(i_r,i_z)
      end do
    end do

    ! Chemical reactions
    pr_c = zero; ls_c = zero; rate_rct = zero
    do i_r = 1, n_rct
      if (chem_type(i_r) == 1) then    ! unimolecular
        ti_sp1 = ireactant(1,i_r); ti_sp_d1 = im_all_diff(ti_sp1)
        ti_sp3 = ireactant(3,i_r); ti_sp_d3 = im_all_diff(ti_sp3)
        ti_sp4 = ireactant(4,i_r); ti_sp_d4 = im_all_diff(ti_sp4)
        ti_sp5 = ireactant(5,i_r); ti_sp_d5 = im_all_diff(ti_sp5)
        do concurrent (i_z = 1:n_z)
          rate_rct(i_r,i_z) = rk(i_r,i_z) * den(i_z,ti_sp1)
          ls_c(i_z,ti_sp1) = ls_c(i_z,ti_sp1) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp3) = pr_c(i_z,ti_sp3) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp4) = pr_c(i_z,ti_sp4) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp5) = pr_c(i_z,ti_sp5) + rate_rct(i_r,i_z)

          rjac(ti_sp_d1,ti_sp_d1,i_z) = rjac(ti_sp_d1,ti_sp_d1,i_z) &
            - rk(i_r,i_z)
          rjac(ti_sp_d3,ti_sp_d1,i_z) = rjac(ti_sp_d3,ti_sp_d1,i_z) &
            + rk(i_r,i_z)
          rjac(ti_sp_d4,ti_sp_d1,i_z) = rjac(ti_sp_d4,ti_sp_d1,i_z) &
            + rk(i_r,i_z)
          rjac(ti_sp_d5,ti_sp_d1,i_z) = rjac(ti_sp_d5,ti_sp_d1,i_z) &
            + rk(i_r,i_z)
        end do
      else
        ti_sp1 = ireactant(1,i_r); ti_sp_d1 = im_all_diff(ti_sp1)
        ti_sp2 = ireactant(2,i_r); ti_sp_d2 = im_all_diff(ti_sp2)
        ti_sp3 = ireactant(3,i_r); ti_sp_d3 = im_all_diff(ti_sp3)
        ti_sp4 = ireactant(4,i_r); ti_sp_d4 = im_all_diff(ti_sp4)
        ti_sp5 = ireactant(5,i_r); ti_sp_d5 = im_all_diff(ti_sp5)
        do concurrent (i_z = 1:n_z)
          rate_rct(i_r,i_z) = rk(i_r,i_z) * den(i_z,ti_sp1) * den(i_z,ti_sp2)
          ls_c(i_z,ti_sp1) = ls_c(i_z,ti_sp1) + rate_rct(i_r,i_z)
          ls_c(i_z,ti_sp2) = ls_c(i_z,ti_sp2) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp3) = pr_c(i_z,ti_sp3) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp4) = pr_c(i_z,ti_sp4) + rate_rct(i_r,i_z)
          pr_c(i_z,ti_sp5) = pr_c(i_z,ti_sp5) + rate_rct(i_r,i_z)

          t_freq1 = rk(i_r,i_z) * den(i_z,ti_sp2)
          t_freq2 = rk(i_r,i_z) * den(i_z,ti_sp1)
          rjac(ti_sp_d1,ti_sp_d1,i_z) = rjac(ti_sp_d1,ti_sp_d1,i_z) - t_freq1
          rjac(ti_sp_d1,ti_sp_d2,i_z) = rjac(ti_sp_d1,ti_sp_d2,i_z) - t_freq2
          rjac(ti_sp_d2,ti_sp_d1,i_z) = rjac(ti_sp_d2,ti_sp_d1,i_z) - t_freq1
          rjac(ti_sp_d2,ti_sp_d2,i_z) = rjac(ti_sp_d2,ti_sp_d2,i_z) - t_freq2
          rjac(ti_sp_d3,ti_sp_d1,i_z) = rjac(ti_sp_d3,ti_sp_d1,i_z) + t_freq1
          rjac(ti_sp_d3,ti_sp_d2,i_z) = rjac(ti_sp_d3,ti_sp_d2,i_z) + t_freq2
          rjac(ti_sp_d4,ti_sp_d1,i_z) = rjac(ti_sp_d4,ti_sp_d1,i_z) + t_freq1
          rjac(ti_sp_d4,ti_sp_d2,i_z) = rjac(ti_sp_d4,ti_sp_d2,i_z) + t_freq2
          rjac(ti_sp_d5,ti_sp_d1,i_z) = rjac(ti_sp_d5,ti_sp_d1,i_z) + t_freq1
          rjac(ti_sp_d5,ti_sp_d2,i_z) = rjac(ti_sp_d5,ti_sp_d2,i_z) + t_freq2
        end do
      end if
    end do

    ! Condensation
    ! ls_cdn = zero
    ! do i_sp1 = 1, n_diff
      ! i_sp2 = im_diff_all(i_sp1)
      ! do i_z = 1, n_z
        ! t_pvap = vapor(sp_list(i_sp2), Tn(i_z))
        ! t_pprs = den(i_z,i_sp2) * kB * Tn(i_z)
        ! if (t_pprs > t_pvap) then
          ! !? What is 1E12? Both instances of 1E12 below are different quantities
          ! ls_cdn(i_z,i_sp2) = 1.E+12_wp * (t_pprs - t_pvap) ! what is 10^12?
          ! rjac(i_sp1,i_sp1,i_z) = rjac(i_sp1,i_sp1,i_z) - 1.0+12_wp * kB * Tn(i_z)
        ! end if
      ! end do
    ! end do

    pr = pr_p + pr_e + pr_c
    ls = ls_p + ls_e + ls_c

    ! if there are no diffusing species, then just calculate production and
    ! loss rates without updating densities
    if (n_diff == 0) exit

    ! --- Test convergence ----------------------------------------------------
    if ((iter > iter_max) .or. (maxval(abs(F_con)) < tol_den) &
      .or. (maxval(abs(delN)) < tol_bal)) then
      ! write(*,'(" DIFF NEWTON ended after iteration ", I0)') iter-1
      exit
    end if

    ! --- Solve system of continuity equations --------------------------------
    ! solve F(N) = 0 using Newton's method: N' = N - F(N)/(dF/dN)
    ! jac = -dF/dN; delN = N' - N; F_con = F(N)
    ! solve (-dF/dN)(N' - N) = F(N)
    flux = zero
    div_flux = zero
    dNdt = zero
    jac_a = zero
    jac_b = zero
    jac_c = zero
    F_con = zero
    delN = zero

    do i_sp1 = 1, n_diff
      i_sp2 = im_diff_all(i_sp1)

      ! Interior
      do i_z = 2, n_z-1
        flux(i_z,i_sp2) = &
          (alpha(i_z,i_sp1) - beta(i_z,i_sp1)) * den(i_z,i_sp2) &
          - (alpha(i_z,i_sp1) + beta(i_z,i_sp1)) * den(i_z+1,i_sp2)
        div_flux(i_z,i_sp2) = a(i_z,i_sp1) * den(i_z-1,i_sp2) &
          + b(i_z,i_sp1) * den(i_z,i_sp2) + c(i_z,i_sp1) * den(i_z+1,i_sp2)
        dNdt(i_z,i_sp2) = (den(i_z,i_sp2) - den_old(i_z,i_sp2)) / tstep_diff

        F_con(i_sp1,i_z) = pr(i_z,i_sp2) - ls(i_z,i_sp2) &
          - div_flux(i_z,i_sp2) - dNdt(i_z,i_sp2)
        jac_a(i_sp1,i_sp1,i_z) = a(i_z,i_sp1)
        jac_b(i_sp1,i_sp1,i_z) = one / tstep_diff + b(i_z,i_sp1)
        do i_sp3 = 1, n_diff
          jac_b(i_sp1,i_sp3,i_z) = jac_b(i_sp1,i_sp3,i_z) &
            - rjac(i_sp1,i_sp3,i_z)
        end do
        jac_c(i_sp1,i_sp1,i_z) = c(i_z,i_sp1)
      end do

      ! Lower Boundary ?fix with updated equations
      i_z = 1
      if (ibnd(i_sp2,1) == 1) then
      ! diffusive equilibrium
      ! only photochemical production and loss, no flux or div(flux)
        dNdt(i_z,i_sp2) = (den(i_z,i_sp2) - den_old(i_z,i_sp2)) / tstep_diff
        F_con(i_sp1,i_z) = pr(i_z,i_sp2) - ls(i_z,i_sp2) - dNdt(i_z,i_sp2)
        jac_b(i_sp1,i_sp1,i_z) = one / tstep_diff
        do i_sp3 = 1, n_diff
          jac_b(i_sp1,i_sp3,i_z) = jac_b(i_sp1,i_sp3,i_z) &
            - rjac(i_sp1,i_sp3,i_z)
        end do
      else if (ibnd(i_sp2,1) == 2) then
      ! fixed velocity; upwards positive
        flux(i_z,i_sp2) = bval(i_sp2,1) * den(i_z,i_sp2)
        div_flux(i_z,i_sp2) = b(i_z,i_sp1) * den(i_z,i_sp2) &
          + c(i_z,i_sp1) * den(i_z+1,i_sp2)
        dNdt(i_z,i_sp2) = (den(i_z,i_sp2) - den_old(i_z,i_sp2)) / tstep_diff
        F_con(i_sp1,i_z) = pr(i_z,i_sp2) - ls(i_z,i_sp2) &
          - div_flux(i_z,i_sp2) - dNdt(i_z,i_sp2)
!?       jac_b(i_sp1,i_sp1,i_z) = one/tstep_diff + b(i_z,i_sp1)
        jac_b(i_sp1,i_sp1,i_z) = b(i_z,i_sp1)
        do i_sp3 = 1, n_diff
          jac_b(i_sp1,i_sp3,i_z) = jac_b(i_sp1,i_sp3,i_z) &
            - rjac(i_sp1,i_sp3,i_z)
        end do
      else if (ibnd(i_sp2,1) == 3) then
      ! fixed mole fraction
        flux(i_z,i_sp2) = &
          (alpha(i_z,i_sp1) - beta(i_z,i_sp1)) * den(i_z,i_sp2) &
          - (alpha(i_z,i_sp1) + beta(i_z,i_sp1)) * den(i_z+1,i_sp2)
        div_flux(i_z,i_sp2) = b(i_z,i_sp1) * den(i_z,i_sp2) &
          + c(i_z,i_sp1) * den(i_z+1,i_sp2)
        dNdt(i_z,i_sp2) = (den(i_z,i_sp2) - bval(i_sp2,1) * den(i_z,0)) &
          / tstep_diff
        F_con(i_sp1,i_z) = -dNdt(i_z,i_sp2)
        jac_b(i_sp1,i_sp1,i_z) = one / tstep_diff
      else if (ibnd(i_sp2,1) == 4) then
      ! fixed flux; upwards positive
        flux(i_z,i_sp2) = bval(i_sp2,1)
        div_flux(i_z,i_sp2) = flux(i_z,i_sp2) &
          / (rz_mid(i_z-1) / rz(i_z))**two / dr_mid(i_z-1) &
          + b(i_z,i_sp1) * den(i_z,i_sp2) + c(i_z,i_sp1) * den(i_z+1,i_sp2)
        dNdt(i_z,i_sp2) = (den(i_z,i_sp2) - den_old(i_z,i_sp2)) / tstep_diff
        F_con(i_sp1,i_z) = -div_flux(i_z,i_sp2) - dNdt(i_z,i_sp2) !? what happened to pr/ls?
        jac_b(i_sp1,i_sp1,i_z) = b(i_z,i_sp1) !? what happened to 1/t? same as above. also, fixed flux term disappears with d/dN
        do i_sp3 = 1, n_diff
          jac_b(i_sp1,i_sp3,i_z) = jac_b(i_sp1,i_sp3,i_z) &
            - rjac(i_sp1,i_sp3,i_z)
        end do
      end if
      jac_a(i_sp1,i_sp1,i_z) = a(i_z,i_sp1)
      jac_c(i_sp1,i_sp1,i_z) = c(i_z,i_sp1)

      ! Upper Boundary ?fix with updated equations
      i_z = n_z
      dNdt(i_z,i_sp2) = (den(i_z,i_sp2) - den_old(i_z,i_sp2)) / tstep_diff
      if (ibnd(i_sp2,2) == 1) then
      ! diffusive equilibrium
        F_con(i_sp1,i_z) = pr(i_z,i_sp2) - ls(i_z,i_sp2) - dNdt(i_z,i_sp2)
      else if ((ibnd(i_sp2,2) == 2) .or. (ibnd(i_sp2,2) == 3)) then
      ! Jeans or fixed velocity
        flux(i_z,i_sp2) = bval(i_sp2,2) * den(i_z,i_sp2)
        div_flux(i_z,i_sp2) = a(i_z,i_sp1) * den(i_z-1,i_sp2) &
          + b(i_z,i_sp1) * den(i_z,i_sp2)
        F_con(i_sp1,i_z) = pr(i_z,i_sp2) - ls(i_z,i_sp2) &
          - div_flux(i_z,i_sp2) - dNdt(i_z,i_sp2)
      else if (ibnd(i_sp2,2) == 4) then
      ! fixed flux; upwards positive
        flux(i_z,i_sp2) = bval(i_sp2,2)
        div_flux(i_z,i_sp2) = flux(i_z,i_sp2) &
          / (rz_mid(i_z-1) / rz(i_z))**two / dr_mid(i_z-1) &
          + a(i_z,i_sp1) * den(i_z-1,i_sp2) + b(i_z,i_sp1) * den(i_z,i_sp2)
        F_con(i_sp1,i_z) = -div_flux(i_z,i_sp2) - dNdt(i_z,i_sp2) !? what happened to pr/ls?
      end if
      jac_a(i_sp1,i_sp1,i_z) = a(i_z,i_sp1)
      jac_b(i_sp1,i_sp1,i_z) = b(i_z,i_sp1) + one / tstep_diff
      do i_sp3 = 1, n_diff
        jac_b(i_sp1,i_sp3,i_z) = jac_b(i_sp1,i_sp3,i_z) - rjac(i_sp1,i_sp3,i_z)
      end do
      jac_c(i_sp1,i_sp1,i_z) = c(i_z,i_sp1)
    end do

    call solve_tridiag(n_diff, n_z, jac_a, jac_b, jac_c, F_con, delN)

    ! --- Update densities ----------------------------------------------------
    do concurrent (i_sp1 = 1:n_diff, i_z = 1:n_z)
      ! limit maximum change to ensure smooth convergence
      if (delN(i_sp1,i_z) < 0) then
        delN(i_sp1,i_z) = &
          max(delN(i_sp1,i_z), -damping * den(i_z,im_diff_all(i_sp1)))
      else
        delN(i_sp1,i_z) = &
          min(delN(i_sp1,i_z), damping * den(i_z,im_diff_all(i_sp1)))
      end if
      den(i_z,im_diff_all(i_sp1)) = &
        max(den(i_z,im_diff_all(i_sp1)) + delN(i_sp1,i_z), eps)
    end do

    ! update total density
    do concurrent (i_z = 1:n_z)
      den(i_z,0) = sum(den(i_z,1:n_sp-1))
    end do

  end do NEWTON

  deallocate(rjac, jac_a, jac_b, jac_c, F_con, delN)

end subroutine solve_diff


subroutine solve_tridiag (n_diff, n_z, a, b, c, d, delN)
! This subroutine implements the tridiagonal matrix algorithm for solving
! [ b(1) c(1)  0    0   ...   0  ]  [ x(1) ]   [ d(1) ]
! [ a(2) b(2) c(2)  0   ...   0  ]  [ x(2) ]   [ d(2) ]
! [  0   a(3) b(3) c(3) ...   0  ]  [ x(3) ] = [ d(3) ]
! [           ...                ]  [ ...  ]   [ ...  ]
! [  0    0    0    0   ... b(n) ]  [ x(n) ]   [ d(n) ]

! a, b, and c have dimensions (n_diff, n_diff, n_z), while d and delN have
! dimensions (n_diff, n_z)

! The algorithm solves for x using Gaussian elimination. This is achieved by
! calculating in the forward direction
! c'(1) = c(1)/b(1)
! c'(i) = c(i)/(b(i) - a(i)*c'(i-1))
! d'(1) = d(1)/b(1)
! d'(i) = (d(i)-a(i)*d'(i-1))/(b(i) - a(i)*c'(i-1))

! The solution is then calculated backwards via
! x(n) = d'(n)
! x(i) = d'(i) - c'(i)*x(i+1)

! delN holds d' until it is converted to x at the last step
  use types, only: wp => dp
  use constants
  use lapack95, only: getrf, getrs

  implicit none

  integer, intent(in) :: n_diff, n_z
  real(wp), intent(in), dimension(:,:,:) :: a, b, c
  real(wp), intent(in), dimension(:,:) :: d
  real(wp), intent(out), dimension(:,:) :: delN

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  ! pre-factor matrix for converting c'->c and d'->d
  real(wp), dimension(n_diff,n_diff) :: prefac
  ! row pivot information from getrf
  integer, dimension(n_diff) :: ipiv
  ! error code for getrf and getrs
  integer :: info
  ! gam(n) is c'(n-1)
  real(wp), dimension(n_diff,n_diff,n_z) :: gam
  ! real(wp), dimension(n_diff,n_z) :: dp, d_sav, delp ! optimization vars

  ! loop variables
  integer :: i_z, i_sp1, i_sp2
  ! temporary / dummy variables
  real(wp), dimension(n_diff) :: sol

  !----------------------------------------------------------------------------
  !  Start Tridiagonal Matrix Algorithm
  !----------------------------------------------------------------------------

  prefac = zero
  gam = zero
  delN = zero

  ! --- 1. Calculate d'(1), i.e., solve b(1)*delN(1)=d(1) ---------------------

  i_z = 1
  prefac = b(:,:,i_z)

  call getrf (prefac, ipiv, info)
  ! The getrf(A, ipiv, info) routine computes the LU factorization of matrix A
  ! into P*L*U, where P is a permutation matrix, L is lower triangular with
  ! unit diagonal elements, and and U is upper triangular. A will be!
  ! overwritten by the factors L and U from the factorization of A = P*L*U; the
  ! unit diagonal elements of L are not stored. ipiv encodes the row swap
  ! operations. info is the error code.
  if (info /= 0) then
    write(*,'("solve_tridiag: Error in getrf at i_z = ", I4, "! Code: ", &
      I3)') i_z, info
    return
  end if
  sol = d(:,i_z)
  call getrs (prefac, ipiv, sol, 'N', info)
  ! The getrs(a, ipiv, b, trans, info) routine, with the 'N' setting, solves
  ! for X in the system of linear equations A*X = B. The array a contains LU
  ! factorization of A returned by getrf, and the ipiv array encodes the row
  ! swap operations returned by getrf. b is overwritten by the solution matrix
  ! X. info is the error code.
  if (info /= 0) then
    write(*,'("solve_tridiag: Error in getrs at i_z = ", I4, "! Code: ", &
      I3)') i_z, info
    return
  end if
  delN(:,i_z) = sol

  ! --- 2. Solve for other altitudes ------------------------------------------
  do i_z = 2, n_z

    ! 2a. Calculate c'

    do i_sp1 = 1, n_diff
      gam(i_sp1,i_sp1,i_z) = c(i_sp1,i_sp1,i_z-1)
      sol = gam(:,i_sp1,i_z) !? should this be (i_sp1,:,i_z) instead?
      call getrs (prefac, ipiv, sol, 'N', info)
      if (info /= 0) then
        write(*,'("solve_tridiag: Error in getrs at i_z = ", I4, "! Code: ", &
          I3)') i_z, info
        return
      end if
      gam(:,i_sp1,i_z) = sol !? should this be (i_sp1,:,i_z) instead?
    end do

    ! 2b. Update prefac for i_z>1

    do concurrent (i_sp1 = 1:n_diff, i_sp2 = 1:n_diff)
      prefac(i_sp1,i_sp2) = b(i_sp1,i_sp2,i_z) &
        - a(i_sp1,i_sp1,i_z) * gam(i_sp1,i_sp2,i_z) ! note that a is diagonal
    end do

    ! 2c. Calculate d'

    do concurrent (i_sp1 = 1:n_diff)
      delN(i_sp1,i_z) = d(i_sp1,i_z) - a(i_sp1,i_sp1,i_z) * delN(i_sp1,i_z-1)
    end do

    call getrf (prefac, ipiv, info)
    if (info /= 0) then
      write(*,'("solve_tridiag: Error in getrf at i_z = ", I4, "! Code: ", &
        I3)') i_z, info
      return
    end if
    sol = delN(:,i_z)
    call getrs (prefac, ipiv, sol, 'N', info)
    if (info /= 0) then
      write(*,'("solve_tridiag: Error in getrs at i_z = ", I4, "! Code: ", &
        I3)') i_z, info
      return
    end if
    delN(:,i_z) = sol

  end do

  ! --- 3. Convert d' to x ----------------------------------------------------

  do i_z = n_z-1, 1, -1
    do i_sp1 = 1, n_diff
      delN(i_sp1,i_z) = delN(i_sp1,i_z) &
        - dot_product(gam(i_sp1,:,i_z+1), delN(:,i_z+1))
    end do
  end do

  !----------------------------------------------------------------------------
  !  Improve Solution, Create error vector (currently disabled)
  !----------------------------------------------------------------------------

  ! if (iopt >= 1) then

  ! d_sav = d
  ! nl = 1
  ! do i = 1, n_diff
    ! sm = zero
    ! do i_sp1 = 1, n_diff
      ! sm = sm + b(i,j,nl)*delN(j,nl)+c(i,j,nl)*delN(j,nl+1)
    ! end do
    ! dp(i,nl)=sm
  ! end do

  ! do nl = 2, n_z-1
  ! do i = 1, n_diff
    ! sm = zero
    ! do j = 1, n_diff
      ! sm = sm + a(i,j,nl)*delN(j,nl-1)+b(i,j,nl)*delN(j,nl)+c(i,j,nl)*delN(j,nl+1)
    ! end do
    ! dp(i,nl)=sm
  ! end do
  ! end do

  ! nl = n_z
  ! do i = 1, n_diff
    ! sm = zero
    ! do j = 1, n_diff
      ! sm = sm + a(i,j,nl)*delN(j,nl-1)+b(i,j,nl)*delN(j,nl)
    ! end do
    ! dp(i,nl)=sm
  ! end do

  ! dp = dp - d_sav

  ! .. solve A*delN=dp

  ! prefac = zero
  ! gam = zero
  ! delp = zero

  ! nl = 1
  ! prefac(1:n_diff,1:n_diff) = b(1:n_diff,1:n_diff,nl)
  ! sol(1:n_diff) = dp(1:n_diff,nl)

  ! call GETRF ( prefac(1:n_diff,1:n_diff), ipiv(1:n_diff))
  ! call GETRS ( prefac(1:n_diff,1:n_diff), ipiv(1:n_diff), sol(1:n_diff) )
  ! delp(1:n_diff,nl) = sol(1:n_diff)

  !
  !  ***** frontwards
  !

  ! do nl = 2, n_z

    ! .. solve dum*gam=c

    ! do j = 1, n_diff
      ! gam(j,j,nl) = c(j,j,nl-1)
    ! end do

    ! do j = 1, n_diff
      ! sol(1:n_diff) = gam(1:n_diff,j,nl)
      ! call GETRS ( prefac(1:n_diff,1:n_diff), ipiv(1:n_diff), sol(1:n_diff) )
      ! gam(1:n_diff,j,nl) = sol(1:n_diff)
    ! end do

    ! .. dum = b(nl)-a(nl)*gam(nl)

    ! do j = 1, n_diff
      ! do i = 1, n_diff
        ! prefac(i,j) = b(i,j,nl) - a(i,i,nl) * gam(i,j,nl)
      ! end do
    ! end do

    ! .. solve dum*delN(nl)=d(nl)-a(nl)*delN(nl-1)

    ! do i = 1, n_diff
      ! delp(i,nl) = dp(i,nl) - a(i,i,nl) * delp(i,nl-1)
    ! end do

    ! call GETRF (prefac(1:n_diff,1:n_diff), ipiv(1:n_diff) )
    ! sol(1:n_diff) = delp(1:n_diff,nl)
    ! call GETRS (prefac(1:n_diff,1:n_diff), ipiv(1:n_diff), sol(1:n_diff) )
    ! delp(1:n_diff,nl) = sol(1:n_diff)

  ! end do

  !
  !  ***** backwards
  !

  ! do nl = n_z-1, 1, -1

    ! do i = 1, n_diff
      ! sm = zero
      ! do m = 1, n_diff
        ! sm = sm + gam(i,m,nl+1) * delp(m,nl+1)
      ! end do
      ! delp(i,nl) = delp(i,nl) - sm
    ! end do

  ! end do

  ! delN = delN - delp

  ! end if


end subroutine solve_tridiag


end module solve
