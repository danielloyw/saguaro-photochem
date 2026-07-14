program saguaro

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------

  use types, only: wp => dp
  use constants
  use global_variables
  use read_inputs
  use diffusion, only: diffusion_setup
  use photo, only: photo_setup, cal_photo_rates
  use electrons, only: electrons_setup, cal_electrons
  use solve, only: solve_chem, solve_diff
  use output, only: write_output

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------

  implicit none
  ! file unit numbers
  integer :: fid_runid, fid_tctl, fid_solar, fid_run

  ! max diffusion loop iterations
  integer :: max_iter_diff
  ! max chemistry loop iterations per diffusion loop iteration
  integer :: max_iter_chem
  ! iteration count for diffusion loop
  integer :: iter_diff
  ! iteration count for chemistry loop
  integer :: iter_chem
  ! number of iterations before recalculating photolysis and photoelectrons
  integer :: n_iter_sol_recal
  ! steady-state test variables for diffusion loop
  real(wp) :: tol_diff, test_diff
  ! steady-state test variables for chemistry loop
  real(wp) :: tol_chem, test_chem
  ! altitude bin and index of species where maximum error occurs
  integer :: test_z, test_sp
  
  ! loop variables
  integer :: i_sp1, i_sp2, i_z
  ! temporary / dummy variables
  character(len=128) :: ts_header, ts_setting
  real(wp) :: t, t_prls

  !----------------------------------------------------------------------------
  !  Run Settings
  !----------------------------------------------------------------------------

  ! get run ID from file created from bash script
  open(newunit=fid_runid, file='runid.out', status='old', action='read')
     read(fid_runid,"(A)") runid
  close(unit=fid_runid)

  ! get time control settings from file created from bash script
  open(newunit=fid_tctl, file='tctl.out', status='old', action='read')
     read(fid_tctl,*) tstep_diff, max_iter_diff, tstep_chem, max_iter_chem
  close(unit=fid_tctl)

  if (tstep_diff <= tstep_chem * max_iter_chem * ten) then
    write(*,'("Warning: Large time step for chemistry loop! t_diff = ", &
      ES8.1, ", t_chem = ", ES8.1, ", max. chemistry iterations = ", I0)') &
      tstep_diff, tstep_chem, max_iter_chem
  end if
  
  ! settings for photolysis
  open(newunit=fid_solar, file='solar.settings', status='old', action='read')
     read(fid_solar,*) ts_header
     read(fid_solar,*) cos_sza
     sin_sza = sqrt(one - cos_sza * cos_sza)
     read(fid_solar,*) ts_header
     read(fid_solar,*) diurnal_average
  close(unit=fid_solar)

  ! other run settings
  open(newunit=fid_run, file='run.settings', status='old', action='read')
    read(fid_run, '(A)') ts_header ! recal_photo_diff
    read(fid_run, '(A)') ts_setting
    if (trim(adjustl(ts_setting)) == 'on') then
      recal_photo_diff = .true.
    else
      recal_photo_diff = .false.
    end if
    
    read(fid_run, '(A)') ts_header ! tol_diff
    read(fid_run, *) tol_diff
    
    read(fid_run, '(A)') ts_header ! recal_photo_chem
    read(fid_run, '(A)') ts_setting
    if (trim(adjustl(ts_setting)) == 'on') then
      recal_photo_chem = .true.
    else
      recal_photo_chem = .false.
    end if
    
    read(fid_run, '(A)') ts_header ! tol_chem
    read(fid_run, *) tol_chem
    
    read(fid_run, '(A)') ts_header ! print_chem
    read(fid_run, '(A)') ts_setting
    if (trim(adjustl(ts_setting)) == 'on') then
      print_chem = .true.
    else
      print_chem = .false.
    end if
    
    read(fid_run, '(A)') ts_header ! n_iter_sol_recal
    read(fid_run, *) n_iter_sol_recal

  close(unit=fid_run)

  !----------------------------------------------------------------------------
  !  Prepare Calculations
  !----------------------------------------------------------------------------
  
  write(*,'("Preparing calculations...")')
  
  ! read in settings for species
  call read_species

  ! read model atmosphere
  call read_atmos

  ! read in reaction rate data
  call read_reactions

  ! calculate diffusion coefficients
  call diffusion_setup

  ! set up photolysis calculation
  call photo_setup

  ! set up suprathermal electron calculation
  call electrons_setup

  ! production and loss from photo processes
  allocate(pr_p(n_z,0:n_sp), ls_p(n_z,0:n_sp))
  ! production and loss from electron processes
  allocate(pr_e(n_z,0:n_sp), ls_e(n_z,0:n_sp))
  ! production and loss from chemical processes
  allocate(pr_c(n_z,0:n_sp), ls_c(n_z,0:n_sp))
  ! net production and loss
  allocate(pr(n_z,0:n_sp), ls(n_z,0:n_sp))
  ! condensation loss
  ! allocate(ls_cdn(n_z,0:n_sp))
  ! flux
  allocate(flux(n_z,0:n_sp))
  ! flux divergence
  allocate(div_flux(n_z,0:n_sp))
  ! net balance
  allocate(dNdt(n_z,0:n_sp))

  !----------------------------------------------------------------------------
  !  Time Loops
  !----------------------------------------------------------------------------
  ! Diffusion is much slower than chemistry. Progress of time in the model is
  ! marked by the iterations of the TIME loop, which calculates the diffusion
  ! that occurs within the time interval spanned by the iteration. The progress
  ! of time is marked by moving den to den_old. After diffusion is calculated
  ! within each iteration of the TIME loop, chemistry is calculated, and is
  ! looped to chemical equilibrium. This second CHEM loop does not represent
  ! progress in the global model time. This approach reduces the numerical
  ! errors arising from sharp changes in composition as a result of the
  ! chemistry. 

  write(*,'("Starting time loop...")')
  
  iter_diff = 0
  TIME: do
    iter_diff = iter_diff + 1
    if (iter_diff > max_iter_diff) exit
    ! currently no equilibrium condition for ending TIME loop

    den_old = den

    ! occasionally recalculate photolysis rate of region B
    if (mod(iter_diff-1, n_iter_sol_recal) == 0) then
      cal_photoA = .true.
      cal_photoB = .true.
      cal_photoC = .true.
      cal_photoJ = .true.

      call cal_photo_rates
      call cal_electrons
    end if

    ! --- Calculate Chemistry -------------------------------------------------
    cal_photoA = .true.
    cal_photoB = .false.
    cal_photoC = .true.
    cal_photoJ = .true.

    iter_chem = 0
    test_chem = ten*tol_chem
    CHEM: do
      iter_chem = iter_chem + 1
      if ((iter_chem > max_iter_chem) .or. (test_chem < tol_chem)) exit

      den_old = den
      call solve_chem

      ! check convergence
      test_chem = zero
      test_sp = 1
      test_z = 1
      do i_sp1 = 1, n_chem
        i_sp2 = im_chem_all(i_sp1)
        do i_z = 1, n_z
          if (den(i_z,i_sp2) > 1._wp) then
            t_prls = max(pr(i_z,i_sp2), ls(i_z,i_sp2))
            if (t_prls > 1.E-6_wp) then
              t = abs(pr(i_z,i_sp2) - ls(i_z,i_sp2) - dNdt(i_z,i_sp2)) / t_prls
            else
              t = zero
            end if
            if (t > test_chem) then
              test_chem = t
              test_z = i_z
              test_sp = i_sp2
            end if
          end if
        end do
      end do

      if (print_chem .and. (n_chem > 0)) then
        write(*,921) iter_chem, sp_list(test_sp), &
          z(test_z) * cm_to_km, den(test_z,test_sp), &
          pr(test_z,test_sp), ls(test_z,test_sp), dNdt(test_z,test_sp), &
          pr(test_z,test_sp) - ls(test_z,test_sp) - dNdt(test_z,test_sp), &
          test_chem
      end if
      
921 format(2X, "CHEM: ITER = ", I3, "   Species = ", A12, &
          "   Z = ", F4.0, "   N =", ES11.3, &
          "   P =", ES11.3, "   L =", ES11.3, "   dN/dt =", ES11.3, &
          "   Err =", ES11.3, &
          "   Relerr =", ES11.3)

    end do CHEM

    ! --- Calculate Diffusion -------------------------------------------------
    cal_photoA = .true.
    cal_photoB = .false.
    cal_photoC = .true.
    cal_photoJ = .true.
    
    call solve_diff
!    call HYDROST

    do concurrent (i_z = 1:n_z)
      den(i_z,0) = sum(den(i_z,1:n_sp-1))
      do concurrent (i_sp1 = 1:n_sp)
        vmr(i_z,i_sp1) = den(i_z,i_sp1) / den(i_z,0)
      end do
      rho(i_z) = dot_product(mmw, den(i_z,1:n_sp)) * amu
      mass(i_z) = rho(i_z) / den(i_z,0) / amu
    end do

    ! check equilibrium
    test_diff = zero
    test_z = 1
    test_sp = 1

    do i_sp1 = 1, n_diff
      i_sp2 = im_diff_all(i_sp1)
      do i_z = 1, n_z
        if (den(i_z,i_sp2) > 1._wp) then
          t = abs(den(i_z,i_sp2) - den_old(i_z,i_sp2)) / den(i_z,i_sp2)
        else
          t = zero
        end if
        if (t > test_diff) then
          test_diff = t
          test_z = i_z
          test_sp = i_sp2
        end if
      end do
    end do

    do i_sp1 = 1, n_chem
      i_sp2 = im_chem_all(i_sp1)
      do i_z = 1, n_z
        if (den(i_z,i_sp2) > 1._wp) then
          t = abs(den(i_z,i_sp2) - den_old(i_z,i_sp2)) / den(i_z,i_sp2)
        else
          t = zero
        end if
        if (t > test_diff) then
          test_diff = t
          test_z = i_z
          test_sp = i_sp2
        end if
      end do
    end do

    write(*,922) iter_diff, sp_list(test_sp), z(test_z) * cm_to_km, &
      den(test_z,test_sp), den(test_z,test_sp) - den_old(test_z,test_sp), &
      test_diff

922 format("TIME Step =", I3, "    Species = ", A12, "Z = ", F4.0, &
          "   N =", ES11.3, "   DelN =", ES11.3, &
          "   Rerr =", ES11.3)

  end do TIME

  call write_output

end program
