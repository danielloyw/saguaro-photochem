module read_atmos_mod

contains

subroutine read_atmos
! This subroutine reads in the model atmosphere from atm1D.in, and calculates 
! other atmospheric state quantities. 

  !----------------------------------------------------------------------------
  !  Modules
  !----------------------------------------------------------------------------
  
  use types, only: wp => dp
  use constants
  use global_variables
  use utils, only : find_name, locate

  !----------------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------------
  
  implicit none
  ! file unit numbers for atm1D.in
  integer :: fid
  ! minimum density
  real(wp), parameter :: den_min = 1.0e-50_wp
  
  ! species names in atm1D.in
  character(len=12), allocatable, dimension(:) :: sp_atm
  ! index mapping from sp_atm to sp_list
  integer, allocatable, dimension(:) :: im_atm_list
  ! does species from settings have density defined in atmosphere?
  logical, allocatable, dimension(:) :: has_den
  
  ! loop variables
  integer :: i_sp, i_z
  ! temporary / dummy variables
  character(len=128) :: ts_header
  character(len=12) :: ts_name
  integer :: tn_sp1, tn_sp2
  
  !----------------------------------------------------------------------------
  !  Read in model atmosphere from atm1D
  !----------------------------------------------------------------------------
  
  ! start with all species undefined, then flip when density found
  allocate(has_den(n_sp))
  has_den = .false.
  
  ! read in model atmosphere file atm1D.in
  open(newunit=fid, file='atm1D.in', status='old', action='read')
    
    read(fid,*) n_z, tn_sp1
    ! allocate variables
    allocate(sp_atm(tn_sp1), im_atm_list(tn_sp1))
    allocate(z(n_z), rz(0:n_z+1), grv(n_z), eK(n_z))
    allocate(Tn(n_z), Ti(n_z), Te(n_z), prs(n_z), mass(n_z), rho(n_z))
    allocate(den(n_z,0:n_sp), den_old(n_z,0:n_sp), vmr(n_z,n_sp))

    ! read species and associate them with sp_list from species settings files
    ! i.e., nmolecules.settings and imolecules.settings
    read(fid,'(A)') ts_header
    read(fid,"(10(2X,A12,1X))") (sp_atm(i_sp),i_sp=1,tn_sp1)
    do concurrent (i_sp = 1:tn_sp1)
      tn_sp2 = find_name(sp_atm(i_sp),sp_list)
      im_atm_list(i_sp) = tn_sp2
      if(tn_sp2 /= 0) has_den(tn_sp2) = .true.
    end do

    ! read in model atmosphere and populate relevant variables
    read(fid,'(A)') ts_header
    read(fid,*) (z(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (rz(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (grv(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (Tn(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (Te(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (prs(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (rho(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (mass(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (eK(i_z), i_z=1, n_z)
    read(fid,'(A)') ts_header
    read(fid,*) (den(i_z,0), i_z=1, n_z) ! den(*,0) is for overall density
    do i_sp = 1, tn_sp1
      tn_sp2 = im_atm_list(i_sp)
      read(fid,'(A)') ts_name
      read(fid,*) (den(i_z,tn_sp2), i_z=1, n_z)
    end do
    
  close(unit=fid)

  Ti = 230._wp !? need to confirm
  
  ! altitude to opaque atmosphere layer (e.g., aerosols)
  z_bot = 0.E5_wp
  iz_bot = locate(z_bot, z)
  
  !----------------------------------------------------------------------------
  !  Recalculate quanitites to ensure consistency
  !----------------------------------------------------------------------------

  ! apply bottom boundary condition of stipulated mole ratio
  do concurrent (i_sp = 1:n_diff)
    tn_sp2 = im_diff_all(i_sp)
    if(ibnd(tn_sp2,1) == 3) then
      den(1,tn_sp2) = bval(tn_sp2,1)*den(1,0)
    end if
  end do

  ! set density to min value if species does not have density defined
  do concurrent (i_sp = 1:n_sp-1)
    if(.not. has_den(i_sp)) then
      den(:,i_sp) = den_min
      write(*,"('Density not defined for ', A, &
        '... Assigning minimum density of ', ES8.1)") sp_list(i_sp), den_min
    end if
  end do

  ! set density to min value if below minimum
  where (den < den_min) den = den_min
  
  ! reset electron density to sum of ion densities
  do concurrent (i_z = 1:n_z)
    den(i_z,n_sp) = sum(den(i_z, n_neu+1:n_sp-1))
  end do
  
  ! reset total density
  do concurrent (i_z = 1:n_z)
    den(i_z,0) = sum(den(i_z,1:n_sp))
  end do
  
  ! set density to fixed mole fraction according to settings
  do concurrent (i_sp = 1:n_sp-1)
    if(istat(i_sp) == 3) then
      den(:,i_sp) = bval(i_sp,1)*den(:,0)
    end if
  end do

  ! calculate mole fraction
  do concurrent (i_z = 1:n_z)
    do concurrent (i_sp = 1:n_sp)
      vmr(i_z,i_sp) = den(i_z,i_sp)/den(i_z,0)
    end do
  end do

  ! reset gravity
  do concurrent (i_z = 1:n_z)
    grv(i_z) = GM/rz(i_z)**2
  end do

  ! reset pressure to be consistent with N & T
  do concurrent (i_z = 1:n_z)
    prs(i_z) = kB*Tn(i_z)*den(i_z,0)
  end do

  ! calculate mass density, mean molecular mass
  do concurrent (i_z = 1:n_z)
    if (n_ion > 0) then
      rho(i_z) = dot_product(mmw,den(i_z,1:n_sp)) * amu
      mass(i_z) = rho(i_z) / (den(i_z,0) + two*den(i_z,n_sp)) / amu !? why
    else
      rho(i_z) = dot_product(mmw,den(i_z,1:n_sp)) * amu
      mass(i_z) = rho(i_z) / den(i_z,0) / amu
    end if
  end do

  !----------------------------------------------------------------------------
  !  Auxillary altitude arrays
  !----------------------------------------------------------------------------

  allocate(rz_mid(0:n_z), dr_mid(0:n_z), dr(n_z))
  allocate(Ht(n_z,0:n_sp), col(n_z))

  rz(0) = rz(1) - (rz(2) - rz(1))
  rz(n_z+1) = rz(n_z) + (rz(n_z) - rz(n_z-1))
  do concurrent (i_z = 0:n_z)
    rz_mid(i_z) = half * (rz(i_z+1) + rz(i_z))
    dr_mid(i_z) = rz(i_z+1) - rz(i_z)
  end do

  do concurrent (i_z = 1:n_z)
    dr(i_z) = rz_mid(i_z) - rz_mid(i_z-1)
  end do

  !----------------------------------------------------------------------------
  !  Scale heights
  !----------------------------------------------------------------------------
  
  ! for each species
  do concurrent (i_z = 1:n_z, i_sp = 1:n_sp-1)
    Ht(i_z,i_sp) = kB * tn(i_z) / grv(i_z) / mmw(i_sp) / amu
  end do
  
  ! mean scale height
  do concurrent (i_z = 1:n_z)
    Ht(i_z,0) = kB * tn(i_z) / grv(i_z) / mass(i_z) / amu
  end do

  !----------------------------------------------------------------------------
  !  Column density above
  !----------------------------------------------------------------------------

  col(n_z) = den(n_z,0) * Ht(n_z,0)
  do i_z = n_z-1, 1, -1
    col(i_z) = col(i_z+1) + half * (den(i_z,0) + den(i_z+1,0)) * dr_mid(i_z)
  end do

  deallocate(sp_atm, im_atm_list, has_den)
end subroutine read_atmos

end module
