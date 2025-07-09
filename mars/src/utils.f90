module utils
  
interface
  subroutine gauleg(x1, x2, x, w)
    use types, only: wp => dp
    real(wp), intent(in) :: x1, x2
    real(wp), dimension(:), intent(out) :: x, w
  end subroutine gauleg

end interface
  
  
contains
  
integer pure function find_name(sp_name, sp_list)
! This function finds the index of a species with name "sp_name" in the list
! of all model species "sp_list". 
! Note: sp_list starts at index 0.
  implicit none
  character(len=12), intent(in) :: sp_name
  character(len=12), intent(in), dimension(0:) :: sp_list
  integer :: n_sp ! size of sp_list
  integer :: i_sp ! loop variable

  n_sp = size(sp_list,1)-1 ! -1 because index starts at 0
  find_name = -1
  do i_sp = 0, n_sp
    if(trim(adjustl(sp_name)) == trim(adjustl(sp_list(i_sp)))) then
      find_name = i_sp
      exit
    end if
  end do
end function find_name

integer pure function locate(x, x_list)
! This function finds the index of the element in x_list that is nearest in
! value to x. x_list is monotonically increasing or decreasing.
! Uses a binary search algorithm. x is tested against the mid-point of the
! search window, and then the upper or lower window is picked. 
! Returns -1 if x_list is empty. Returns extreme values if x is out of range.
  use types, only: wp => dp
  implicit none
  real(wp), intent(in) :: x
  real(wp), intent(in), dimension(:) :: x_list
  logical :: is_increasing ! is list ascending?
  integer :: low, high, mid ! index bounds of search window

  ! initialize result to 'not found' and handle empty list
  locate = -1
  if (size(x_list) == 0) return

  ! get actual bounds of the array
  low = lbound(x_list, 1)
  high = ubound(x_list, 1)

  ! determine if increasing
  is_increasing = (x_list(low) <= x_list(high))

  ! handle x when out of x_list range (also handles single-element x_list)
  if (is_increasing) then
    if (x <= x_list(low)) then
      locate = low
      return
    end if
    if (x >= x_list(high)) then
      locate = high
      return
    end if
  else
    if (x >= x_list(low)) then
      locate = low
      return
    end if
    if (x <= x_list(high)) then
      locate = high
      return
    end if
  end if

  ! binary search to find indices bracketing x
  do while (high - low > 1) ! low and high not adjacent indices
    mid = (high + low) / 2

    ! exact match found
    if (x_list(mid) == x) then
      locate = mid 
      return
    end if

    if (is_increasing) then
      if (x_list(mid) < x) then
        low = mid
      else
        high = mid
      end if
    else
      if (x_list(mid) > x) then
        low = mid
      else
        high = mid
      end if
    end if
  end do

  ! now x is between x_list(low) and x_list(high)
  ! find whether x_list(low) and x_list(high) is closer to x
  if (abs(x - x_list(low)) <= abs(x - x_list(high))) then
    locate = low
  else
    locate = high
  end if

end function locate
  

subroutine intrp(xp, yp, x, y)
! This subroutine calculate the values y at sample points x using a linear
! interpolation of points (xp, yp). xp needs to be monotonically increasing. 
! 0 is returned when extrapolated beyond the range spanned by xp.
  use types, only: wp => dp
  implicit none
  real(wp), intent(in), dimension(:) :: xp, yp
  real(wp), intent(in), dimension(:) :: x
  real(wp), intent(out), dimension(:) :: y
  integer :: n_xp, n_x
  ! indices of points on the two sides of x, xp(n1)<x<xp(n2)
  integer :: n1, n2 
  ! loop variables
  integer :: i
  
  n_xp = size(xp)
  n_x = size(x)
  
  do i = 1, n_x
    ! outside low range
    if (x(i) < xp(1)) then
      y(i) = 0
    ! outside high range
    else if (x(i) > xp(n_xp)) then
      y(i) = 0
    else
      ! determine adjacent xp values
      n1 = locate(x(i), xp)
      if (x(i) - xp(n1) >= 0) then
        n2 = n1 + 1
      else
        n2 = n1
        n1 = n2 - 1
      end if
      y(i) = yp(n1) + (yp(n2)-yp(n1)) / (xp(n2)-xp(n1)) * (x(i)-xp(n1))
    end if
  end do

end subroutine intrp



  subroutine inDEXX(arrin, indx)
    use types, only: wp => dp
    implicit none
    real(wp), intent(in), dimension(:) :: arrin
    integer, intent(out), dimension(:) :: indx
    integer :: i, j, l, n, ir, indxt
    real(wp) :: q
    n = SIZE(arrin)
    do 11 j = 1, n
    indx(j) = j
 11 continue
    l = (n / 2) + 1
    ir = n
 10 continue
    if (l .gt. 1) then
      l = l - 1
      indxt = indx(l)
      q = arrin(indxt)
    else
      indxt = indx(ir)
      q = arrin(indxt)
      indx(ir) = indx(1)
      ir = ir - 1
    if (ir .eq. 1) then
      indx(1) = indxt
      return 
    end if
  end if
      i = l
      j = l + l
   20 if (j .le. ir) then
      if (j .lt. ir) then
    if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
    end if
    if (q .lt. arrin(indx(j))) then
      indx(i) = indx(j)
      i = j
      j = j + j
    else
      j = ir + 1
    end if
    goto 20
    end if
    indx(i) = indxt
    goto 10
  end subroutine inDEXX




! Calculates the vapor pressure of sp_list at temperature T.

  function VAPOR( sp_list, T)
    use types, only: wp => dp
    implicit none
    real(wp), intent(in) :: T
    CHARACTER(LEN=*), intent(in) :: sp_list
    real(wp) :: VAPOR

    if(TRIM(sp_list) == 'C2H2') then
!  based on Tickner & Losing 1951 data
       VAPOR = 1.333_wp*10._wp**(9.25_wp-1201.75_wp/T)
       return
    else if(TRIM(sp_list) == 'C2H4') then
       if (T < 104._wp) then
       VAPOR = 1333._wp*10._wp**(8.724_wp - 901.6_wp/(T-2.555_wp))
       else if ((T >= 104._wp) .AND. (T < 120._wp)) then
       VAPOR = 1333._wp*10._wp**(50.79_wp - 1703._wp/T - 17.141_wp*LOG10(T))
       else if (T >= 120._wp) then
       VAPOR = 1333._wp*10._wp**(6.74756_wp-585._wp/(T-18.16_wp))
       end if
       return
    else if(TRIM(sp_list) == 'C2H6') then
       if (T < 90._wp) then
       VAPOR = 1333._wp*10._wp**(10.01_wp - 1085._wp/(T-0.561_wp))
       else if (T > 90.) then
       VAPOR = 1333._wp*10._wp**(5.9366_wp-1086.17_wp/T+3.83464_wp*LOG10(1000._wp/T))
       end if
       return
    else if(TRIM(sp_list) == 'C3H8') then
!  based on Tickner & Losing 1951 data
       VAPOR = 1333._wp*10._wp**(8.16173_wp-1176._wp/T)
       return
    else if(TRIM(sp_list) == 'C4H2') then
!  from Moses 1992
       VAPOR = 1333._wp*10._wp**(5.3817_wp-3300.5_wp/T + 16.63415_wp*LOG10(1000._wp/T))
       return
    else if(TRIM(sp_list) == 'C4H6') then
!  from Moses 1992
!       WRITE(*,*) ' CALCULATE C4H10 VP, T = ',T
!       WRITE(*,*) ' CALCULATE C4H10 VP, ARG = ',8.446_wp-1461.2_wp/T
       VAPOR = 1333._wp*10._wp**(8.032581_wp-1441.42_wp/T)
       return
    else if(TRIM(sp_list) == 'C4H10') then
!  from Moses 1992
       VAPOR = 1333._wp*10._wp**(8.446_wp-1461.2_wp/T)
       return
    else if(TRIM(sp_list) == 'C6H6') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_list) == 'C7H8') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_list) == 'C8H10') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_list) == 'RinG') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_list) == 'H2O') then
!  from Moses 1992
!       VAPOR = 1333._wp*10._wp**(9.184_wp-0.2185*10999.398_wp/T)
       VAPOR = 1333.22368_wp*10._wp**(-2445.5646_wp/T+8.2312_wp*LOG10(T)-0.01677006_wp*T+1.20514E-5_wp*(T**2)-6.757169_wp)
       return
    else if(TRIM(sp_list) == 'CO2') then
!  from Moses 1992
       VAPOR = 1333._wp*EXP(2.13807649E+01_wp-2.57064700E+03_wp/T-7.78129489E+04_wp/T**2  &
            +4.32506256E+06_wp/T**3-1.20671368E+08_wp/T**4+1.34966306E+09_wp/T**5)
       return
    else if(TRIM(sp_list) == 'HCN') then
       VAPOR = 10._wp**(12.54747_wp-(1893.068_wp/(T+0.309_wp)))
       return
    else if(TRIM(sp_list) == 'HNC') then
       VAPOR = 10._wp**(12.54747_wp-(1893.068_wp/(T+0.309_wp)))
       return
    else if(TRIM(sp_list) == 'HC3N') then
       VAPOR =10._wp**(13.305_wp-(2210._wp/(T)))
       return
    else if(TRIM(sp_list) == 'HC5N') then
       VAPOR = 10._wp**(13.305_wp-(2210._wp/(T)))  ! assumed same as HC3N
       return
    else if(TRIM(sp_list) == 'CH3CN') then
       VAPOR = 10._wp**(10.52111_wp-(1492.375_wp/(T-24.208_wp)))
       return
    else if(TRIM(sp_list) == 'C3H3N') then
       VAPOR = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))
       return
    else if(TRIM(sp_list) == 'C4H3N') then
       VAPOR = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))  !  assumed same as C3H3N
       return
    else if(TRIM(sp_list) == 'C5H5N') then
       VAPOR = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))  !  assumed same as C3H3N
       return
    else if(TRIM(sp_list) == 'C2N2') then
       VAPOR = 10._wp**(12.53784_wp-(1566.647_wp/(T-10.461_wp)))
       return
    else if(TRIM(sp_list) == 'C4N2') then
       VAPOR = 10._wp**(14.73702_wp-(3722.003_wp/(T+3.036_wp)))
       return
    else
       VAPOR = 1.0E30_wp
    end if
  end function VAPOR


! This function calculates the differential cross-section for secondary electron production, with secondary electron energy Es, primary electron energy E and ionization potential I
  function DifCS_SEC(Es,E,I)
    use types, only: wp => dp
    implicit none
    real(wp) :: E, Es, E0, I, DifCS_SEC
    E0 = 13.0 ! shape parameter for N2
!    I = 15.6
    DifCS_SEC = ((1. + (Es/E0))**(-2.1))/(E0*tanh((E-I)/(2.*E0)))
  end function DifCS_SEC


! This function finds the index of the bin containing Etest in an ascending vector elctreV with bin size elctDeV.
  function find_bin(elctreV,elctDeV,Etest)
    use types, only: wp => dp
    implicit none
    real(wp), intent(in), dimension(:) :: elctreV, elctDeV
    real(wp), intent(in) :: Etest
    integer :: FinD_Bin
    integer :: i, nelb
    nelb = SIZE(elctreV)
    FinD_Bin = 1
    do i = 1, nelb-1
       if((elctreV(i)-0.5_wp*elctDeV(i) <= Etest) .and.   &
            (Etest <= elctreV(i+1)-0.5_wp*elctDeV(i+1)) ) then 
          FinD_Bin = i 
          EXIT
       endif
    enddo
    if(Etest > elctreV(nelb-1)+0.5_wp*elctDeV(nelb-1)) FinD_Bin=nelb
!    if(FinD_Bin == 0) then
!       WRITE(*,*) ' ERROR in FinD_Bin, SET TO 1'
!       FinD_Bin = 1
!    end if


!   do ne = 1, nelb
!       egrid(ne) = elctreV(ne)-half*elctDeV(ne)
!   end do
!   FinD_Bin = LOCATE(egrid,Etest)
!

  end function find_bin


end module utils
