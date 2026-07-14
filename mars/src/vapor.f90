real(wp) pure function vapor(sp, T)
! Calculates the saturation vapor pressure of chemical species sp at
! temperature T. 
! Condensation currently deactivated. 
!? Need to update for organics
  use types, only: wp => dp
  implicit none
  real(wp), intent(in) :: T
  character(len=*), intent(in) :: sp

  if (trim(sp) == 'H2O') then
    if (T > 273.16._wp) then
      ! Buck [1996] for liquid water
      vapor = 611.21 * exp((18.678 - (T + 273.15) / 234.5) &
        * ((T + 273.15) / (257.14 + T + 273.15)))
    else
      ! Murphy and Koop [2005] for H2O ice (valid down to 110 K)
      vapor = Pa_to_cgs * exp(9.550426 − 5723.265 / T &
        + 3.53068 * log(T) − 0.00728332 * T)
      return
    end if
  else if (trim(sp) == 'CO2') then
    ! Giauque & Egan [1937]
    vapor = ten * mmHg_to_cgs * ten ** (-1354.210 / T + 8.69903 + &
      0.0015880 * T - 4.5107E-6 * T**2)
!  from Moses 1992
!    vapor = 1333._wp*exp(2.13807649E+01_wp-2.57064700E+03_wp/T-7.78129489E+04_wp/T**2  &
!         +4.32506256E+06_wp/T**3-1.20671368E+08_wp/T**4+1.34966306E+09_wp/T**5)
    return
  else if (trim(sp) == 'C2H2') then
!  based on Tickner & Losing 1951 data
    vapor = 1.333_wp*10._wp**(9.25_wp-1201.75_wp/T)
    return
  else if (trim(sp) == 'C2H4') then
    if (T < 104._wp) then
      vapor = 1333._wp*10._wp**(8.724_wp - 901.6_wp/(T-2.555_wp))
    else if ((T >= 104._wp) .AND. (T < 120._wp)) then
      vapor = 1333._wp*10._wp**(50.79_wp - 1703._wp/T - 17.141_wp*log10(T))
    else if (T >= 120._wp) then
      vapor = 1333._wp*10._wp**(6.74756_wp-585._wp/(T-18.16_wp))
    end if
    return
  else if (trim(sp) == 'C2H6') then
    if (T < 90._wp) then
      vapor = 1333._wp*10._wp**(10.01_wp - 1085._wp/(T-0.561_wp))
    else if (T > 90._wp) then
      vapor = 1333._wp*10._wp**(5.9366_wp-1086.17_wp/T+3.83464_wp*log10(1000._wp/T))
    end if
    return
  else if (trim(sp) == 'C3H8') then
!  based on Tickner & Losing 1951 data
    vapor = 1333._wp*10._wp**(8.16173_wp-1176._wp/T)
    return
  else if (trim(sp) == 'C4H2') then
!  from Moses 1992
    vapor = 1333._wp*10._wp**(5.3817_wp-3300.5_wp/T + 16.63415_wp*log10(1000._wp/T))
    return
  else if (trim(sp) == 'C4H6') then
!  from Moses 1992
!       write(*,*) ' CALCULATE C4H10 VP, T = ',T
!       write(*,*) ' CALCULATE C4H10 VP, ARG = ',8.446_wp-1461.2_wp/T
    vapor = 1333._wp*10._wp**(8.032581_wp-1441.42_wp/T)
    return
  else if (trim(sp) == 'C4H10') then
!  from Moses 1992
    vapor = 1333._wp*10._wp**(8.446_wp-1461.2_wp/T)
    return
  else if (trim(sp) == 'C6H6') then
    vapor = 10._wp*exp(26._wp-7640._wp/(T+30._wp))
    return
  else if (trim(sp) == 'C7H8') then
    vapor = 10._wp*exp(26._wp-7640._wp/(T+30._wp))
    return
  else if (trim(sp) == 'C8H10') then
    vapor = 10._wp*exp(26._wp-7640._wp/(T+30._wp))
    return
  else if (trim(sp) == 'RING') then
    vapor = 10._wp*exp(26._wp-7640._wp/(T+30._wp))
    return

  else if (trim(sp) == 'HCN') then
    vapor = 10._wp**(12.54747_wp-(1893.068_wp/(T+0.309_wp)))
    return
  else if (trim(sp) == 'HNC') then
    vapor = 10._wp**(12.54747_wp-(1893.068_wp/(T+0.309_wp)))
    return
  else if (trim(sp) == 'HC3N') then
    vapor =10._wp**(13.305_wp-(2210._wp/(T)))
    return
  else if (trim(sp) == 'HC5N') then
    ! assumed same as HC3N
    vapor = 10._wp**(13.305_wp-(2210._wp/(T)))
    return
  else if (trim(sp) == 'CH3CN') then
    vapor = 10._wp**(10.52111_wp-(1492.375_wp/(T-24.208_wp)))
    return
  else if (trim(sp) == 'C3H3N') then
    vapor = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))
    return
  else if (trim(sp) == 'C4H3N') then
    ! assumed same as C3H3N
    vapor = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))
    return
  else if (trim(sp) == 'C5H5N') then
    ! assumed same as C3H3N
    vapor = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))
    return
  else if (trim(sp) == 'C2N2') then
    vapor = 10._wp**(12.53784_wp-(1566.647_wp/(T-10.461_wp)))
    return
  else if (trim(sp) == 'C4N2') then
    vapor = 10._wp**(14.73702_wp-(3722.003_wp/(T+3.036_wp)))
    return
  else
    vapor = 1.0E30_wp    ! 1E24 bar, i.e. no condensation
  end if
end function vapor