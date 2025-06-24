module subs
  
  interface
    subroutine gauleg(x1, x2, x, w)
      use types, only: wp => dp
      real(wp), intent(in) :: x1, x2
      real(wp), dimension(:), intent(out) :: x, w
    end subroutine gauleg
  
    function difcs_rees(w,ep)
      use types, only: wp => dp
      implicit none
      real(wp) :: ep, w, term1, term2, difcs_rees
    end function difcs_rees
  
    subroutine wait
    end subroutine wait
  end interface
  
  
  contains
  
  integer pure function find_name(xname, sp_name)
  ! This functions finds the index of a species with name "xname" in the list 
  ! of all model species "sp_name". 
  ! Note: sp_name starts at index 0.
    implicit none
    character(len=12), intent(in) :: xname ! species to search for
    character(len=12), intent(in), dimension(0:) :: sp_name ! list of species
    integer :: n_sp ! size of sp_name
    integer :: i_sp ! loop variable

    n_sp = size(sp_name,1)-1 ! -1 because index starts at 0
    find_name = 0
    do i_sp = 0, n_sp
      if(trim(adjustl(xname)) == trim(adjustl(sp_name(i_sp)))) then
        find_name = i_sp
        exit
      end if
    end do
  end function find_name


  function locate(xtab, x)
  ! This function locates the value x in a sorted list xtab. The earlier index will be picked.
    use types, only: wp => dp
    implicit none
    real(wp), dimension(:), intent(in) :: xtab
    real(wp), intent(in) :: x
    integer :: locate
    integer :: n,jl,jm,ju
    logical :: ascnd ! is list ascending?
    n=size(xtab)
    ascnd = (xtab(n) >= xtab(1))
    jl=0
    ju=n+1
    do
      if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd .eqv. (x >= xtab(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do
    if (x == xtab(1)) then
      locate=1
    else if (x == xtab(n)) then
      locate=n-1 !?
    else
      locate=jl
    end if
  end function locate


  subroutine LUBKSB(a, indx, b)
    use types, only: wp => dp
    implicit none

    real(wp), intent(in), dimension(:,:) :: a
    integer, intent(in), dimension(:) :: indx
    real(wp), intent(inout), dimension(:) :: b

    integer :: i,ii,ll,j,n
    real(wp) :: sum

    n = SIZE(b)
    ii = 0
    do 12 i = 1, n
    ll = indx(i)
    sum = b(ll)
    b(ll) = b(i)
    if (ii .ne. 0) then
    do 11 j = ii, i - 1
    sum = sum - (a(i,j) * b(j))
 11 continue
    else if (sum .ne. 0._wp) then
    ii = i
    end if
    b(i) = sum
 12 continue
    do 14 i = n, 1, -1
    sum = b(i)
    if (i .lt. n) then
    do 13 j = i + 1, n
    sum = sum - (a(i,j) * b(j))
 13 continue
    end if
    b(i) = sum / a(i,i)
 14 continue
    return 
  end subroutine LUBKSB



  subroutine LUDCMP(a, indx, d)
    use types, only: wp => dp
    implicit none
    real(wp), intent(inout), dimension(:,:) :: a
    integer, intent(out) :: indx(:)
    real(wp), intent(out) :: d

    real(wp), parameter :: tiny = 1.E-60_wp
    integer i,j,k,imax,n
    real(wp) vv(SIZE(a,1)),aamax,sum,dum

    intrinsic ABS
    
    n = SIZE(a,1)
    d = 1._wp
    do 12 i = 1, n
    aamax = 0._wp
    do 11 j = 1, n
      if (abs(a(i,j)) > aamax) aamax = ABS(a(i,j))
 11  continue
     if (aamax == 0._wp) stop 'Singular matrix.'
     vv(i) = 1._wp / aamax
 12  continue
       do 19 j = 1, n
         if (j .gt. 1) then
            do 14 i = 1, j - 1
               sum = a(i,j)
               if (i > 1) then
                  do 13 k = 1, i - 1
                     sum = sum - (a(i,k) * a(k,j))
13                continue
                  a(i,j) = sum
               end if
14          continue
         end if
         aamax = 0._wp
         do 16 i = j, n
            sum = a(i,j)
            if (j .gt. 1) then
               do 15 k = 1, j - 1
                  sum = sum - (a(i,k) * a(k,j))
15             continue
               a(i,j) = sum
            end if
            dum = vv(i) * abs(sum)
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            end if
16       continue
         if (j .ne. imax) then
            do 17 k = 1, n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
17          continue
            d = - d
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (j .ne. n) then
            if (a(j,j) == 0._wp) a(j,j) = tiny
            dum = 1._wp / a(j,j)
            do 18 i = j + 1, n
               a(i,j) = a(i,j) * dum
18          continue
          end if
19     continue
       if (a(n,n) == 0._wp) a(n,n) = tiny
       return
  end subroutine LUDCMP
                         

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



  subroutine TRIDAG (a, b, c, r, u)
    use types, only: wp => dp
    implicit none

    real(wp), intent(in), dimension(:) ::  a, b, c, r
    real(wp), intent(out), dimension(:) ::  u
    integer :: n, j
    real(wp), dimension(SIZE(a)) :: gam
    real(wp) :: bet
    real(wp), parameter :: zero=0.

    n = SIZE(a)
    u = zero

    if (b(1) == zero) STOP 'tridag: error'
    bet=b(1)
    u(1)=r(1)/bet
    do j = 2, n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if (bet == zero) STOP 'tridag failed'
      u(j)=(r(j)-a(j)*u(j-1))/bet
    end do
    do j = n-1, 1, -1
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do

    return
  end subroutine TRIDAG


! This function interpolates linearly to find yval values at various xval values given the ascending list (xtab, ytab).
! Values beyond the range spanned by xtab are given the extreme ytab values.
  subroutine intrp(xtab,ytab,xval,yval)
    use types, only: wp => dp
    implicit none
    real(wp), intent(in), dimension(:) :: xtab, ytab, xval
    real(wp), intent(out), dimension(:) :: yval
    integer :: ntab, nval, nv, l1, l2
    ntab = SIZE(xtab)
    nval = SIZE(xval)
    do nv = 1, nval
      if (xval(nv) <= xtab(1)) then
            yval(nv) = ytab(1)
        else if (xval(nv) >= xtab(ntab)) then
            yval(nv) = ytab(ntab)
        else
           l1 = LOCATE(xtab,xval(nv))
           l2 = l1 + 1
           yval(nv)=ytab(l1)+(xval(nv)-xtab(l1))*(ytab(l2)-ytab(l1))/   &
                (xtab(l2)-xtab(l1))
        end if
     end do
     return
   end subroutine intrp


! Calculates the vapor pressure of sp_name at temperature T.

  function VAPOR( sp_name, T)
    use types, only: wp => dp
    implicit none
    real(wp), intent(in) :: T
    CHARACTER(LEN=*), intent(in) :: sp_name
    real(wp) :: VAPOR

    if(TRIM(sp_name) == 'C2H2') then
!  based on Tickner & Losing 1951 data
       VAPOR = 1.333_wp*10._wp**(9.25_wp-1201.75_wp/T)
       return
    else if(TRIM(sp_name) == 'C2H4') then
       if (T < 104._wp) then
       VAPOR = 1333._wp*10._wp**(8.724_wp - 901.6_wp/(T-2.555_wp))
       else if ((T >= 104._wp) .AND. (T < 120._wp)) then
       VAPOR = 1333._wp*10._wp**(50.79_wp - 1703._wp/T - 17.141_wp*LOG10(T))
       else if (T >= 120._wp) then
       VAPOR = 1333._wp*10._wp**(6.74756_wp-585._wp/(T-18.16_wp))
       end if
       return
    else if(TRIM(sp_name) == 'C2H6') then
       if (T < 90._wp) then
       VAPOR = 1333._wp*10._wp**(10.01_wp - 1085._wp/(T-0.561_wp))
       else if (T > 90.) then
       VAPOR = 1333._wp*10._wp**(5.9366_wp-1086.17_wp/T+3.83464_wp*LOG10(1000._wp/T))
       end if
       return
    else if(TRIM(sp_name) == 'C3H8') then
!  based on Tickner & Losing 1951 data
       VAPOR = 1333._wp*10._wp**(8.16173_wp-1176._wp/T)
       return
    else if(TRIM(sp_name) == 'C4H2') then
!  from Moses 1992
       VAPOR = 1333._wp*10._wp**(5.3817_wp-3300.5_wp/T + 16.63415_wp*LOG10(1000._wp/T))
       return
    else if(TRIM(sp_name) == 'C4H6') then
!  from Moses 1992
!       WRITE(*,*) ' CALCULATE C4H10 VP, T = ',T
!       WRITE(*,*) ' CALCULATE C4H10 VP, ARG = ',8.446_wp-1461.2_wp/T
       VAPOR = 1333._wp*10._wp**(8.032581_wp-1441.42_wp/T)
       return
    else if(TRIM(sp_name) == 'C4H10') then
!  from Moses 1992
       VAPOR = 1333._wp*10._wp**(8.446_wp-1461.2_wp/T)
       return
    else if(TRIM(sp_name) == 'C6H6') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_name) == 'C7H8') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_name) == 'C8H10') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_name) == 'RinG') then
       VAPOR = 10._wp*EXP(26._wp-7640._wp/(T+30._wp))
       return
    else if(TRIM(sp_name) == 'H2O') then
!  from Moses 1992
!       VAPOR = 1333._wp*10._wp**(9.184_wp-0.2185*10999.398_wp/T)
       VAPOR = 1333.22368_wp*10._wp**(-2445.5646_wp/T+8.2312_wp*LOG10(T)-0.01677006_wp*T+1.20514E-5_wp*(T**2)-6.757169_wp)
       return
    else if(TRIM(sp_name) == 'CO2') then
!  from Moses 1992
       VAPOR = 1333._wp*EXP(2.13807649E+01_wp-2.57064700E+03_wp/T-7.78129489E+04_wp/T**2  &
            +4.32506256E+06_wp/T**3-1.20671368E+08_wp/T**4+1.34966306E+09_wp/T**5)
       return
    else if(TRIM(sp_name) == 'HCN') then
       VAPOR = 10._wp**(12.54747_wp-(1893.068_wp/(T+0.309_wp)))
       return
    else if(TRIM(sp_name) == 'HNC') then
       VAPOR = 10._wp**(12.54747_wp-(1893.068_wp/(T+0.309_wp)))
       return
    else if(TRIM(sp_name) == 'HC3N') then
       VAPOR =10._wp**(13.305_wp-(2210._wp/(T)))
       return
    else if(TRIM(sp_name) == 'HC5N') then
       VAPOR = 10._wp**(13.305_wp-(2210._wp/(T)))  ! assumed same as HC3N
       return
    else if(TRIM(sp_name) == 'CH3CN') then
       VAPOR = 10._wp**(10.52111_wp-(1492.375_wp/(T-24.208_wp)))
       return
    else if(TRIM(sp_name) == 'C3H3N') then
       VAPOR = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))
       return
    else if(TRIM(sp_name) == 'C4H3N') then
       VAPOR = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))  !  assumed same as C3H3N
       return
    else if(TRIM(sp_name) == 'C5H5N') then
       VAPOR = 10._wp**(8.9178_wp-(706.474_wp/(T-109.392_wp)))  !  assumed same as C3H3N
       return
    else if(TRIM(sp_name) == 'C2N2') then
       VAPOR = 10._wp**(12.53784_wp-(1566.647_wp/(T-10.461_wp)))
       return
    else if(TRIM(sp_name) == 'C4N2') then
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


end module subs
