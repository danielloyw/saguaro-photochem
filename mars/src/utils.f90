module utils
  
contains
  
integer pure function find_name(sp_name, sp_list)
! This function finds the index of a species with name "sp_name" in the list
! of all model species "sp_list". 
! Note: sp_list starts at index 0.
  implicit none
  character(len=12), intent(in) :: sp_name
  character(len=12), intent(in), dimension(0:) :: sp_list
  integer :: n_sp    ! size of sp_list
  integer :: i_sp    ! loop variable

  n_sp = size(sp_list,1)-1    ! -1 because index starts at 0
  find_name = -1
  do i_sp = 0, n_sp
    if (trim(adjustl(sp_name)) == trim(adjustl(sp_list(i_sp)))) then
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
  logical :: is_increasing    ! is list ascending?
  integer :: low, high, mid    ! index bounds of search window

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
  do while (high - low > 1)    ! low and high not adjacent indices
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


integer pure function find_bin(enrg, e_enrg, e_denrg)
! This function finds the index of the bin containing the value enrg in an
! monotically increasing vector e_enrg with bin sizes e_denrg.
! Returns extreme bins if outside of the range spanned by e_enrg. 
  use types, only: wp => dp
  use constants
  implicit none
  real(wp), intent(in) :: enrg
  real(wp), intent(in), dimension(:) :: e_enrg, e_denrg
  integer :: i, n_e_enrg
  
  n_e_enrg = size(e_enrg)
  
  ! initialize result to first bin
  find_bin = 1
  
  do i = 2, n_e_enrg-1
    if ((e_enrg(i) - half * e_denrg(i) <= enrg) .and. &
      (enrg <= e_enrg(i+1) - half * e_denrg(i+1))) then
      find_bin = i
      return
    end if
  end do
  
  ! case for last bin
  if (enrg > e_enrg(n_e_enrg) - half * e_denrg(n_e_enrg)) then
    find_bin = n_e_enrg
  end if
  
end function find_bin

subroutine heapsort(x_list, sort_order)
! This subroutine performs an indirect heapsort on x_list and returns
! sort_order, which contains the indices of the elements of x_list in ascending
! order. 
! Note: Starting index of x_list is 1. 
  use types, only: wp => dp
  implicit none
  real(wp), intent(in), dimension(:) :: x_list
  integer, intent(out), dimension(:) :: sort_order
  ! length of x_list
  integer :: n
  ! first and last indices spanned by heap in sort_order
  integer :: i_start, i_end
  ! swap value
  real(wp) :: t_value
  ! index of swap value in x_list
  integer :: t_index
  ! indices of active parent and child to be tested against swap value
  integer :: i_parent, i_child
  ! loop variables
  integer :: i
  
  ! initialize variables
  n = size(x_list)
  do concurrent (i = 1:n)
    sort_order(i) = i
  end do
  i_start = (n / 2) + 1    ! lowest parent in heap
  i_end = n
  
  do
    if (i_start > 1) then
    ! initial heap construction, progressing upwards from lowest parent
      i_start = i_start - 1
      t_index = sort_order(i_start)
      t_value = x_list(t_index)
    else
      ! take out last value for swap
      t_index = sort_order(i_end)
      t_value = x_list(t_index)
      ! move first and largest element to the back out of heap
      sort_order(i_end) = sort_order(1)
      i_end = i_end - 1
      if (i_end == 1) then    ! heap size 1 => all elements sorted
        sort_order(1) = t_index
        return 
      end if
    end if
    
    ! Sifting down algorithm
    i_parent = i_start
    i_child = 2 * i_parent    ! left child

    do while (i_child <= i_end)
      ! 1. Select the larger child
      if (i_child < i_end) then    ! if right child exists
        if (x_list(sort_order(i_child)) < x_list(sort_order(i_child+1))) then
          i_child = i_child + 1
        end if
      end if
      
      ! 2. If larger child > parent, promote to parent and proceed down heap.
      if (t_value < x_list(sort_order(i_child))) then
        sort_order(i_parent) = sort_order(i_child)
        i_parent = i_child
        i_child = 2 * i_parent
      else
        exit
      end if
    end do
    
    ! 3. Insert swap value at final position in current sift down
    sort_order(i_parent) = t_index
  end do
end subroutine heapsort


end module
