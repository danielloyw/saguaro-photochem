! Defines precision for reals and other types
module types
  use iso_fortran_env, only: real32, real64

  implicit none
   
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
   
end module