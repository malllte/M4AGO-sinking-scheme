module mo_m4ago_kind

use iso_fortran_env, only:  int32, int64, real32, real64, real128

#ifdef HAMOCC
  use mod_types, only: r8
#endif

  implicit none

  private

  public :: wp


#ifdef REAL64
  integer,parameter :: rk = real64
#elif REAL128
  integer,parameter :: rk = real128
#else
  integer,parameter :: rk = real32
#endif

!#ifdef INT64
!  integer, parameter :: ik = int64
!#else
!  integer, parameter :: ik = int32
!#endif




#ifdef HAMOCC
  integer, parameter :: wp = r8 !selected_real_kind()
#else
  integer, parameter :: wp = real64
#endif

end module mo_m4ago_kind
