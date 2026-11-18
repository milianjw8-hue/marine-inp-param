! Writes a table for tests/test_fortran.py:
!   T  w15_nm  m18_ns  ww_window  ww_nm_add  nm_w15warm  rate_test
! Warm-window parameters are read from stdin (a b t_warm t_cold t_ref delta).
program test_driver
  use marine_inp_mod
  implicit none
  type(warm_window_t) :: ww
  real(r8) :: t, rate
  integer :: i
  read(*, *) ww%a, ww%b, ww%t_warm, ww%t_cold, ww%t_ref, ww%delta
  do i = 0, 60
     t = -30.0_r8 + 0.5_r8*i
     rate = immersion_rate(inp_number(scheme_w15warm, t, 1.0e-7_r8, 0.0_r8, ww, .false.), &
                           5.0e7_r8, 1.0e3_r8, 1800.0_r8)
     write(*, '(f8.2, 6es24.15)') t, w15_nm(t, .false.), m18_ns(t, .false.), ww_window(t, ww), &
          ww_nm_add(t, ww), nm_w15warm(t, ww, .false.), rate
  end do
  if (scheme_from_name('w15warm') /= scheme_w15warm .or. scheme_from_name('bogus') /= -1) stop 2
end program test_driver
