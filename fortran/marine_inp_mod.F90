!-----------------------------------------------------------------------------------------------
! marine_inp_mod: marine ice-nucleating particle (INP) parameterizations for host models.
!
! Pure Fortran (no host-model dependencies). Same formulas as the Python package marinp
! (src/marinp/params.py and enhancement.py); tests/test_fortran.py checks that both agree.
!
!   W15  Wilson et al. (2015), Nature 525, 234-238, doi:10.1038/nature14986
!        INPs per gram TOC: n_m(T) = exp(11.2186 - 0.4459 T),  T in degC, -27 <= T <= -6.5
!   M18  McCluskey et al. (2018), JGR Atmos. 123, 6196-6212, doi:10.1029/2017JD028033
!        active site density per SSA surface: n_s(T) = exp(-0.545 (T_K - 273.15) + 1.0125) m-2
!   Warm window (this project): additive term in the same per-gram-TOC space as W15
!        n_add(T) = 10**(a + b (T - T_ref)) * W(T),  W = two logistic ramps (T_cold .. T_warm)
!        a and b must be fitted to observations (marinp.fit_warm_enhancement); no default is physical.
!
! Difference from the Python package: outside the nominal temperature range the functions return
! 0 (not NaN) unless extrapolate is .true., so a host model never receives NaN.
!
! Units: T in degC, mass concentrations in g m-3, surface area in m2 m-3, INP number in m-3.
!-----------------------------------------------------------------------------------------------
module marine_inp_mod

  implicit none
  private

  integer, parameter, public :: r8 = selected_real_kind(12)

  real(r8), parameter, public :: w15_a = 11.2186_r8, w15_b = 0.4459_r8
  real(r8), parameter, public :: w15_tmin = -27.0_r8, w15_tmax = -6.5_r8
  real(r8), parameter, public :: m18_slope = 0.545_r8, m18_intercept = 1.0125_r8
  real(r8), parameter, public :: m18_tmin = -30.0_r8, m18_tmax = -5.0_r8

  integer, parameter, public :: scheme_none = 0, scheme_m18 = 1, scheme_w15 = 2, scheme_w15warm = 3

  type, public :: warm_window_t
     real(r8) :: a = 7.0_r8          ! log10 g-1 at T_ref (placeholder until fitted)
     real(r8) :: b = -0.15_r8        ! log10 g-1 per degC
     real(r8) :: t_warm = -10.0_r8   ! degC
     real(r8) :: t_cold = -22.0_r8   ! degC
     real(r8) :: t_ref = -16.0_r8    ! degC
     real(r8) :: delta = 1.0_r8      ! degC, width of the logistic ramps
  end type warm_window_t

  public :: w15_nm, m18_ns, ww_window, ww_nm_add, nm_w15warm
  public :: inp_number, immersion_rate, moa_profile, scheme_from_name

contains

  elemental function w15_nm(t_c, extrapolate) result(nm)
    ! W15 cumulative INPs per gram TOC [g-1]
    real(r8), intent(in) :: t_c
    logical,  intent(in) :: extrapolate
    real(r8) :: nm
    if (extrapolate .or. (t_c >= w15_tmin .and. t_c <= w15_tmax)) then
       nm = exp(w15_a - w15_b*t_c)
    else
       nm = 0.0_r8
    end if
  end function w15_nm

  elemental function m18_ns(t_c, extrapolate) result(ns)
    ! M18 active site density of sea spray aerosol [m-2]
    real(r8), intent(in) :: t_c
    logical,  intent(in) :: extrapolate
    real(r8) :: ns
    if (extrapolate .or. (t_c >= m18_tmin .and. t_c <= m18_tmax)) then
       ns = exp(-m18_slope*t_c + m18_intercept)
    else
       ns = 0.0_r8
    end if
  end function m18_ns

  elemental function ww_window(t_c, ww) result(w)
    real(r8), intent(in) :: t_c
    type(warm_window_t), intent(in) :: ww
    real(r8) :: w
    w = 1.0_r8/(1.0_r8 + exp(-(t_c - ww%t_cold)/ww%delta)) * &
        1.0_r8/(1.0_r8 + exp( (t_c - ww%t_warm)/ww%delta))
  end function ww_window

  elemental function ww_nm_add(t_c, ww) result(nm)
    ! additive warm-window INPs per gram TOC [g-1]
    real(r8), intent(in) :: t_c
    type(warm_window_t), intent(in) :: ww
    real(r8) :: nm
    nm = 10.0_r8**(ww%a + ww%b*(t_c - ww%t_ref)) * ww_window(t_c, ww)
  end function ww_nm_add

  elemental function nm_w15warm(t_c, ww, extrapolate) result(nm)
    real(r8), intent(in) :: t_c
    type(warm_window_t), intent(in) :: ww
    logical,  intent(in) :: extrapolate
    real(r8) :: nm
    nm = w15_nm(t_c, extrapolate) + ww_nm_add(t_c, ww)
  end function nm_w15warm

  elemental function inp_number(scheme, t_c, toc_g_m3, ssa_m2_m3, ww, extrapolate) result(n)
    ! marine INP number concentration [m-3]
    !   W15 and W15+warm use TOC mass [g m-3]; M18 uses SSA surface area [m2 m-3]
    integer,  intent(in) :: scheme
    real(r8), intent(in) :: t_c, toc_g_m3, ssa_m2_m3
    type(warm_window_t), intent(in) :: ww
    logical,  intent(in) :: extrapolate
    real(r8) :: n
    select case (scheme)
    case (scheme_m18)
       n = m18_ns(t_c, extrapolate)*max(ssa_m2_m3, 0.0_r8)
    case (scheme_w15)
       n = w15_nm(t_c, extrapolate)*max(toc_g_m3, 0.0_r8)
    case (scheme_w15warm)
       n = nm_w15warm(t_c, ww, extrapolate)*max(toc_g_m3, 0.0_r8)
    case default
       n = 0.0_r8
    end select
  end function inp_number

  elemental function immersion_rate(n_inp, nc, ni, dt) result(rate)
    ! Deterministic immersion freezing rate [m-3 s-1] in the liquid-cloud volume.
    ! The INP spectrum is cumulative: at temperature T at most n_inp droplets can be frozen.
    ! Ice already present (ni) is taken to come from earlier freezing, so only the deficit freezes,
    ! limited by the droplet number nc, within one time step dt [s].
    ! All number concentrations in m-3 (in-cloud).
    real(r8), intent(in) :: n_inp, nc, ni, dt
    real(r8) :: rate
    rate = max(min(n_inp, nc) - max(ni, 0.0_r8), 0.0_r8)/dt
  end function immersion_rate

  elemental function moa_profile(z_m, moa_sfc, scale_height_m) result(moa)
    ! prescribed marine organic aerosol mass profile: moa_sfc * exp(-z/H)  (same unit as moa_sfc)
    real(r8), intent(in) :: z_m, moa_sfc, scale_height_m
    real(r8) :: moa
    moa = moa_sfc*exp(-max(z_m, 0.0_r8)/max(scale_height_m, 1.0_r8))
  end function moa_profile

  function scheme_from_name(name) result(scheme)
    character(len=*), intent(in) :: name
    integer :: scheme
    select case (trim(adjustl(name)))
    case ('m18', 'M18')
       scheme = scheme_m18
    case ('w15', 'W15')
       scheme = scheme_w15
    case ('w15warm', 'W15warm', 'W15WARM')
       scheme = scheme_w15warm
    case ('none', 'off', '')
       scheme = scheme_none
    case default
       scheme = -1
    end select
  end function scheme_from_name

end module marine_inp_mod
