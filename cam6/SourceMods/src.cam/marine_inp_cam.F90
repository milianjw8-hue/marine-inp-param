!-----------------------------------------------------------------------------------------------
! marine_inp_cam: CAM6 (CESM2.1) interface for marine_inp_mod.
!
! Adds a deterministic marine INP immersion-freezing rate to the CAM6 heterogeneous freezing
! (hetfrz_classnuc_cam.F90, which passes FRZIMM to MG2). Dust and black carbon stay on the
! default classical-nucleation-theory scheme; the marine term is added on top.
!
! Settings are read from the file 'marinp_in' in the run directory (namelist /marinp_nl/), so no
! change to CAM build-namelist is needed. Without the file the scheme is off.
!
! Marine organic aerosol (MOA) source, because CESM2.1 MAM4 has no marine organic species:
!   moa_source = 'profile'  MOA(z) = moa_sfc_ugm3 * exp(-z / moa_scale_height_m)
!   moa_source = 'seasalt'  MOA    = seasalt_om_ratio * (sea-salt mass, accumulation [+ coarse])
! TOC = MOA / om_oc.  M18 uses the sea-salt surface area from the sea-salt mass and an assumed
! volume-to-surface (Sauter) diameter ssa_dvs_um.
! All of these are sensitivity choices, not established values: report them with the results.
!
! Written for cam_cesm2_1_rel (commit 405a3f2). Not compiled inside CAM by the author: build and
! check in your SCAM6 case before use.
!-----------------------------------------------------------------------------------------------
module marine_inp_cam

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun
  use marine_inp_mod,  only: warm_window_t, inp_number, immersion_rate, moa_profile, scheme_from_name, &
                             scheme_none, scheme_m18, mr8 => r8

  implicit none
  private

  public :: marinp_cam_readnl, marinp_active, marinp_frzimm
  logical, public :: marinp_disable_hm = .false.   ! read by the micro_mg2_0.F90 patch (SIP experiment)

  integer  :: scheme = scheme_none
  logical  :: use_coarse_seasalt = .false.
  logical  :: extrapolate = .false.
  character(len=16) :: moa_source = 'profile'
  real(r8) :: moa_sfc_ugm3 = 0.1_r8
  real(r8) :: moa_scale_height_m = 1000.0_r8
  real(r8) :: seasalt_om_ratio = 0.1_r8
  real(r8) :: om_oc = 2.0_r8
  real(r8) :: ssa_dvs_um = 1.0_r8
  real(r8), parameter :: rho_seasalt = 1900.0_r8  ! kg m-3
  type(warm_window_t) :: ww

contains

  subroutine marinp_cam_readnl()
    use mpishorthand
    character(len=16) :: marinp_scheme = 'none'
    character(len=16) :: marinp_moa_source = 'profile'
    real(r8) :: moa_sfc = 0.1_r8, moa_h = 1000.0_r8, ss_om = 0.1_r8, omoc = 2.0_r8, dvs = 1.0_r8
    real(r8) :: ww_a = 7.0_r8, ww_b = -0.15_r8, ww_t_warm = -10.0_r8, ww_t_cold = -22.0_r8
    real(r8) :: ww_t_ref = -16.0_r8, ww_delta = 1.0_r8
    logical  :: disable_hm = .false., coarse = .false., extrap = .false.
    integer  :: unitn, ierr, ibuf(4)
    real(r8) :: rbuf(11)
    logical  :: exists
    namelist /marinp_nl/ marinp_scheme, marinp_moa_source, moa_sfc, moa_h, ss_om, omoc, dvs, &
         ww_a, ww_b, ww_t_warm, ww_t_cold, ww_t_ref, ww_delta, disable_hm, coarse, extrap

    if (masterproc) then
       inquire(file='marinp_in', exist=exists)
       if (exists) then
          open(newunit=unitn, file='marinp_in', status='old', action='read')
          read(unitn, marinp_nl, iostat=ierr)
          close(unitn)
          if (ierr /= 0) call endrun('marinp_cam_readnl: error reading namelist marinp_nl from marinp_in')
       end if
       ibuf = (/ scheme_from_name(marinp_scheme), merge(1, 0, disable_hm), merge(1, 0, coarse), &
                 merge(1, 0, extrap) /)
       if (trim(marinp_moa_source) /= 'profile' .and. trim(marinp_moa_source) /= 'seasalt') ibuf(1) = -1
       if (ibuf(1) < 0) call endrun('marinp_cam_readnl: unknown marinp_scheme or marinp_moa_source')
       rbuf = (/ moa_sfc, moa_h, ss_om, omoc, dvs, ww_a, ww_b, ww_t_warm, ww_t_cold, ww_t_ref, ww_delta /)
    end if
#ifdef SPMD
    call mpibcast(ibuf, 4, mpiint, 0, mpicom)
    call mpibcast(rbuf, 11, mpir8, 0, mpicom)
    call mpibcast(marinp_moa_source, len(marinp_moa_source), mpichar, 0, mpicom)
#endif
    scheme = ibuf(1)
    marinp_disable_hm = ibuf(2) == 1
    use_coarse_seasalt = ibuf(3) == 1
    extrapolate = ibuf(4) == 1
    moa_source = marinp_moa_source
    moa_sfc_ugm3 = rbuf(1); moa_scale_height_m = rbuf(2); seasalt_om_ratio = rbuf(3)
    om_oc = rbuf(4); ssa_dvs_um = rbuf(5)
    ww%a = rbuf(6); ww%b = rbuf(7); ww%t_warm = rbuf(8); ww%t_cold = rbuf(9)
    ww%t_ref = rbuf(10); ww%delta = rbuf(11)

    if (masterproc) then
       write(iulog, *) 'marinp: scheme=', scheme, ' moa_source=', trim(moa_source), &
            ' moa_sfc_ugm3=', moa_sfc_ugm3, ' H=', moa_scale_height_m, ' ss_om=', seasalt_om_ratio, &
            ' om_oc=', om_oc, ' ssa_dvs_um=', ssa_dvs_um, ' disable_hm=', marinp_disable_hm
       write(iulog, *) 'marinp: warm window a,b,Twarm,Tcold,Tref,delta=', ww
    end if
  end subroutine marinp_cam_readnl

  logical function marinp_active()
    marinp_active = scheme /= scheme_none
  end function marinp_active

  subroutine marinp_frzimm(t_k, zm, seasalt_ac_kgm3, seasalt_co_kgm3, nc_m3, ni_m3, dt, rate_m3s, ninp_m3)
    ! t_k       temperature [K]
    ! zm        height above surface [m]
    ! seasalt_* sea-salt mass concentration, accumulation and coarse mode (interstitial + cloud-borne) [kg m-3]
    ! nc_m3     in-cloud droplet number [m-3];  ni_m3  in-cloud ice number [m-3];  dt  time step [s]
    ! rate_m3s  marine immersion freezing rate [m-3 s-1];  ninp_m3  marine INP number at T [m-3]
    real(r8), intent(in)  :: t_k, zm, seasalt_ac_kgm3, seasalt_co_kgm3, nc_m3, ni_m3, dt
    real(r8), intent(out) :: rate_m3s, ninp_m3
    real(r8) :: t_c, ss, moa_g_m3, toc_g_m3, ssa_m2_m3

    rate_m3s = 0.0_r8
    ninp_m3 = 0.0_r8
    if (scheme == scheme_none) return
    t_c = t_k - 273.15_r8
    ss = seasalt_ac_kgm3
    if (use_coarse_seasalt) ss = ss + seasalt_co_kgm3
    if (trim(moa_source) == 'seasalt') then
       moa_g_m3 = seasalt_om_ratio*ss*1.0e3_r8
    else
       moa_g_m3 = moa_profile(zm, moa_sfc_ugm3, moa_scale_height_m)*1.0e-6_r8
    end if
    toc_g_m3 = moa_g_m3/om_oc
    ssa_m2_m3 = 0.0_r8
    if (scheme == scheme_m18) ssa_m2_m3 = 6.0_r8*ss/(rho_seasalt*ssa_dvs_um*1.0e-6_r8)
    ninp_m3 = inp_number(scheme, real(t_c, mr8), real(toc_g_m3, mr8), real(ssa_m2_m3, mr8), ww, extrapolate)
    rate_m3s = immersion_rate(real(ninp_m3, mr8), real(nc_m3, mr8), real(ni_m3, mr8), real(dt, mr8))
  end subroutine marinp_frzimm

end module marine_inp_cam
