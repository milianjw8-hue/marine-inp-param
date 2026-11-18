! test driver for cam6/SourceMods/src.cam/marine_inp_cam.F90 (with stubs.F90 in place of CAM modules)
program drv
  use marine_inp_cam
  implicit none
  real(8) :: r, n
  call marinp_cam_readnl()
  ! T = -15 degC, z = 200 m, sea salt 1e-9 / 2e-9 kg m-3, nc = 5e7 m-3, ni = 1e3 m-3, dt = 1800 s
  call marinp_frzimm(258.15d0, 200d0, 1d-9, 2d-9, 5d7, 1d3, 1800d0, r, n)
  print '(a,2es24.15,l3)', 'rate,ninp,hm_off=', r, n, marinp_disable_hm
end program
