# Changelog

## 0.2.0 (2026-11-18)

- `fortran/marine_inp_mod.F90`: W15, M18, warm-window term, deterministic immersion freezing rate and
  MOA profile in pure Fortran; checked against the Python package (relative difference < 1e-12)
- `cam6/`: CESM2.1 CAM6 interface (`marine_inp_cam.F90`, namelist file `marinp_in`), patches for
  `hetfrz_classnuc_cam.F90` (marine term added to FRZIMM, outputs MARFRZIMM and MARINP) and
  `micro_mg2_0.F90` (Hallett-Mossop switch), SCAM6 M-PACE 8-case script
- `marinp.scam`: comparison of SCAM6 cases against the control and the Hallett-Mossop twin
