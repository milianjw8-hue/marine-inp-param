#!/usr/bin/env bash
# SCAM6 (CESM2.1) M-PACE experiment matrix for marine INP sensitivity:
#   marine INP scheme {none, m18, w15, w15warm}  x  Hallett-Mossop rime splintering {on, off}  = 8 cases
#
# Requirements: a working CESM2.1 checkout (CAM tag cam_cesm2_1_rel, checked at commit 405a3f2) with
# the SCAM IOP data (atm/cam/scam/iop/MPACE_4scam.nc) under DIN_LOC_ROOT.
# The create_newcase options follow the CESM2.1 SCAM user-mods (scam_mpace); check them for your
# machine (MACH, compiler, project) before running.
#
# Usage:  CESMROOT=~/cesm2.1 CASEROOT=~/scam_cases MACH=mymachine bash scam6_mpace_matrix.sh [--submit]
set -euo pipefail

: "${CESMROOT:?set CESMROOT to the CESM2.1 checkout}"
: "${CASEROOT:?set CASEROOT to the folder for the cases}"
MACH="${MACH:-}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
CAMSRC="$CESMROOT/components/cam/src/physics/cam"
SUBMIT="${1:-}"

# warm-window parameters: replace with the values fitted to your observations (marinp fit)
WW_A="${WW_A:-6.5}"; WW_B="${WW_B:--0.12}"
# MOA profile (sensitivity choice, see README): surface MOA (ug m-3) and scale height (m)
MOA_SFC="${MOA_SFC:-0.1}"; MOA_H="${MOA_H:-1000.}"

mkdir -p "$CASEROOT"
for scheme in none m18 w15 w15warm; do
  for hm in on off; do
    case_name="mpace_${scheme}_hm${hm}"
    case_dir="$CASEROOT/$case_name"
    if [[ -d "$case_dir" ]]; then echo "skip existing $case_dir"; continue; fi
    mach_opt=(); [[ -n "$MACH" ]] && mach_opt=(--mach "$MACH")
    "$CESMROOT/cime/scripts/create_newcase" --case "$case_dir" --compset FSCAM --res T42_T42 \
        --user-mods-dir scam_mpace "${mach_opt[@]}" --run-unsupported
    pushd "$case_dir" >/dev/null

    # source modifications
    sm=SourceMods/src.cam
    cp "$REPO/fortran/marine_inp_mod.F90" "$REPO/cam6/SourceMods/src.cam/marine_inp_cam.F90" "$sm/"
    cp "$CAMSRC/hetfrz_classnuc_cam.F90" "$CAMSRC/micro_mg2_0.F90" "$sm/"
    (cd "$sm" && patch -p4 < "$REPO/cam6/patches/hetfrz_classnuc_cam.F90.patch" \
              && patch -p4 < "$REPO/cam6/patches/micro_mg2_0.F90.patch")

    cat >> user_nl_cam <<EOF
fincl1 = 'MARFRZIMM','MARINP','DSTFREZIMM','BCFREZIMM','FREQIMM','NUMICE','NUMLIQ','CLDICE','CLDLIQ',
         'AWNI','AWNC','TGCLDLWP','TGCLDIWP','CLOUD','MSACWIO','T','Q'
EOF
    ./case.setup
    disable_hm=.false.; [[ "$hm" == off ]] && disable_hm=.true.
    rundir="$(./xmlquery --value RUNDIR)"
    mkdir -p "$rundir"
    cat > "$rundir/marinp_in" <<EOF
&marinp_nl
 marinp_scheme = '$scheme'
 marinp_moa_source = 'profile'
 moa_sfc = $MOA_SFC
 moa_h = $MOA_H
 omoc = 2.0
 ww_a = $WW_A
 ww_b = $WW_B
 disable_hm = $disable_hm
/
EOF
    cp "$rundir/marinp_in" ./marinp_in.used   # keep a copy with the case
    ./case.build
    [[ "$SUBMIT" == "--submit" ]] && ./case.submit
    popd >/dev/null
    echo "ready: $case_dir"
  done
done
