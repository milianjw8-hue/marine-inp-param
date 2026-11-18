# CAM6 / SCAM6 연결 (CESM2.1)

해양 INP 파라미터화(W15, M18, W15 + 고온 가산항)를 CESM2.1 CAM6의 이질 결빙에 넣고,
SCAM6 M-PACE 사례로 민감도 실험을 하는 데 필요한 파일입니다.

> **확인 필요:** 이 연결부는 CAM 원본 소스(`cam_cesm2_1_rel`, 커밋 405a3f2)에 패치가 깨끗이 적용되는지,
> 그리고 대역 모듈(stub)과 함께 컴파일되는지까지만 확인했습니다. 실제 CESM 빌드와 SCAM 실행은 사용하는 기계에서 확인하세요.

## 구성

| 파일 | 내용 |
|---|---|
| `../fortran/marine_inp_mod.F90` | 순수 Fortran 물리 모듈. CAM에 의존하지 않음. Python 구현과 같은 값인지 시험함 |
| `SourceMods/src.cam/marine_inp_cam.F90` | CAM 연결 모듈. 설정 읽기, MOA 결정, 결빙률 계산 |
| `patches/hetfrz_classnuc_cam.F90.patch` | 이질 결빙 루틴에 해양 INP 결빙률을 더하고 출력 변수 2개 추가 |
| `patches/micro_mg2_0.F90.patch` | Hallett-Mossop 상고대 파편(rime splintering)을 끌 수 있게 함 |
| `scam6_mpace_matrix.sh` | M-PACE 실험 8개(INP 방식 4 × Hallett-Mossop 켬/끔)를 만들고 빌드 |

## 어떻게 동작하나

1. CAM6 기본 이질 결빙(먼지·블랙카본, 고전 핵생성 이론)은 그대로 둡니다.
2. −38 ~ −4 °C의 액체 구름에서 해양 INP 수 n_INP(T)를 계산합니다.
   - W15, W15 + 고온 가산항: TOC 질량 × 그램당 INP 수
   - M18: 해염 표면적 × 활성점 밀도
3. 결빙률은 결정론적으로 정합니다.
   **결빙률 = max(min(n_INP, N_c) − N_i, 0) / Δt**
   - N_c: 구름 안 물방울 수
   - N_i: 이미 있는 얼음 수(격자 평균 / 액체 구름 비율로 근사)
   - INP 스펙트럼은 누적값이므로 같은 온도에서 같은 INP가 거듭 얼지 않게 합니다.
4. 이 값을 CAM의 `FRZIMM`(cm⁻³ s⁻¹)에 더해 MG2 미세물리로 넘깁니다.

### MOA를 정하는 방법

CESM2.1 기본 MAM4에는 해양 유기 에어로졸(MOA) 성분이 없습니다. 그래서 두 방법 중 하나를 고릅니다.

| `marinp_moa_source` | MOA | 설정 |
|---|---|---|
| `profile` | 지표값 × exp(−z/H) | `moa_sfc` (µg m⁻³), `moa_h` (m) |
| `seasalt` | 해염 질량 × 비율 | `ss_om`, `coarse` (조대 모드 포함 여부) |

TOC = MOA / `omoc`. 이 값들은 모두 **민감도 선택값**입니다. 관측(예: 선상 에어로졸 화학 관측)으로 정하고, 결과와 함께 보고하세요.

## 설정 파일 `marinp_in`

실행 폴더(RUNDIR)에 두면 읽습니다. 파일이 없으면 해양 INP 항은 꺼집니다. 그래서 CAM `build-namelist`를 고칠 필요가 없습니다.

```fortran
&marinp_nl
 marinp_scheme = 'w15warm'      ! none | m18 | w15 | w15warm
 marinp_moa_source = 'profile'  ! profile | seasalt
 moa_sfc = 0.1                  ! ug m-3 (profile)
 moa_h = 1000.                  ! m (profile)
 ss_om = 0.1                    ! MOA / sea salt (seasalt)
 coarse = .false.               ! include coarse-mode sea salt (seasalt)
 omoc = 2.0                     ! OM:OC
 dvs = 1.0                      ! um, Sauter diameter for the M18 surface area
 ww_a = 6.5                     ! warm window, fitted with marinp.fit_warm_enhancement
 ww_b = -0.12
 ww_t_warm = -10., ww_t_cold = -22., ww_t_ref = -16., ww_delta = 1.
 disable_hm = .false.           ! .true.: Hallett-Mossop off (SIP experiment)
 extrap = .false.               ! .true.: use the formulas outside their temperature range
/
```

## 실험 설계 (M-PACE, 2004년 10월, 북극 혼합상 구름)

| 요인 | 수준 |
|---|---|
| 해양 INP | none (대조), m18, w15, w15warm |
| Hallett-Mossop 상고대 파편 | on (기본), off |

이렇게 8개 실험이 됩니다. Hallett-Mossop을 요인으로 넣는 이유가 있습니다. SCAM6 M-PACE 연구에서 이차빙정생성(SIP)이 얼음 생성의 대부분을 차지했고,
강한 빙정핵 생성이 SIP를 줄이는 상호작용도 보고되었기 때문입니다(Zhao et al., 2021; ACP 22, 2585, 2022).
CESM2.1 공개판의 SIP는 Hallett-Mossop 하나뿐입니다. 위 연구들에서 쓴 다른 SIP 과정(빗방울 파쇄, 얼음끼리 충돌 파쇄)은 공개판에 없습니다(확인 필요).

```bash
CESMROOT=~/cesm2.1 CASEROOT=~/scam_cases MACH=<기계 이름> \
  WW_A=6.5 WW_B=-0.12 MOA_SFC=0.1 MOA_H=1000. \
  bash cam6/scam6_mpace_matrix.sh            # --submit 을 붙이면 제출까지
```

## 결과 비교

```python
from marinp.scam import compare, find_cases
t = compare(find_cases("~/scam_cases"))
t.to_csv("mpace_summary.csv", index=False)
```

사례마다 액체·얼음 수경로, 액체 비율, 얼음 수 최대값, 해양·먼지·BC 결빙률의 연직 적분을 냅니다.
같은 Hallett-Mossop 설정의 대조 실험(`none`)과의 차이(`d_*_vs_none`)로 INP 효과를 봅니다.
Hallett-Mossop 끔과 켬의 차이(`d_*_hm_off_minus_on`)로 SIP 효과를 봅니다.
관측(ARM NSA 원격탐사, M-PACE 항공 관측)과의 비교는 관측 자료를 따로 준비해 하세요.

## 주의할 점

- 결빙률은 단순 결정론 방식입니다. INP의 예측 변수(소모, 수송)는 없습니다. 같은 기단에서 INP가 다 쓰였다는 효과는 N_i로만 근사합니다.
- `state%zm`은 지표 위 높이로 씁니다. MOA 연직 분포의 기준 높이가 다르면 바꾸세요.
- M18 표면적은 해염 질량과 가정한 Sauter 지름으로 계산하므로 지름 선택에 민감합니다.
- 논문에 쓰기 전에 원래 CAM과 비트 단위로 같은지 확인하세요. `marinp_in`이 없거나 `marinp_scheme='none'`이고 `disable_hm=.false.`이면 결과가 원래와 같아야 합니다.

## 참고

- Zhao, X., Liu, X., Phillips, V. T. J., Patade, S. (2021). Impacts of secondary ice production on Arctic mixed-phase clouds based on ARM observations and CAM6 single-column model simulations. *ACP*, 21, 5685–5703. https://doi.org/10.5194/acp-21-5685-2021
- Zhao, X., Liu, X., Burrows, S. M., Shi, Y. (2021). Effects of marine organic aerosols as sources of immersion-mode ice-nucleating particles on high-latitude mixed-phase clouds. *ACP*, 21, 2305–2327. https://doi.org/10.5194/acp-21-2305-2021
- ACP 22, 2585 (2022): M-PACE 단일층 혼합상 구름에서 빙정핵 생성과 SIP의 상호작용. https://acp.copernicus.org/articles/22/2585/2022/
- McCluskey, C. (2023). CESM AMWG 발표: D15 + M18 결정론적 침적 결빙을 CAM6/MG3에 적용. https://www.cesm.ucar.edu/sites/default/files/2023-03/2023-AMWG-C-McCluskey.pdf
