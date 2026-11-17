# marine-inp-param (marinp)

해양 생물기원 빙정핵(Ice Nucleating Particle, INP)을 계산하는 Python 도구입니다.
문헌 파라미터화 두 가지와, 관측으로 맞추는 고온 생물기원 가산항, 물리 제약 머신러닝 에뮬레이터를 함께 제공합니다.

*English summary below.*

## 무엇을 하나

| 모듈 | 내용 |
|---|---|
| `marinp.params` | Wilson et al. (2015, W15) 총유기탄소(TOC) 기반 식, McCluskey et al. (2018, M18) 해양 비말 에어로졸(SSA) 표면적 기반 식 |
| `marinp.enhancement` | -10 ~ -22 °C 구간에서 W15를 넘는 INP를 더해 주는 **고온 생물기원 가산항**과, 관측 스펙트럼에서 그 계수를 맞추는 함수 |
| `marinp.emulator` | 온도와 에어로졸·해양 변수로 INP를 예측하는 **물리 제약 ML 에뮬레이터** (온도가 낮을수록 INP가 줄지 않도록 단조 제약) |
| `marinp.spectra` | 액적 동결 실험의 동결 비율을 대기 INP 농도로 바꾸는 Vali (1971) 변환 |

## 왜 필요한가

- 북극해 등 해양 관측에서 -10 ~ -22 °C의 비교적 따뜻한 온도에서 기존 파라미터화보다 INP가 많이 나오는 사례가 보고됩니다.
  열에 약한 생물기원 INP가 원인으로 꼽힙니다.
- 기후모델에 넣으려면 이 초과분을 **기존 식과 분리된 항**으로 표현하고, 관측으로 계수를 정하고, 불확실성을 함께 내야 합니다.
- 변수가 많아지면 회귀식 대신 ML이 낫지만, 물리적으로 말이 안 되는 예측(따뜻할수록 INP 증가)을 막아야 합니다.

## 설치

```bash
git clone https://github.com/milianjw8-hue/marine-inp-param.git
cd marine-inp-param
pip install -e ".[dev]"
```

## 바로 해 보기 (합성 예제)

```bash
python scripts/make_synthetic_data.py                     # data/synthetic_spectra.csv (가짜 자료)
python scripts/run_demo.py --data data/synthetic_spectra.csv \
    --features T_C log10_toc log10_chl --out outputs/demo
```

결과(`outputs/demo/`):

| 파일 | 내용 |
|---|---|
| `warm_term_fit.json` | 가산항 계수 a, b와 표준오차 |
| `spectrum_per_toc.png` | TOC 1 g당 INP 스펙트럼: 관측, W15, W15 + 가산항 |
| `cv_scores.json`, `parity.png` | 표본 단위 교차검증 점수(log10 공간)와 W15 대비 비교 |

합성 자료는 알려진 계수(a = 8.5, b = -0.12)로 만들었습니다. 맞춘 값이 이를 되찾는지 확인하는 용도이며 실제 관측 결과가 아닙니다.

## 내 관측에 적용

CSV 한 파일에 다음 열을 넣습니다. 한 행은 한 표본의 한 온도입니다.

| 열 | 단위 | 설명 |
|---|---|---|
| `sample_id` | - | 필터(표본) 식별자. 교차검증을 표본 단위로 나눌 때 씀 |
| `T_C` | °C | 동결 온도 |
| `inp_L` | L⁻¹ | INP 농도 |
| `toc_g_m3` | g C m⁻³ | 같은 시료의 총유기탄소 질량농도 (MOA만 있으면 `moa_to_toc`로 변환, OM:OC 기본 2) |
| 그 밖 | - | 예측 변수로 쓸 열(엽록소, 해수면온도, 풍속 등) |

```bash
python scripts/run_demo.py --data my_spectra.csv --features T_C log10_toc log10_chl
```

```python
from marinp import WarmWindow, fit_warm_enhancement, inp_enhanced
ww = fit_warm_enhancement(df.T_C, df.inp_L, df.toc_g_m3, template=WarmWindow(T_warm=-10, T_cold=-22))
inp = inp_enhanced(-15.0, 2e-7, ww)    # L^-1
```

## 주의할 점

- W15는 해수표면미세층(SML) 시료로 만든 식입니다. 대기 중 SSA에 대해서는 과대 추정 경향이 보고돼 있습니다.
  가산항을 해석할 때 기준선의 이 편향을 함께 고려해야 합니다.
- W15 유효 범위는 약 -27 ~ -6.5 °C로 두었고, 밖에서는 기본적으로 NaN을 돌려줍니다(`extrapolate=True`로 해제).
- M18의 유효 온도 범위는 원문으로 확인한 뒤 쓰세요. 코드의 기본 범위(-30 ~ -5 °C)는 넓게 잡은 값입니다.
- INP 농도는 측정 조건 그대로의 L⁻¹입니다. 표준상태 환산은 하지 않습니다.
- 이 저장소에는 **코드와 합성 자료만** 있습니다. 실제 관측 자료는 해당 논문의 자료 공개 경로를 따릅니다.

## 방법

식, 단위, 피팅과 검증 방법은 [`docs/METHODS.md`](docs/METHODS.md)에 있습니다.

## 시험

```bash
pytest
```

## 라이선스와 인용

MIT 라이선스입니다. 인용 정보는 `CITATION.cff`에 있습니다.

---

## English summary

`marinp` provides (1) the Wilson et al. (2015) TOC-based and McCluskey et al. (2018)
SSA-surface-based marine INP parameterizations, (2) an additive warm-temperature
(-10 to -22 °C) biogenic term on top of W15 with least-squares fitting in log10 space and
standard errors, (3) a physics-constrained gradient-boosting emulator (monotone in
temperature) in direct or hybrid residual-to-W15 mode with sample-grouped cross-validation,
and (4) Vali (1971) conversion of droplet-freezing frozen fractions to INP per litre of air.
The repository contains code and synthetic data only.
