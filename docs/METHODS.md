# 방법 (Methods)

## 1. 기본 파라미터화

### W15 — Wilson et al. (2015)

해수표면미세층 시료의 침적 동결(immersion freezing) 실험에서 얻은, 총유기탄소(TOC) 1 g당 누적 INP 수입니다.

$$n_{m}(T) = \exp(11.2186 - 0.4459\,T) \quad [\mathrm{g^{-1}}],\; T\,[^{\circ}\mathrm{C}]$$

$$N_{\mathrm{INP}}(T) = n_m(T)\, M_{\mathrm{TOC}} \times 10^{-3} \quad [\mathrm{L^{-1}}]$$

$M_{\mathrm{TOC}}$는 g C m⁻³입니다. 해양 유기에어로졸 질량(MOA)만 있으면 $M_{\mathrm{TOC}} = M_{\mathrm{MOA}} / (\mathrm{OM{:}OC})$로 바꿉니다.
OM:OC는 기본 2이며, 문헌에서는 1.2 ~ 2.6 범위로 민감도를 봅니다.
전구모델 적용 문헌을 따라 유효 범위를 -27 ~ -6.5 °C로 둡니다.

### M18 — McCluskey et al. (2018)

청정 북대서양 기단의 관측에서 얻은, SSA 표면적당 빙정핵 활성점 밀도입니다. Niemand et al. (2012)와 같은 함수형입니다.

$$n_s(T) = \exp[-0.545\,(T_K - 273.15) + 1.0125] \quad [\mathrm{m^{-2}}]$$

$$N_{\mathrm{INP}} = n_s(T)\, S_{\mathrm{SSA}} \times 10^{-3} \quad [\mathrm{L^{-1}}],\; S_{\mathrm{SSA}}\,[\mathrm{m^2\,m^{-3}}]$$

## 2. 고온 생물기원 가산항

W15와 같은 "TOC 1 g당" 공간에서 더하는 항으로 정의합니다. 같은 유기물 질량에 일관되게 작용하게 하기 위해서입니다.

$$n_{\mathrm{tot}}(T) = n_m^{\mathrm{W15}}(T) + 10^{\,a + b\,(T - T_{\mathrm{ref}})}\, W(T)$$

$$W(T) = \sigma\!\left(\frac{T - T_{\mathrm{cold}}}{\delta}\right)\, \sigma\!\left(\frac{T_{\mathrm{warm}} - T}{\delta}\right),\quad \sigma(x) = \frac{1}{1+e^{-x}}$$

- 기본값: $T_{\mathrm{warm}} = -10$ °C, $T_{\mathrm{cold}} = -22$ °C, $T_{\mathrm{ref}} = -16$ °C, $\delta = 1$ °C
- 곱셈 보정이 아니라 가산항이므로, 창 밖에서는 W15로 돌아가고 창 안에서만 별도 소스로 작동합니다.
- 계수 $a$, $b$는 관측으로 정합니다. 목적함수는 log10 공간 잔차입니다.

$$\min_{a,b} \sum_i \left[\log_{10} n_{\mathrm{tot}}(T_i) - \log_{10}\frac{N_i \times 10^3}{M_{\mathrm{TOC},i}}\right]^2$$

창 경계에서 $2\delta$ 안쪽까지의 양수 관측만 씁니다. Levenberg–Marquardt로 풀고, 최적점의 야코비안으로 표준오차를 구합니다.
같은 필터의 온도별 값은 서로 독립이 아니므로, 논문에 쓸 때는 표본 단위 부트스트랩으로 신뢰구간을 다시 구하는 것을 권합니다.

## 3. 물리 제약 ML 에뮬레이터

- 모델: 히스토그램 기반 그래디언트 부스팅 회귀(scikit-learn `HistGradientBoostingRegressor`)
- 입력: 온도 `T_C`, `log10_toc`, 선택 변수(예: `log10_chl`)
- 단조 제약: 온도에 대해 감소(차가울수록 INP가 줄지 않음). `direct` 모드에서는 TOC에 대해 증가.
- 목표 변수
  - `direct`: $\log_{10} N_{\mathrm{INP}}$
  - `residual`: $\log_{10}(N_{\mathrm{obs}} / N_{\mathrm{W15}})$. 예측은 $N_{\mathrm{W15}} \times 10^{\hat{r}}$ (물리식 + ML 보정의 하이브리드)
- 검증: 표본(필터) 단위 GroupKFold. 같은 필터의 다른 온도 값이 학습과 시험에 동시에 들어가지 않게 합니다.

## 4. 검증 지표 (log10 공간)

INP는 대략 로그정규 분포를 따르므로 모든 비교를 log10 공간에서 합니다. $d_i = \log_{10}\hat{N}_i - \log_{10}N_i$.

| 지표 | 정의 |
|---|---|
| `log10_bias` | $\overline{d}$ |
| `log10_rmse` | $\sqrt{\overline{d^2}}$ |
| `within_x2`, `within_x10` | $\lvert d_i \rvert \le \log_{10}2$, $\le 1$인 비율 |
| `r_log10` | log10 값의 피어슨 상관 |

## 5. 액적 동결 실험 변환

Vali (1971)의 특이 기술(singular description)을 따릅니다.

$$K(T) = -\frac{\ln[1 - f(T)]}{V_{\mathrm{drop}}},\qquad N_{\mathrm{INP}}(T) = K(T)\, D\, \frac{V_{\mathrm{liq}}}{V_{\mathrm{air}}}$$

$f$는 누적 동결 비율, $D$는 희석배수입니다. $f = 1$이면 값이 정의되지 않으므로 NaN(하한)으로 둡니다.

## 6. 수치모델 연결: 결정론적 침적 결빙률

Fortran 모듈(`fortran/marine_inp_mod.F90`)은 1~3절과 같은 식을 씁니다. 단, 유효 온도 범위 밖에서는 NaN 대신 0을 돌려줍니다.
모델 안의 결빙률은 다음과 같습니다.

$$R = \frac{\max\left(\min(n_{\mathrm{INP}}(T),\, N_c) - N_i,\ 0\right)}{\Delta t}$$

$n_{\mathrm{INP}}(T)$는 누적 INP 수, $N_c$는 구름 안 물방울 수, $N_i$는 이미 있는 얼음 수, $\Delta t$는 시간 간격입니다(모두 구름 안, m⁻³).
INP 스펙트럼이 누적값이므로, 이미 언 만큼은 다시 세지 않습니다. INP 소모와 수송을 예측하지 않는 단순한 방식이며, CAM6 연결 방법은 `cam6/README.md`에 있습니다.

## 참고 문헌

- Wilson, T. W., et al. (2015). A marine biogenic source of atmospheric ice-nucleating particles. *Nature*, 525, 234–238. https://doi.org/10.1038/nature14986
- McCluskey, C. S., et al. (2018). Marine and terrestrial organic ice-nucleating particles in pristine marine to continentally influenced northeast Atlantic air masses. *J. Geophys. Res. Atmos.*, 123, 6196–6212. https://doi.org/10.1029/2017JD028033
- McCluskey, C. S., DeMott, P. J., Ma, P.-L., & Burrows, S. M. (2019). Numerical representations of marine ice-nucleating particles in remote marine environments evaluated against observations. *Geophys. Res. Lett.* https://doi.org/10.1029/2018GL081861
- Niemand, M., et al. (2012). A particle-surface-area-based parameterization of immersion freezing on desert dust particles. *J. Atmos. Sci.*, 69, 3077–3092.
- Vali, G. (1971). Quantitative evaluation of experimental results on the heterogeneous freezing nucleation of supercooled liquids. *J. Atmos. Sci.*, 28, 402–409.
