## Force generation[^Rice2000]

With ATP consumption from Cortassa, 2006[^Cortassa2006]
tTe rate of ATP hydrolysis associated with force generation through actomyosin ATPase depends explicitly on both ATP and ADP, as previously demonstrated.

$$
\begin{aligned}
f_{01} &= 3f_{XB} \\
f_{12} &= 10f_{XB} \\
f_{23} &= 7f_{XB} \\
g_{01} &= g_{XB}^{min} \\
g_{12} &= 2g_{XB}^{min} \\
g_{23} &= 3g_{XB}^{min} \\
g_{01,SL} &= \phi  \cdot g_{01} \\
g_{12,SL} &= \phi  \cdot g_{12} \\
g_{23,SL} &= \phi  \cdot g_{23} \\
g_{01,SL, off} &= \phi  \cdot g_{off} \\
\phi &= 1 + \frac{2.3-SL}{(2.3-1.7)^{1.6}} \\
K_{Ca}^{trop} &= \frac{k^-_{ltrpn}}{k^+_{ltrpn}} \\
K_{1/2}^{trop} &= \left( 1 + \frac{K_{Ca}^{trop}}{1.7  \cdot 10^{-3} - 0.8  \cdot 10^{-3}\frac{(SL-1.7)}{0.6}} \right)^{-1} \\
N_{trop} &= 3.5 \cdot SL - 2.0 \\
k_{np}^{trop} &= k_{pn}^{trop} \left( \frac{[LTRPNCa]}{K_{1/2}^{trop}[LTRPN]_{tot}} \right) ^{N_{trop}} \\
\Sigma PATHS &= g_{01}g_{12}g_{23} + f_{01}g_{12}g_{23} + f_{01}f_{12}g_{23} + f_{01}f_{12}f_{23} \\
P1_{max} &= \frac{f_{01}g_{12}g_{23}}{\Sigma PATHS} \\
P2_{max} &= \frac{f_{01}f_{12}g_{23}}{\Sigma PATHS} \\
P3_{max} &= \frac{f_{01}f_{12}f_{23}}{\Sigma PATHS} \\
Force &= \zeta \frac{[P_1] + 2[P_2] + 3[P_3] + [N_1]}{P1_{max} + 2P2_{max} + 3P3_{max}} \\
Force_{norm} &= \frac{[P_1] + [P_2] + [P_3] + [N_1]}{P1_{max} + P2_{max} + P3_{max}} \\
v_{01} &= f_{01} [P_0] - g_{01(SL)} [P_1] \\
v_{12} &= f_{12} [P_1] - g_{21(SL)} [P_2]  \\
v_{23} &= f_{23} [P_2] - g_{23(SL)} [P_3] \\
v_{04} &= k_{pn}^{trop} [P_0] - k_{np}^{trop} [N_0] \\
 [N_0] &= 1 - [P_0] - [P_1] - [P_2] - [P_3] - [N_1]  \\
v_{15} &= k_{pn}^{trop} [P_1] - k_{np}^{trop} [N_1] \\
v_{54} &= g_{01,off} [N_1]  \\
 [HTRPN] &=  [HTRPN]_{tot} - [HTRPNCa]  \\
 [LTRPN] &=  [LTRPN]_{tot} - [LTRPNCa]  \\
f_{ATP}^{AM} &= Hill([ATP]_i  \cdot Hill(K_{i,AM}^{ADP}, [ADP]_i, 1), K_{m,AM}^{ATP}, 1) \\
V_{AM} &= V_{max}^{AM}  \cdot f_{ATP}^{AM}  \cdot  \frac{f_{01}[P_0] + f_{12}[P_1] + f_{23}[P_2]}{f_{01} + f_{12} + f_{23}} \\
J_{trpn} &= \frac{d[HTRPNCa]}{dt} + \frac{d[LTRPNCa]}{dt} \\
\frac{d[HTRPNCa]}{dt} &= k^{+}_{htrpn}[Ca^{2+}]_i[HTRPN] - k^{-}_{htrpn}[HTRPNCa]  \\
\frac{d[LTRPNCa]}{dt} &= k^{+}_{ltrpn}[Ca^{2+}]_i[LTRPN] - k^{-}_{ltrpn}(1-\frac{2}{3}Force_{norm})[LTRPNCa]  \\
\frac{d[P_0]}{dt} &= - v_{01} - v_{04}  \\
\frac{d[P_1]}{dt} &= v_{01} - v_{12} - v_{15}  \\
\frac{d[P_2]}{dt} &= v_{12} - v_{23}  \\
\frac{d[P_3]}{dt} &= v_{23}  \\
\frac{d[N_1]}{dt} &= v_{15} - v_{54} \\
\end{aligned}
$$


[^Rice2000]: Rice JJ, Jafri MS, Winslow RL. Modeling short-term interval-force relations in cardiac muscle. Am J Physiol Heart Circ Physiol. 2000 Mar;278(3):H913-31. [APS](https://www.physiology.org/doi/full/10.1152/ajpheart.2000.278.3.H913)

[^Cortassa2006]: Cortassa S, Aon MA, O'Rourke B, et al. A computational model integrating electrophysiology, contraction, and mitochondrial bioenergetics in the ventricular myocyte. Biophys J. 2006;91(4):1564-89. [PMC1518641](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1518641/)

## Force generation parameters

| Symbol          | Value  | Units              | Description                                                  |
| --------------- | ------ | ------------------ | ------------------------------------------------------------ |
| $k_{pn}^{trop}$ | $40$   | $\text{Hz}$        | Transition rate from tropomyosin permissive to nonpermissive |
| $\text{SL}$     | $2.15$ | $\mu \text{m}$     | Sarcomere length                                             |
| $f_{XB}$        | $50$   | $\text{Hz}$        | Transition rate from weak to strong crossbridge              |
| $g_{XB}^{min}$  | $100$  | $\text{Hz}$        | Minimum transition rate from strong to weak crossbridge      |
| $\zeta$         | $0.1$  | $\text{N mm}^{-2}$ | Conversion factor normalizing to physiological force         |
| $V_{AM}^{max}$  | $7.2$  | $\text{mM s}^{-1}$ | Conversion factor normalizing to physiological force         |
| $K_{ATP}^{AM}$  | $0.03$ | $\text{mM}$        | ATP half-saturation constant of AM ATPase                    |
| $K_{ADP}^{AM}$  | $0.26$ | $\text{mM}$        | ADP inhibition constant of AM ATPase                         |
