# Cytosolic calcium dynamics

## L-type Ca current (ICa & ICaK)[^Cortassa2006]

Common pool of subspace calcium model

$$
\begin{aligned}
\alpha &= 0.4 e^{(V_m+2) / 10}  \\
\beta &= 0.4 e^{-(V_m+2) / 13}  \\
\alpha^\prime  &=  a \alpha \\
\beta^\prime  &=  \beta / b \\
\gamma &= 0.1875 [Ca^{2+}]_{ss}  \\
C_0 &= 1 - C_0 - C_1 - C_2 - C_3 - C_4 - O - C_{Ca0} - C_{Ca1} - C_{Ca2} - C_{Ca3} - C_{Ca4}   \\
v_{01} &= 4\alpha C_0 - \beta C_1   \\
v_{12} &= 3\alpha C_1 - 2\beta C_2   \\
v_{23} &= 2\alpha C_2 - 3\beta C_3   \\
v_{34} &= \alpha C_3 - 4\beta C_4   \\
v_{45} &= f C_4 - g O   \\
v_{67} &= 4\alpha^\prime C_{Ca0} - \beta^\prime C_{Ca1}   \\
v_{78} &= 3\alpha^\prime C_{Ca1} - 2\beta^\prime C_{Ca2}   \\
v_{89} &= 2\alpha^\prime C_{Ca2} - 3\beta^\prime C_{Ca3}   \\
v_{910} &= \alpha^\prime C_{Ca3} - 4\beta^\prime C_{Ca4}   \\
v_{06} &= \gamma C_0 - \omega C_{Ca0}  \\
v_{17} &= a \gamma C_1 - \omega C_{Ca1} / b \\
v_{28} &= a^2 \gamma C_2 - \omega C_{Ca2} / b^2  \\
v_{39} &= a^3 \gamma C_3 - \omega C_{Ca3} / b^3  \\
v_{410} &= a^4 \gamma C_4 - \omega C_{Ca4} / b^4  \\
\frac{dC_0}{dt}  &=  -v_{01} -v_{06}  \\
\frac{dC_1}{dt}  &=  v_{01} - v_{12} - v_{17}  \\
\frac{dC_2}{dt}  &=  v_{12} - v_{23} - v_{28}  \\
\frac{dC_3}{dt}  &=  v_{23} - v_{34} - v_{34}  \\
\frac{dC_4}{dt}  &=  v_{34} - v_{45} - v_{410} \\
\frac{dO}{dt}  &=  v_{45}  \\
\frac{dC_{Ca0}}{dt}  &=  v_{06} - v_{67}  \\
\frac{dC_{Ca1}}{dt}  &=  v_{17} + v_{67} - v_{78}  \\
\frac{dC_{Ca2}}{dt}  &=  v_{28} + v_{78} - v_{89}  \\
\frac{dC_{Ca3}}{dt}  &=  v_{39} + v_{89} - v_{910}  \\
I_{Ca}^{max} &= \Phi_{Ca}(P_{Ca}, z_{Ca}, V_m, 0.001, 0.341[Ca^{2+}]_o)  \\
I_{Ca} &= 6 I_{Ca}^{max}  \cdot y_{Ca}  \cdot O  \\
I_{Ca,K} &= y_{Ca}  \cdot O  \cdot  \Phi_{Ca}(P_{K}, z_{K}, V_m, [K^+]_i, [K^+]_o)  \\
P_{K}  &= P_{K}^{max}  \cdot Hill(I_{Ca}^{half}, I_{Ca}^{max}, 1)  \\
y_\infty &= \frac{1}{1 + e^{(V_m + 55) / 7.5}} + \frac{0.5}{1 + e^{(-V_m + 21) / 6}}  \\
\tau_y &= 20 + \frac{600}{1 + e^{(V_m + 30) / 9.5}}  \\
\frac{dy_{Ca}}{dt}  &=  \frac{y_\infty - y_{Ca}}{\tau_y}  \\
\end{aligned}
$$

| Parameter      | Value                 | Units            | Description                                |
| -------------- | --------------------- | ---------------- | ------------------------------------------ |
| $A$            | $2$                   |                  | Mode transition parameter                  |
| $B$            | $2$                   |                  | Mode transition parameter                  |
| $\omega$       | $10$                  | Hz               | Mode transition parameter                  |
| $f$            | $300$                 | Hz               | Transition rate into open state            |
| $g$            | $2000$                | Hz               | Transition rate into open state            |
| $f^\prime$     | $0$                   | Hz               | Transition rate into open state            |
| $g^\prime$     | $0$                   | Hz               | Transition rate into open state            |
| $P_{Ca}^{LCC}$ | $1.24 \cdot 10^{-3}$  | $cm/s$           | L-type Ca2+ channel permeability to Ca2+   |
| $P_{K}^{LCC}$  | $1.11 \cdot 10^{-11}$ | $cm/s$           | L-type Ca2+ channel permeability to K+     |
| $I_{Ca, half}$ | $-0.4583$             | $\mu A \ cm^{2}$ | ICa level that reduces equation Pk by half |

## Ryanodine receptor (calcium release, Jrel)[^Cortassa2006]

With optimization from Plank et al. (2008)

$$
\begin{aligned}
P_{C1} &= 1 - P_{O1} - P_{O2} - P_{C2}  \\
\\
If \ [Ca^{2+}]_{ss} &\ge [Ca^{2+}]_{ss}^{*} : \\
P_{O1} &:= (P_{O1} + P_{C1}) Hill(k_a^+[Ca^{2+}]_{ss}^n, k_a^-, 1)  \\
v_{o1c1} &= 0  \\
If \ [Ca^{2+}]_{ss} &< [Ca^{2+}]_{ss}^{*} : \\
v_{o1c1} &= -k_a^-P_{O1} + k_a^+[Ca^{2+}]_{ss}^n P_{C1} \\
\\
v_{o1o2} &= k_b^+ [Ca^{2+}]_{ss}^m P_{O1} - k_b^- P_{O2} \\
v_{o1c2} &= k_c^+ P_{O1} - k_c^- P_{C2} \\
\dot{P_{O1}}  &= -v_{o1c1} - v_{o1o2} - v_{o1c2}  \\
\dot{P_{O2}}  &= v_{o1o2}  \\
\dot{P_{C2}}  &= v_{o1c2}  \\
J_{rel} &= r_{ryr} (P_{O1} + P_{O2})([Ca^{2+}]_{JSR} - [Ca^{2+}]_{ss}) \\
\end{aligned}
$$

| Parameter | Value                  | Units               | Description               |
| --------- | ---------------------- | ------------------- | ------------------------- |
| $v_1$     | $3600$                 | Hz                  | RyR flux channel constant |
| $n$       | $4$                    |                     | Cooperativity parameter   |
| $m$       | $3$                    |                     | Cooperativity parameter   |
| $k_a^+$   | $1.215 \cdot 10^{13} $ | $\text{Hz mM}^{-4}$ | RyR rate constant         |
| $k_a^-$   | $576$                  | $\text{Hz}$         | RyR rate constant         |
| $k_b^+$   | $4.05 \cdot 10^{6} $   | $\text{Hz mM}^{-3}$ | RyR rate constant         |
| $k_b^-$   | $1930$                 | $\text{Hz}$         | RyR rate constant         |
| $k_c^+$   | $100$                  | $\text{Hz}$         | RyR rate constant         |
| $k_c^-$   | $0.8$                  | $\text{Hz}$         | RyR rate constant         |

## Plama membrane calcium ATPase (PMCA) current (IpCa)[^Cortassa2006]

Modified rate expression incorporating the ATP-dependence of pump activity.
Plama membrane calcium ATPase (PMCA) rate exhibits two different K0.5 values for ATP

$$
\begin{aligned}
f_{ATP} &= Hill([ATP]_i  \cdot Hill(K_{i,ADP}^{PMCA}, [ADP]_i, 1), K_{M1,ATP}^{PMCA}, 1) + Hill([ATP]_i, K_{M2,ATP}^{PMCA}, 1)  \\
f_{Ca} &= Hill([Ca^{2+}]_i, \ K_{M, Ca}^{PMCA}, 1) \\
I_{pCa} &= I_{max}^{PMCA}  \cdot f_{Ca}  \cdot f_{ATP} \\
\end{aligned}
$$

| Parameter         | Value   | Units                 | Description                                                   |
| ----------------- | ------- | --------------------- | ------------------------------------------------------------- |
| $I_{max}^{PMCA}$  | $0.575$ | $\mu A \cdot cm^{-2}$ | Maximum sarcolemmal Ca2+ pump current                         |
| $K_{Ca}^{PMCA}$   | $0.5$   | $\mu M$               | Ca2+ half-saturation constant for sarcolemmal Ca2+ pump       |
| $K_{ATP1}^{PMCA}$ | $0.012$ | $mM$                  | First ATP half-saturation constant for sarcolemmal Ca2+ pump  |
| $K_{ATP2}^{PMCA}$ | $0.23$  | $mM$                  | Second ATP half-saturation constant for sarcolemmal Ca2+ pump |
| $K_{ADP}^{PMCA}$  | $1.0$   | $mM$                  | ADP inhibition constant for sarcolemmal Ca2+ pump             |

## SERCA calcium pump (Jup)[^Cortassa2006]

Michaelis-Menten dependence of enzyme activity with respect to ATP and mixed-type inhibition of the enzyme by ADP.

Reversible during diastole (low cytoplasmic calcium).

$$
\begin{aligned}
J_{up} &= \frac{V_{f}^{up}f_b-V_{r}^{up}r_b}{(1 + f_b + r_b)f_{ATP}^{SERCA}} \\
f_b &= \left( \frac{[Ca^{2+}]_i}{K_{fb}} \right)^{N_{fb}} \\
r_b &= \left( \frac{[Ca^{2+}]_{NSR}}{K_{rb}} \right)^{N_{rb}} \\
f_{ATP}^{SERCA} &= K_{m,up}^{ATP} / ([ATP]_i  \cdot Hill(K_{i1, up}, [ADP]_i, 1)) + Hill(K_{i2, up}, [ADP]_i, )^{-1}  \\
\end{aligned}
$$

[^Cortassa2006]: Cortassa S, Aon MA, O'Rourke B, et al. A computational model integrating electrophysiology, contraction, and mitochondrial bioenergetics in the ventricular myocyte. Biophys J. 2006;91(4):1564-89. [PMC1518641](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1518641/)

| Parameter            | Value     | Units   | Description                                    |
| -------------------- | --------- | ------- | ---------------------------------------------- |
| $V_{max, f}^{SERCA}$ | $0.2989$  | $mM/s$  | SERCA forward rate parameter                   |
| $V_{max, b}^{SERCA}$ | $0.3179$  | $mM/s$  | SERCA reverse  rate parameter                  |
| $K_{f}^{SERCA}$      | $0.24$    | $\mu M$ | Forward Ca2+ half-saturation constant of SERCA |
| $K_{r}^{SERCA}$      | $1.64269$ | $mM$    | Reverse Ca2+ half-saturation constant of SERCA |
| $N_{f}^{SERCA}$      | $1.4$     |         | Forward cooperativity constant of SERCA        |
| $N_{r}^{SERCA}$      | $1.0$     |         | Reverse  cooperativity constant of SERCA       |
| $K_{ATP}^{SERCA}$    | $0.01$    | $mM$    | ATP half-saturation constant for SERCA         |
| $K_{ADP1}^{SERCA}$   | $0.14$    | $mM$    | ADP first inhibition constant for SERCA        |
| $K_{ADP2}^{SERCA}$   | $5.1$     | $mM$    | ADP second  inhibition constant for SERCA      |

## Ca2+ transport and buffering parameters

| Symbol          | Value    | Units               | Description                                          |
| --------------- | -------- | ------------------- | ---------------------------------------------------- |
| $\tau_{tr}$     | $574.7$  | Hz                  | Time constant for transfer from subspace to myoplasm |
| $\tau_{xfer}$   | $9090$   | Hz                  | Time constant for transfer from NSR to JSR           |
| $K_{m}^{CMDN}$  | $2.38$   | $\mu M$             | Ca2+ half saturation constant for calmodulin         |
| $K_{m}^{CSQN}$  | $0.8$    | $mM$                | Ca2+ half saturation constant for calsequestrin      |
| $h_{trpn}^{+}$  | $100000$ | $\text{Hz mM}^{-1}$ | Ca2+ on-rate for troponin high-affinity sites        |
| $h_{trpn}^{-}$  | $0.33$   | $\text{Hz}$         | Ca2+ off-rate for troponin high-affinity sites       |
| $l_{trpn}^{+}$  | $100000$ | $\text{Hz mM}^{-1}$ | Ca2+ on-rate for troponin low-affinity sites         |
| $l_{trpn}^{-}$  | $40$     | $\text{Hz}$         | Ca2+ off-rate for troponin low-affinity sites        |
| $\Sigma[HTRPN]$ | $0.14$   | $mM$                | Total troponin high-affinity sites                   |
| $\Sigma[LTRPN]$ | $0.07$   | $mM$                | Total troponin low-affinity sites                    |
| $\Sigma[CMDN]$  | $0.05$   | $mM$                | Total myoplasmic calmodulin concentration            |
| $\Sigma[CQSN]$  | $15$     | $mM$                | Total NSR calsequestrin concentration                |

## ODE for cytosolic calcium

$$
\begin{aligned}
β_i &= Hill((K_m^{CMDN} + [Ca^{2+}]_i)^2, K_m^{CMDN}  \cdot  [CMDN]_{tot}, 1)  \\
β_{SR} &= Hill((K_m^{CSQN} + [Ca^{2+}]_{SR})^2, K_m^{CSQN}  \cdot  [CSQN]_{tot}, 1)  \\
\frac{d[Ca^{2+}]_i}{dt} &= \beta_i(J_{xfer}\frac{V_{ss}}{V_{myo}} - J_{up} - J_{trpn} - (I_{Ca,b} -2I_{NaCa} + I_{pCa})\frac{A_{cap}}{2V_{myo}F} + (V_{NaCa} - V_{uni})\frac{V_{mito}}{V_{myo}}) \\
\frac{d[Ca^{2+}]_{SR}}{dt} &= \beta_{SR}(J_{up}\frac{V_{myo}}{V_{SR}} - J_{rel}\frac{V_{ss}}{V_{SR}}) \\
\end{aligned}
$$