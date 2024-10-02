## Phosphate carrier[^Wei2011]

Follows equilibrium random Bi:Bi reaction kinetics 

$$
\begin{aligned}
J_{PiC} &= c_{PiC}(k_f AB - k_b PQ) / \Delta  \\
A &= [HPO_4^{2-}]_i / K_{pi,i} \\
P &= [HPO_4^{2-}]_m / K_{pi,m} \\
Q &= [OH^-]_i / K_{OH,i} \\
B &= [OH^-]_m / K_{OH,m} \\
\Delta &= 1 + A + B + P + Q + AB + PQ \\
k_b &= \frac{k_f K_{pi,m} K_{OH,i}}{K_{eq} K_{pi,i} K_{OH,m}}
\end{aligned}
$$

| Parameter  | Value   | Unit | Desc.                                     |
| ---------- | ------- | ---- | ----------------------------------------- |
| $K_{pi,i}$ | 11.06   | mM   | Extra-matrix Pi binding constant          |
| $K_{pi,m}$ | 11.06   | mM   | Mitochondrial matrix Pi binding constant  |
| $K_{OH,i}$ | 4.08E-5 | mM   | Extra-matrix OH- binding constant         |
| $K_{OH,m}$ | 4.08E-5 | mM   | Mitochondrial matrix OH- binding constant |
| $k_f$      | 1.5     | Hz   | Forward rate                              |
| $c_{PiC}$  | 4.9     | mM   | PiC activity                              |
| $K_{eq}$   | 1       | -    | Equilibrium constant of PiC               |

## Mitochondrial Sodium-hydrogen exchanger (mNHE)[^Wei2011]

Following Smith and Crampin's model of counterpart on the plasma membrane

$$
\begin{aligned}
J_{NHE} &= c_{NHE} f_h \frac{β^+_1β^+_2-β^-_1β^-_2}{β^+_1 + β^+_2 + β^-_1 + β^-_2}\\
f_h &= \frac{([H^+]_i)^n}{([H^+]_i)^n + K_i^n} \\
A &= [Na^+]_m / K_{Na}  \\
B &= [H^+]_i / K_{H}  \\
P &= [Na^+]_i / K_{Na}  \\
Q &= [H^+]_m / K_{H}  \\
β^+_1 &= \frac{k_1^+ A}{1 + A + Q}  \\
β^-_1 &= \frac{k^-_1 P}{1 + P + B}   \\
β^+_2 &= \frac{k_4^+ B}{1 + P + B}  \\
β^-_2 &= \frac{k^-_4 Q}{1 + A + Q}   \\
k_4^- &= \frac{k_1^+ k_4^+}{K_{eq} k_1^-}
\end{aligned}
$$

| Parameter   | Value    | Unit | Desc.                           |
| ----------- | -------- | ---- | ------------------------------- |
| $c_{NHE}$   | 0.00785  | mM   | NHE concentration               |
| $K_{Na}$    | 24       | mM   | Na Dissociation constant        |
| $K_H$       | 158.5E-6 | mM   | H Dissociation constant         |
| $K_i$       | 3.02E-6  | mM   | Proton binding constant         |
| $n$         | 3        | -    | Hill coefficient for H+ binding |
| $k_{1}^{+}$ | 25.2     | Hz   | NHE forward rate constant       |
| $k_{1}^{-}$ | 42.9     | Hz   | NHE backward rate constant      |
| $k_{4}^{+}$ | 160      | Hz   | NHE forward rate constant       |
| $K_{eq}$    | 1        | -    | Equilibrium constant of NHE     |

## Adenine Nucleotide translocator (ANT) [^Wei2011]

$$
\begin{aligned}
J_{ANT} &= V_{max}\frac{AB - \delta PQ}{(B + \delta^h P)(A + Q)}  \\
A &= [ATP^{4-}]_m  \\
B &= [ADP^{3-}]_i  \\
P &= [ATP^{4-}]_i  \\
Q &= [ADP^{3-}]_m  \\
\delta &= \text{exp}(-\Delta\Psi_m / V_T)
\end{aligned}
$$


| Parameter | Value | Unit    | Desc.            |
| --------- | ----- | ------- | ---------------- |
| $V_{max}$ | 3150  | mM * Hz | Maximal rate     |
| $h$       | 0.5   | -       | Fraction of dpsi |

## Mitochondrial calcium uniporter (MCU)[^Wei2011]

$$
\begin{aligned}
J_{uni} &= V_{max} \frac{S (1+S)^3}{(1+S)^4 + L(1 + A)^n} \frac{\delta}{e^\delta-1}  \\
S &= [Ca^{2+}]_i / K_{trans}  \\
A &= [Ca^{2+}]_i / K_{act}    \\
\delta &= -Z_{Ca} (\Delta\Psi_m - \Delta\Psi_0) / V_T
\end{aligned}
$$

| Parameter      | Value  | Unit    | Desc.                              |
| -------------- | ------ | ------- | ---------------------------------- |
| $V_{max}$      | 4.46   | mM * Hz | Maximal rate                       |
| $\Delta\Psi_0$ | 91     | mV      | Offset potential                   |
| $K_{act}$      | 3.8E-4 | mM      | Activation constant for calcium    |
| $K_{trans}$    | 0.019  | mM      | Dissociation constant for calcium  |
| n              | -2.8   | -       | Activation cooperativity           |
| L              | 110    | -       | Keq for conformational transitions |

## Mitochondrial sodium-calcium exchanger (NCLX)[^Wei2011]

$$
\begin{aligned}
J_{NCLX} &= V_{max}\text{exp}(b\Delta\Psi_m/V_T)\frac{[Ca^{2+}]_m}{[Ca^{2+}]_i} (\frac{A}{1+A})^n \frac{B}{1+B}  \\
A &= [Na^+]_i / K_{Na}  \\
B &= [Ca^{2+}]_m / K_{Ca}  \\
\end{aligned}
$$

| Parameter | Value   | Unit    | Desc.                             |
| --------- | ------- | ------- | --------------------------------- |
| $V_{max}$ | 0.183   | mM * Hz | Maximal rate                      |
| b         | 0.5     | -       | Ffraction of $\Delta\Psi_m$       |
| $K_{Na}$  | 9.4     | mM      | Dissociation constant for sodium  |
| $K_{Ca}$  | 3.75E-4 | mM      | Dissociation constant for calcium |
| $n$       | 3       |         |                                   |

Mitochondrial proton leak
$$
J_{hleak} = g_H\Delta\Psi_m
$$

| Parameter | Value | Unit         | Desc.                                   |
| --------- | ----- | ------------ | --------------------------------------- |
| $g_{H}$   | 2     | mM / (V * s) | Ionic conductance of the inner membrane |


## Mitochondrial hydrogen flux balance[^Wei2011]

- $J_H$: Proton influx to mitochondrial matrix by pumps / transporters
- $J_{Hn}$: Proton flux due to enzyme stoichiometry
- $J_{HL}$: Proton flux due to ligand binding / unbinding

$$
\begin{aligned}
J_H &= -J_{h, Res}+ J_{hu} + J_{NHE} + J_{PiC} + J_{Hleak}   \\
J_{Hn} &= -(J_{IDH3} + J_{KGDH} + J_{MDH} - J_{F1Fo})   \\
J_{HL} &= \frac{[H^+]_m}{K_{a, ATP}P_{ATP}}\frac{d[ATP]_m}{dt} + \frac{[H^+]_m}{K_{a, ADP}P_{ADP}}\frac{d[ADP]_m}{dt} + \frac{[H^+]_m}{K_{a, Pi}P_{Pi}}\frac{d[Pi]_m}{dt} + \frac{[H^+]_m}{K_{a, SUC}P_{SUC}}\frac{d[SUC]_m}{dt}   \\
\frac{d[H^+]_m}{dt} &= δ_H(J_H - J_{Hn} - J_{HL})   \\
\end{aligned}
$$

## ODEs for mitochondrial ions

$$
\begin{aligned}
\frac{d [Ca^{2+}]_m}{dt} &=\delta_{Ca}( J_{uni} - J_{NCLX}) \\
\frac{d [Na^+]_m}{dt} &= J_{NCLX} - J_{NaH} \\
\frac{d [Pi]_m}{dt} &= -J_{F1Fo} + J_{PiC} - J_{SL}  \\
C_{m}\frac{d \Delta \Psi_m}{dt} &= J_{Hres} - J_{Hu} - J_{ANT} - J_{Hleak} -J_{NCLX} - J_{uni} - J_{IMAC} \\
\end{aligned}
$$

[^Wei2011]: Wei AC, Aon MA, O'Rourke B, Winslow RL, Cortassa S. Mitochondrial energetics, pH regulation, and ion dynamics: a computational-experimental approach. Biophys J. 2011;100(12):2894-903. [PMC3123977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3123977/)
