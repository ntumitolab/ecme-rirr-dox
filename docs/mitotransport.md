## Adenine Nucleotide translocator (ANT)

$$
\begin{aligned}
J_{ANT} &= V_{max}^{ANT}\frac{AB - \delta PQ}{(B + \delta^{h_{ANT}} P)(A + Q)}  \\
A &= [ATP^{4-}]_m = 0.025 [ATP]_m  \\
B &= [ADP^{3-}]_i = 0.45 [ADP]_i \\
P &= [ATP^{4-}]_i = 0.25 [ATP]_i \\
Q &= [ADP^{3-}]_m = 0.17 [ADP]_m \\
\delta &= \text{exp}(-\Delta\Psi_m F / RT)
\end{aligned}
$$


| Parameter | Value | Unit    | Desc.            |
| --------- | ----- | ------- | ---------------- |
| $V_{max}^{ANT}$ | 5  | mM * Hz | Maximal rate of ANT |
| $h_{ANT}$       | 0.5   | -       | Fraction of MMP |

## Mitochondrial calcium uniporter (MCU)

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

## Mitochondrial sodium-calcium exchanger (NCLX)

$$
J_{NCLX} = V_{max}^{NCLX}\exp(b\Delta\Psi_m F/RT)\frac{[Ca^{2+}]_m}{[Ca^{2+}]_i} (\frac{[Na^+]_i}{[Na^+]_i + K_{Na}^{NCLX}})^n \frac{[Ca^{2+}]_m}{[Ca^{2+}]_m + K_{Ca}^{NCLX}}
$$

| Parameter | Value   | Unit    | Desc.                             |
| --------- | ------- | ------- | --------------------------------- |
| $V_{max}^{NCLX}$ | 0.04665 | mM/s    | Maximal rate of NCLX              |
| b         | 0.5     | -       | Fraction of MMP                   |
| $K_{Na}^{NCLX}$  | 9.4     | mM      | Dissociation constant for sodium  |
| $K_{Ca}^{NCLX}$  | 3.75E-4 | mM      | Dissociation constant for calcium |
| $n$       | 3       |         |                                   |

## Mitochondrial proton leak

$$
J_{hleak} = g_H\Delta\Psi_m
$$

| Parameter | Value | Unit         | Desc.                                   |
| --------- | ----- | ------------ | --------------------------------------- |
| $g_{H}$   | 2     | mM / (Volt * s) | Ionic conductance of the inner mitochondrial membrane |

## ODEs for mitochondrial ion transport

$$
\begin{aligned}
\frac{d [Ca^{2+}]_m}{dt} &=\delta_{Ca}( J_{uni} - J_{NCLX}) \\
\frac{d [Na^+]_m}{dt} &= J_{NCLX} - J_{NaH} \\
\frac{d [Pi]_m}{dt} &= -J_{F1Fo} + J_{PiC} - J_{SL}  \\
C_{m}\frac{d \Delta \Psi_m}{dt} &= J_{Hres} - J_{Hu} - J_{ANT} - J_{Hleak} -J_{NCLX} - J_{uni} - J_{IMAC} \\
\end{aligned}
$$
