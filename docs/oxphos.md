# OXPHOS

## Complex I model[^Gauthier2013A]

Assuming single electron transfer for each cycle.

$$
\begin{aligned}
\nu &= \text{exp}((\Delta\Psi_m - \Delta\Psi_B) / V_T) \\
a_{12} &= k_{12} ([H^+]_m)^2 \\
a_{21} &= k_{21} \\
a_{65} &= k_{65} ([H^+]_i)^2  \\
a_{56} &= k_{56} \\
a_{61} &= k_{61} / \nu \\
a_{16} &= k_{16} \nu  \\
a_{23} &= k_{23} \sqrt{[NADH]} \\
a_{32} &= k_{32}  \\
a_{34} &= k_{34}  \\
a_{43} &= k_{43} \sqrt{[NAD^+]}  \\
a_{47} &= C1_{inhib}  \cdot K_{47} \sqrt{[Q_n][H^+]_m} \\
a_{74} &= k_{74} \\
a_{57} &=  C1_{inhib}  \cdot K_{57} \sqrt{[QH_2]} \\
a_{75} &= k_{75} \\
k_{42}^\prime &= k_{42} \\
a_{42} &= k_{42}^\prime [O_2] \\
K_{eq}^{ROS} &= \text{exp}((E_{FMN} - E_{sox}) / V_T) \\
a_{24} &= a_{42} K_{eq}^{ROS} [O_2^{ \cdot  -}]_m  \\
a_{25} &= a_{52} = 0  \\
\end{aligned}
$$

$$
\begin{aligned}
e_{1} &=
a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{74} +
a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{75} +
a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{61} \cdot a_{74}  \\ &+
a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74} +
a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{74} +
a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{75}  \\ &+
a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{61} \cdot a_{74} +
a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74} +
a_{21} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75}  \\ &+
a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{74} +
a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{75} +
a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{61} \cdot a_{74}  \\ &+
a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74} +
a_{21} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75} +
a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75}  \\ &+
a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75} +
a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75}  \\
e_{2} &=
a_{12} \cdot a_{32} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{74}  +
a_{12} \cdot a_{32} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{75}  +
a_{12} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{61} \cdot a_{74}  \\ &+
a_{12} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74}  +
a_{12} \cdot a_{32} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{74}  +
a_{12} \cdot a_{32} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{75}  \\ &+
a_{12} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{61} \cdot a_{74}  +
a_{12} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}  +
a_{12} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75}  \\ &+
a_{12} \cdot a_{34} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{74}  +
a_{12} \cdot a_{34} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{75}  +
a_{12} \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{61} \cdot a_{74}  \\ &+
a_{12} \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74}  +
a_{12} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75}  +
a_{16} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74}  \\ &+
a_{16} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}  +
a_{16} \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74}  \\
e_{3} &=
a_{12} \cdot a_{23} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{23} \cdot a_{42} \cdot a_{56} \cdot a_{61} \cdot a_{75}   +
a_{12} \cdot a_{23} \cdot a_{42} \cdot a_{57} \cdot a_{61} \cdot a_{74}   \\ &+
a_{12} \cdot a_{23} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{12} \cdot a_{23} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{23} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{75}   \\ &+
a_{12} \cdot a_{23} \cdot a_{43} \cdot a_{57} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{23} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{12} \cdot a_{23} \cdot a_{47} \cdot a_{56} \cdot a_{61} \cdot a_{75}   \\ &+
a_{12} \cdot a_{24} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{24} \cdot a_{43} \cdot a_{56} \cdot a_{61} \cdot a_{75}   +
a_{12} \cdot a_{24} \cdot a_{43} \cdot a_{57} \cdot a_{61} \cdot a_{74}   \\ &+
a_{12} \cdot a_{24} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{23} \cdot a_{42} \cdot a_{57} \cdot a_{65} \cdot a_{74}   \\ &+
a_{16} \cdot a_{23} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{24} \cdot a_{43} \cdot a_{57} \cdot a_{65} \cdot a_{74}   \\
e_{4} &=
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{56} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{56} \cdot a_{61} \cdot a_{75}   +
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{57} \cdot a_{61} \cdot a_{74}   \\ &+
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{56} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{56} \cdot a_{61} \cdot a_{75}   \\ &+
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{57} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{56} \cdot a_{61} \cdot a_{74}   \\ &+
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{56} \cdot a_{61} \cdot a_{75}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{57} \cdot a_{61} \cdot a_{74}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{57} \cdot a_{65} \cdot a_{74}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{23} \cdot a_{34} \cdot a_{57} \cdot a_{65} \cdot a_{74}   \\ &+
a_{16} \cdot a_{24} \cdot a_{32} \cdot a_{57} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{24} \cdot a_{34} \cdot a_{57} \cdot a_{65} \cdot a_{74}   \\
e_{5} &=  
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{61} \cdot a_{75}   +
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{65} \cdot a_{75}   +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{61} \cdot a_{75}   \\ &+
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{65} \cdot a_{75}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{61} \cdot a_{75}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{65} \cdot a_{75}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{65} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{65} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{65} \cdot a_{74}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{65} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{47} \cdot a_{65} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{65} \cdot a_{74}   \\ &+
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{65} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{47} \cdot a_{65} \cdot a_{75}   +
a_{16} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{65} \cdot a_{75}   \\ &+
a_{16} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{65} \cdot a_{75}   +
a_{16} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{65} \cdot a_{75}   \\
e_{6} &=  
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{75} +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{75}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{75}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{56} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{56} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{74}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{56} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{56} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{74}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{75}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{56} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{56} \cdot a_{75}   \\ &+
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{74}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{75}   +
a_{16} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{75}   \\ &+
a_{16} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{75}   +
a_{16} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{75}   \\
e_{7} &=  
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{61} +
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{61}   +
a_{12} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{65}   \\ &+
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{56} \cdot a_{61}   +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{57} \cdot a_{61}   +
a_{12} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{57} \cdot a_{65}   \\ &+
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{56} \cdot a_{61}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{61}   +
a_{12} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{65}   \\ &+
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{42} \cdot a_{57} \cdot a_{65}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{43} \cdot a_{57} \cdot a_{65}   +
a_{16} \cdot a_{21} \cdot a_{32} \cdot a_{47} \cdot a_{57} \cdot a_{65}   \\ &+
a_{16} \cdot a_{21}choco  \cdot a_{34} \cdot a_{42} \cdot a_{57} \cdot a_{65}   +
a_{16} \cdot a_{21} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{65}   +
a_{16} \cdot a_{23} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{65}   \\ &+
a_{16} \cdot a_{24} \cdot a_{32} \cdot a_{47} \cdot a_{57} \cdot a_{65}   +
a_{16} \cdot a_{24} \cdot a_{34} \cdot a_{47} \cdot a_{57} \cdot a_{65}   \\
\end{aligned}
$$

$$
\begin{aligned}
Δ &= e_{1} + e_{2} + e_{3} + e_{4} + e_{5} + e_{6} + e_{7}  \\
\rho_{C1}^\prime &= \rho_{C1}  \cdot mt_{prot} / \Delta  \\
J_{Hres}^{C1} &= 2\rho_{C1}^\prime (e_{6}a_{61} - e_{1}a_{16}) \\
J_{Q}^{C1} &= 0.5\rho_{C1}^\prime (e_{4}a_{47} - e_{7}a_{74}) \\
J_{NADH}^{C1} &= 0.5\rho_{C1}^\prime (e_{3}a_{34} - e_{4}a_{43}) \\
J_{ROS}^{C1} &= \rho_{C1}^\prime (e_{4}a_{42} - e_{2}a_{24})  \\
\end{aligned}
$$

| Parameter      | Value     | Units       | Desc.                                        |
| -------------- | --------- | ----------- | -------------------------------------------- |
| $\rho_{C1}$    | 8.849     | mM          | Concentration of complex I<br />(Adjustable) |
| $\Delta\Psi_B$ | 50        | mV          | Phase boundary potential                     |
| $k_{12}$       | 6.3396E11 | Hz/mM^2     |                                              |
| $k_{21}$       | 5         | Hz          |                                              |
| $k_{56}$       | 100       | Hz          |                                              |
| $k_{65}$       | 2.5119E13 | Hz/mM^2     |                                              |
| $k_{61}$       | 1E7       | Hz          |                                              |
| $k_{16}$       | 130       | Hz          |                                              |
| $k_{23}$       | 3886.7    | Hz/mM^{1/2} |                                              |
| $k_{32}$       | 9.1295E6  | Hz          |                                              |
| $k_{34}$       | 639.1364  | Hz          |                                              |
| $k_{43}$       | 3.2882    | Hz/mM^{1/2} |                                              |
| $k_{47}$       | 1.5962E7  | Hz/mM       |                                              |
| $k_{74}$       | 65.2227   | Hz          |                                              |
| $k_{75}$       | 24615     | Hz          |                                              |
| $k_{57}$       | 1166.7    | Hz/mM^{1/2} |                                              |
| $k_{42}$       | 6.0318    | Hz/mM       |                                              |
| $E_{FMN}$      | -375      | mV          | Midpoint potential of flavin mononucleotide  |
| $E_{sox}$      | -150      | mV          | Midpoint potential of superoxide             |

## Complex II (Succinate dehydrogenase)[^Gauthier2013A]

$$
\begin{aligned}
f_Q &= \frac{[Q]_n}{[Q]_n + [QH_2]_n}  \\
f_{OAA} &= \frac{K_i}{[OAA] + K_i} \\
f_{SUC} &= 0.085 \sqrt{[SUC] / [FUM]} \\
J_{SDH} &= V_{SDH} C2_{inhib} f_{SUC} f_{OAA} \frac{f_Q}{f_Q + K_m}  \\
J_{c2} &= J_{SDH} \\
\end{aligned}
$$

| Parameter | Value | Units   | Desc.                                |
| --------- | ----- | ------- | ------------------------------------ |
| $V_{SDH}$ | 4.167 | mM * Hz | Maximum rate of SDH                  |
| $K_i$     | 0.150 | mM      | Inhibition constant for oxaloacetate |
| $K_m$     | 0.6   | -       | Michaelis constant for CoQ           |

## Complex III[^Gauthier2013A]

$$
\begin{aligned}
f_{hi} & = [H^+]_{i}  / 10^{-7}M   \\
v_{1} &= v_{Q}^{C1} + v_{Q}^{C2}   \\
v_2 &= k_d([QH_2]_{n} - [QH_2]_{p})  \\
k_{3} &= k_{03}K_{eq3}f_{hi} \\
k_{-3} &= k_{03} \\
v_{3} &= k_3[QH_2]_{p} [FeS]_{ox} - k_{-3}[Q^{ \cdot  -}]_{p} [FeS]_{rd}  \\
k_{4, ox} &= k_{04}K_{eq4, ox} \text{exp}(-\alpha\delta_1\Delta\Psi_m / V_T) \\
k_{4, rd} &= k_{04}K_{eq4, rd} \text{exp}(-\alpha\delta_1\Delta\Psi_m / V_T) \\
k_{-4, ox} &= k_{04} \text{exp}(\alpha(1-\delta_1)\Delta\Psi_m / V_T)  \\
k_{-4, rd} &= k_{04} \text{exp}(\alpha(1-\delta_1)\Delta\Psi_m / V_T)  \\
v_{4, ox} &= k_{4, ox}[Q^{ \cdot -}]_{p} [b1] - k_{-4, ox}[Q]_{p} [b2]  \\
v_{4, rd} &= k_{4, rd}[Q^{ \cdot -}]_{p} [b3] - k_{-4, rd}[Q]_{p} [b4]  \\
v_{5} &= k_d([Q]_{p} - [Q]_{n})  \\
k_{6} &= K_{06}K_{eq6}\text{exp}(-\beta\delta_2\Delta\Psi_m / V_T) \\
k_{-6} &= k_{06} \text{exp}(\beta(1-\delta_2)\Delta\Psi_m / V_T)  \\
v_{6} &= k_{6} [b2] - k_{-6} [b3]  \\
k_{7, ox} &= k_{07, ox}K_{eq7, ox}\text{exp}(-\gamma\delta_3\Delta\Psi_m / V_T) \\
k_{7, rd} &= k_{07, rd}K_{eq7, rd}\text{exp}(-\gamma\delta_3\Delta\Psi_m / V_T) \\
k_{-7, ox} &= k_{07, ox} \text{exp}(\gamma(1-\delta_3)\Delta\Psi_m / V_T)  \\
k_{-7, rd} &= k_{07, rd} \text{exp}(\gamma(1-\delta_3)\Delta\Psi_m / V_T)  \\
v_{7, ox} &= (k_{7, ox}[Q]_{n}[b3] - k_{-7, ox}[Q^{ \cdot  -}]_{n}[b1])C3_{inhib} \\
v_{7, rd} &= (k_{7, rd}[Q]_{n}[b4] - k_{-7, rd}[Q^{ \cdot  -}]_{n}[b2])C3_{inhib} \\
f_{hm} & = [H^+]_{m}   / 10^{-7}M  \\
k_{8, ox} &= k_{08, ox}K_{eq8, ox}\text{exp}(-\gamma\delta_3\Delta\Psi_m / V_T)(f_{hm})^2  \\
k_{8, rd} &= k_{08,rd}K_{eq8, rd}\text{exp}(-\gamma\delta_3\Delta\Psi_m / V_T)(f_{hm})^2  \\
k_{-8, ox} &= k_{08, ox} \text{exp}(\gamma(1-\delta_3)\Delta\Psi_m / V_T)  \\
k_{-8, rd} &= k_{08, rd} \text{exp}(\gamma(1-\delta_3)\Delta\Psi_m / V_T)  \\
v_{8, ox} &= (k_{8, ox}[Q^{ \cdot  -}]_{n}[b3] - k_{-8, ox}[QH_2]_{n}[b1])C3_{inhib} \\
v_{8, rd} &= (k_{8, rd}[Q^{ \cdot  -}]_{n}[b4] - k_{-8, rd}[QH_2]_{n}[b2])C3_{inhib} \\
k_9 &= k_{09}K_{eq9}  \\
k_{-9} &= k_{09} \\
v_{9} &= k_{9}[FeS]_{rd}[cytc1]_{ox} - k_{-9}[FeS]_{ox}[cytc1]_{rd}\\
k_{10} &= k_{010}K_{eq10}  \\
k_{-10} &= k_{010} \\
v_{10} &= k_{10}[Q^{ \cdot  -}]_p[O_2] - k_{-10}[Q]_p[O_2^{ \cdot -}]  \\
v_{10b} &= v_{10}  \\
v_{33} &= k_{33}(K_{eq}[cytc1]_{rd}[cytc]_{ox} - [cytc]_{rd}[cytc1]_{ox})  \\
\rho_{C3}^\prime &= \rho_{C3} \cdot mt_{prot}  \\
\rho_{C4}^\prime &= \rho_{C4} \cdot mt_{prot}  \\
FeS_{rd} &= \rho_{C3}^\prime - FeS_{ox} \\
cytc1_{rd} &= \rho_{C3}^\prime - cytc1_{ox}  \\
cytc_{rd} &= \rho_{C4}^\prime - cytc_{ox}  \\
 [b4] &= \rho_{C3}^\prime - [b1] - [b2] - [b3]  \\
 [QH_2]_p &= \Sigma[Q] - [Q]_n - [Q]_p - [QH_2]_n - [Q^{ \cdot  -}]_p - [Q^{ \cdot  -}]_n  \\
J_{hRes}^{C3} &= 2v_{3}     \\
J_{ROS, m}^{C3} &= v_{10}   \\
J_{ROS, i}^{C3} &= v_{10b}  \\
\end{aligned}
$$

| Parameter    | Value    | Unit  | Desc.                                                     |
| ------------ | -------- | ----- | --------------------------------------------------------- |
| $k_{03}$     | 1,666.63 | Hz/mM | Reverse rate constant for reaction 3                      |
| $K_{eq3}$    | 0.6877   | -     | Equilibrium constant for reaction 3                       |
| $k_{04}$     | 60.67    | Hz/mM | Reverse rate constant for reaction 4                      |
| $K_{eq4,ox}$ | 129.9853 | -     | Equilibrium constant for reaction 4 <br />(bH oxidized)   |
| $K_{eq4,rd}$ | 13.7484  | -     | Equilibrium constant for reaction 4 <br />(bH reduced)    |
| $\delta_1$   | 0.5      | -     |                                                           |
| $\alpha$     | 0.2497   | -     |                                                           |
| $k_d$        | 22000    | Hz    | Rate of diffusion across the membrane <br />for Q and QH2 |
| $k_{06}$     | 166.67   | Hz/mM | Reverse rate constant for reaction 6                      |
| $K_{eq6}$    | 9.4596   | -     | Equilibrium constant for reaction 6                       |
| $\delta_2$   | 0.5      | -     |                                                           |
| $\beta$      | 0.5006   | -     |                                                           |
| $k_{07,ox}$  | 13.33    | Hz/mM | Reverse rate constant for reaction 7 <br />(bL oxidized)  |
| $K_{eq7,ox}$ | 3.0748   | -     | Equilibrium constant for reaction 7 <br />(bL oxidized)   |
| $k_{07,rd}$  | 1.667    | Hz/mM | Reverse rate constant for reaction 7 <br />(bL reduced)   |
| $K_{eq7,rd}$ | 29.0714  | -     | Equilibrium constant for reaction 7 <br />(bL reduced)    |
| $\delta_3$   | 0.5      | -     |                                                           |
| $\gamma$     | 0.2497   | -     | $\alpha + \beta + \gamma = 1$                             |
| $k_{08,ox}$  | 83.33    | Hz/mM | Reverse rate constant for reaction 8 <br />(bL oxidized)  |
| $K_{eq8,ox}$ | 129.9853 | -     | Equilibrium constant for reaction 8 <br />(bL oxidized)   |
| $k_{08,rd}$  | 8.333    | Hz/mM | Reverse rate constant for reaction 8 <br />(bL reduced)   |
| $K_{eq8,rd}$ | 9.4596   | -     | Equilibrium constant for reaction 8 <br />(bL reduced)    |
| $k_{09}$     | 833      | Hz/mM | Reverse rate constant for reaction 9                      |
| $K_{eq9}$    | 0.2697   | -     | Equilibrium constant for reaction 9                       |
| $k_{010}$    | 0.8333   | Hz/mM | Reverse rate constant for reaction 10                     |
| $K_{eq10}$   | 1.4541   | -     | Equilibrium constant for reaction 10                      |
| $k_{33}$     | 2469.13  | Hz/mM | Reverse rate constant for reaction 33                     |
| $K_{eq33}$   | 2.1145   | -     | Equilibrium constant for reaction 33                      |
| $\rho_{C3}$  | 0.325    | mM    | Total complex III protein                                 |


## Complex IV[^Gauthier2013A]

$$
\begin{aligned}
f_{H_{m}} &= \text{exp}(-\delta_5\Delta\Psi_m / V_T) ([H^+]_m /10^{-7}M) \\
f_{H_{i}} &= \text{exp}((1-\delta_5)\Delta\Psi_m / V_T) ([H^+]_i /10^{-7}M)  \\
f_{C_{rd}} &= [cytc]_{rd} \\
f_{C_{ox}} &= \text{exp}((1-\delta_5)\Delta\Psi_m / V_T) [cytc]_{ox} \\
a_{12} &= k_{34} f_{C_{rd}}^3 f_{H_{m}}^4  \\
a_{14} &= k_{-37} f_{H_{i}}   \\
a_{21} &= k_{-34} f_{C_{ox}}^3 f_{H_{i}}  \\
a_{23} &= k_{35} [O_2] C4_{inhib}  \\
a_{34} &= k_{36} f_{C_{rd}} f_{H_{m}}^3  \\
a_{41} &= k_{37} f_{H_{m}}  \\
a_{43} &= k_{-36} f_{C_{ox}} f_{H_{i}}^2 \\
e_1 &= a_{21}a_{41}a_{34} + a_{41}a_{34}a_{23}  \\
e_2 &= a_{12}a_{41}a_{34}  \\
e_3 &= a_{23}a_{12}a_{41} + a_{43}a_{14}a_{21} + a_{23}a_{43}a_{12} + a_{23}a_{43}a_{14}  \\
e_4 &= a_{14}a_{34}a_{21} + a_{34}a_{23}a_{12} + a_{34}a_{23}a_{14}  \\
\Delta &= e_1 + e_2 + e_3+ e_4  \\
Y &= e_1 / \Delta  \\
Yr &= e_2 / \Delta  \\
YO &= e_3 / \Delta  \\
YOH &= e_4 / \Delta  \\
v_{34} &= \rho_{C4}^\prime (Y  \cdot a_{12} - Yr  \cdot a_{21})  \\
v_{35} &= \rho_{C4}^\prime Yr  \cdot a_{23}  \\
v_{36} &= \rho_{C4}^\prime (YO  \cdot a_{34} - YOH  \cdot a_{43})  \\
v_{37} &= \rho_{C4}^\prime (YOH  \cdot a_{41} - Y  \cdot a_{14})  \\
V_e &= 3v_{34} + v_{35}  \\
J_{hRes}^{C4}  &= v_{34} + 2v_{36} + v_{37}  \\
J_{O_2} &= v_{35}  \\
J_{hRes} &= J_{hRes}^{C1} + J_{hRes}^{C3} + J_{hRes}^{C4} \\
\rho_{C4}^\prime &= \rho_{C4} \cdot mt_{prot}
\end{aligned}
$$

| Parameter     | Value     | Unit    | Desc.                    |
| ------------- | --------- | ------- | ------------------------ |
| $\Sigma cytc$ | 0.325     | mM      | Cytochrome c pool        |
| $\rho_{C4}$   | 0.325     | mM      | Complex IV concentration |
| $k_{34}$      | 2.9445E10 | Hz/mM^3 | @ pH = 7                 |
| $k_{-34}$     | 290.03    | Hz/mM^3 | @ pH = 7                 |
| $k_{35}$      | 45000     | Hz/mM   |                          |
| $k_{36}$      | 4.826E11  | Hz/mM   | @ pH = 7                 |
| $k_{-36}$     | 4.826     | Hz/mM   | @ pH = 7                 |
| $k_{37}$      | 1.7245E8  | Hz      | @ pH = 7                 |
| $k_{-37}$     | 17.542    | Hz      | @ pH = 7                 |


## Complex V (ATP synthase) [^Wei2011]

$$
\begin{aligned}
J_{F1Fo} &= -\rho^{F1} ((100 p_a + p_{c1} v_B) v_a - (p_a + p_{c2} v_a) v_h)  / \Delta \\
J_H^{F1Fo} &= -3\rho^{F1} (100p_a(1 + v_a) - (p_a + p_b)v_h) / \Delta \\
\Delta &= (1 + p_1 v_a)v_B + (p_2 + p_3 v_a)v_h \\
v_B &= \text{exp}(3\Delta\Psi_B / V_T)   \\
v_h &= \text{exp}(3\Delta p / V_T)  \\
v_a &= \frac{K_{eq}^{'} \cdot \Sigma[ATP]_m}{ \Sigma[Pi]_m \cdot \Sigma[ADP]_m } \\
\end{aligned}
$$

| Parameter      | Value     | Unit | Desc.                                                                                   |
| -------------- | --------- | ---- | --------------------------------------------------------------------------------------- |
| $\rho_{F1}$    | 5         | mM   | Concentration of F1-Fo ATPase                                                           |
| $K_{eq}^{'}$   | 6.47E5    | M    | Apparent equilibrium constant for ATP hydrolysis<br />From Golding's work[^Golding1995] |
| $\Delta\Psi_B$ | 50        | mV   | Phase boundary potential                                                                |
| $p_{a}$        | 1.656E-5  | Hz   | Sum of products of rate constants                                                       |
| $p_{b}$        | 3.373E-7  | Hz   | Sum of products of rate constants                                                       |
| $p_{c1}$       | 9.651E-14 | Hz   | Sum of products of rate constants                                                       |
| $p_{c2}$       | 4.585E-14 | Hz   | Sum of products of rate constants                                                       |
| $p_{1}$        | 1.346E-4  | -    | Sum of products of rate constants                                                       |
| $p_{2}$        | 7.739E-7  | -    | Sum of products of rate constants                                                       |
| $p_{3}$        | 6.65E-15  | -    | Sum of products of rate constants                                                       |

## ODE system for the Q cycle

$$
\begin{aligned}
\frac{d[Q]_n}{dt} &= v_5 - v_{7,ox}- v_{7,rd} - v_1  \\
\frac{d[Q^{ \cdot -}]_n}{dt} &= v_{7,ox} + v_{7,rd} - v_{8,ox}- v_{8,rd}  \\
\frac{d[QH_2]_n}{dt} &= v_{8,ox} + v_{8,rd} + v_1 - v_2   \\
\frac{d[QH_2]_p}{dt} &= v_2 -v_3 \\
\frac{d[Q^{ \cdot -}]_p}{dt} &= v_3 - v_{10} - v_{10b} - v_{4,ox} - v_{4,rd}   \\
\frac{d[Q]_p}{dt} &= v_{10} + v_{10b} + v_{4,ox} + v_{4,rd} - v_5   \\
\frac{d[b1]}{dt} &= v_{7,ox} + v_{8,ox} - v_{4,ox}    \\
\frac{d[b2]}{dt} &= v_{4,ox} + v_{7,rd} - v_{8,rd} - v_6   \\
\frac{d[b3]}{dt} &= v_6 - v_{4,rd} + v_{7,ox} - v_{8,ox}    \\
\frac{d[b4]}{dt} &= v_{4,rd} - v_{7,rd} - v_{8,rd}   \\
\frac{d[FeS]_{ox}}{dt} &= v_9 - v_3      \\
\frac{d[cytc1]_{ox}}{dt} &= v_{33} - v_9   \\
\frac{d[cytc]_{ox}}{dt} &= V_e - v_{33}   \\
\end{aligned}
$$


[^Wei2011]: Wei AC, Aon MA, O'Rourke B, Winslow RL, Cortassa S. Mitochondrial energetics, pH regulation, and ion dynamics: a computational-experimental approach. Biophys J. 2011;100(12):2894-903. [PMC3123977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3123977/)

[^Golding1995]: Golding, E. M., Teague, W. E., & Dobson, G. P. (1995). Adjustment of K’ to varying pH and pMg for the creatine kinase, adenylate kinase and ATP hydrolysis equilibria permitting quantitative bioenergetic assessment. The Journal of Experimental Biology, 198(Pt 8), 1775–1782.

[^Gauthier2013A]: Gauthier LD, Greenstein JL, O’Rourke B, Winslow RL. An Integrated Mitochondrial ROS Production and Scavenging Model: Implications for Heart Failure. Biophysical Journal. 2013;105(12):2832-2842. doi:10.1016/j.bpj.2013.11.007. [PMC3882515](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882515)
