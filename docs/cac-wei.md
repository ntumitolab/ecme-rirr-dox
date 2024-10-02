# TCA cycle rates[^Wei2011]

## Citrate synthase (CS)

$$
\begin{aligned}
J_{CS} &= \frac{k_{cat} E_T AB}{(1+A)(1+B)}  \\
A &= [AcCoA] / K_m^{AcCoA} \\
B &= [OAA] / K_m^{OAA}
\end{aligned}
$$

| Parameter        | Value   | Unit | Description                         |
| :--------------- | ------- | ---- | ----------------------------------- |
| $k_{cat}$        | 0.23523 | Hz   | Catalytic constant                  |
| $E_T$            | 0.4     | mM   | Enzyme concentration of CS          |
| $K_m^{AcCoA}$    | 0.0126  | mM   | Michaelis constant for AcCoA        |
| $K_m^{OAA}$      | 6.4E-4  | mM   | Michaelis constant for OAA          |
| $[AcCoA]$        | 1       | mM   | Acetyl CoA concentration            |
| $k_{cat}$ (cell) | 0.15891 | Hz   | Catalytic constant (cellular model) |

## Aconitase (ACO)

$$
\begin{aligned}
J_{ACO} &= k_f ([CIT] - [ISOC] / K_{eq})  \\
\ [CIT] &= \Sigma_{CAC} - [ISOC] - [\alpha KG]-[SCoA] - [SUC] - [FUM] - [MAL] - [OAA]
\end{aligned}
$$

| Parameter      | Value    | Unit | Description                            |
| -------------- | -------- | ---- | -------------------------------------- |
| $k_f$          | 0.11688  | Hz   | Forward rate constant of ACO           |
| $K_{eq}$       | 2.22     | -    | Equilibrium constant of ACO            |
| $\Sigma_{CAC}$ | 1.300    | mM   | Sum of TCA cycle intermediates         |
| $k_f$ (cell)   | 0.078959 | Hz   | Forward rate constant (cellular model) |

## Isocitrate dehydrogenase, NADH-producing (IDH3)

$$
\begin{aligned}
J_{IDH3} &= \frac{k_{cat} E_T AB}{f_H AB + f_i B + f_a A + f_a f_i} \\
f_H & = 1 + \frac{[H^+]_m}{K_{H1}} + \frac{K_{H2}}{[H^+]_m}  \\
A &= [NAD] / K_{NAD} \\
B &= ([ISOC] / K_{ISOC})^n  \\
f_a &= \frac{K_A}{K_A + [ADP]_m} \frac{K_{CA}}{K_{CA} + [Ca^{2+}]_m}  \\
f_i &= 1 + \frac{[NADH]}{K_{NADH}}  \\
\end{aligned}
$$

| Parameter        | Value | Unit | Description                       |
| ---------------- | ----- | ---- | --------------------------------- |
| $k_{cat}$        | 11880 | Hz   | Rate constant of IDH3             |
| $E_T$            | 0.109 | mM   | Concentration of IDH3             |
| $K_{H1}$         | 1E-6  | mM   | Ionization constant of IDH3       |
| $K_{H2}$         | 9E-4  | mM   | Ionization constant of IDH3       |
| $K_{NAD}$        | 0.923 | mM   | Michaelis constant for NAD        |
| $K_{ISOC}$       | 1.520 | mM   | Michaelis constant for isocitrate |
| $n$              | 2     | -    | Cooperativity for isocitrate      |
| $K_A$            | 0.62  | mM   | Activation constant by ADP        |
| $K_{CA}$         | 5E-4  | mM   | Activation constant for calcium   |
| $K_{NADH}$       | 0.19  | mM   | Inhibition constant by NADH       |
| $k_{cat}$ (cell) | 535   | Hz   | Rate constant (cellular model)    |

## Alpha-ketoglutarate dehydrogenase (KGDH)

$$
\begin{aligned}
J_{KGDH} &= \frac{k_{cat} E_T AB}{f_H AB + f_a (A + B)} \\
f_H & = 1 + \frac{[H^+]_m}{K_{H1}} + \frac{K_{H2}}{[H^+]_m}  \\
A &= [NAD] / K_{NAD} \\
B &= ([\alpha KG] / K_{AKG})^n  \\
f_a &= \frac{K_{MG}}{K_{MG} + [Mg^{2+}]_m} \frac{K_{CA}}{K_{CA} + [Ca^{2+}]_m}  \\
\end{aligned}
$$

| Parameter        | Value  | Unit | Description                    |
| ---------------- | ------ | ---- | ------------------------------ |
| $k_{cat}$        | 13.2   | Hz   | Rate constant of KGDH          |
| $E_T$            | 0.5    | mM   | Concentration of KGDH          |
| $K_{H1}$         | 4E-5   | mM   | Ionization constant of KGDH    |
| $K_{H2}$         | 7E-5   | mM   | Ionization constant of KGDH    |
| $K_{NAD}$        | 38.7   | mM   | Michaelis constant for NAD     |
| $K_{AKG}$        | 30     | mM   | Michaelis constant for αKG     |
| $n$              | 1.2    | -    | Hill coefficient for αKG       |
| $K_{MG}$         | 0.0308 | mM   | Activation constant for Mg     |
| $K_{CA}$         | 1.5E-4 | mM   | Activation constant for Ca     |
| $k_{cat}$ (cell) | 17.9   | Hz   | Rate constant (cellular model) |

## Succinate-CoA ligase (SL)

$$
\begin{aligned}
J_{SL} &= k_f ([SCoA][ADP]_m[Pi]_m - [SUC][ATP]_m[CoA]/K_{eq}^{app}) \\
K_{eq}^{app} &= K_{eq} \frac{P_{SUC}P_{ATP}}{P_{Pi}P_{ADP}}
\end{aligned}
$$

| Parameter    | Value   | Unit    | Description                            |
| ------------ | ------- | ------- | -------------------------------------- |
| $k_f$        | 2.8E-5  | mM * Hz | Forward rate constant of SL            |
| $K_{eq}$     | 3.115   | -       | Equilibrium constant of SL             |
| [CoA]        | 0.020   | mM      | Coenzyme A concentration               |
| $k_f$ (cell) | 2.84E-5 | mM * Hz | Forward rate constant (cellular model) |

## Succinate dehydrogenase (SDH)

See OXPHOS part: complex II (Succinate dehydrogenase).

## Fumarate hydratase (FH)

$$
J_{FH} = k_f ([FUM] - [MAL] / K_{eq})
$$

| Parameter    | Value | Unit | Description                            |
| ------------ | ----- | ---- | -------------------------------------- |
| $k_f$        | 8.3   | Hz   | Forward rate constant                  |
| $K_{eq}$     | 1.0   | -    | Equilibrium constant                   |
| $k_f$ (cell) | 8.4   | Hz   | Forward rate constant (cellular model) |

## Malate dehydrogenase (MDH)

$$
\begin{aligned}
J_{MDH} &= \frac{k_{cat} E_T AB f_a f_i}{(1+A)(1+B)}  \\
A &= \frac{[MAL]}{K_{MAL}}\frac{K_{OAA}}{K_{OAA} + [OAA]}  \\
B &= [NAD] / K_{NAD}  \\
f_a &= k_{offset} + \left( 1 + \frac{[H^+]_m}{K_{H1}} (1 + \frac{[H^+]_m}{K_{H2}})    \right)^{-1}  \\
f_i &= \left( 1 + \frac{K_{H3}}{[H^+]_m} (1 + \frac{K_{H4}}{[H^+]_m})    \right)^{2}  \\
\end{aligned}
$$

| Parameter        | Value    | Units | Description                          |
| ---------------- | -------- | ----- | ------------------------------------ |
| $k_{cat}$        | 124.2    | Hz    | Rate constant                        |
| $E_T$            | 0.154    | mM    |                                      |
| $K_{H1}$         | 1.131E-5 | mM    | Ionization constant                  |
| $K_{H2}$         | 26.7     | mM    | Ionization constant                  |
| $K_{H3}$         | 6.68E-9  | mM    | Ionization constant                  |
| $K_{H4}$         | 5.62E-6  | mM    | Ionization constant                  |
| $k_{offset}$     | 0.0399   |       | Offset of MDH pH activation factor   |
| $K_{NAD}$        | 0.2244   | mM    | Michaelis constant for NAD           |
| $K_{MAL}$        | 1.493    | mM    | Michaelis constant for malate        |
| $K_{OAA}$        | 0.031    | mM    | Inhibition constant for oxaloacetate |
| $k_{cat}$ (cell) | 125.9    | Hz    | Rate constant for cellular model     |

## Aspartate aminotransferase (AAT)

$$
J_{AAT} = k_f [OAA][GLU] \frac{k_{ASP}K_{eq}}{k_{ASP}K_{eq} + k_f[\alpha KG]}
$$

| Parameter    | Value  | Units | Description                            |
| ------------ | ------ | ----- | -------------------------------------- |
| $k_f$        | 21.4   | Hz    | Forward rate constant                  |
| $k_{ASP}$    | 0.0015 | Hz    | Rate constant of aspartate consumption |
| $K_{eq}$     | 6.6    |       | Equilibrium constant                   |
| [GLU]        | 30.000 | mM    | Glutamate concentration                |
| $k_f$ (cell) | 21.7   | Hz    | Forward rate constant (cellular model) |

## ODEs specifically to the Citric acid cycle (CAC)

$$
\begin{aligned}
\frac{d [ISOC]}{dt} &= J_{ACO} -J_{IDH3} -J_{IDH2}  \\
\frac{d [\alpha KG]}{dt} &= J_{IDH3} + J_{IDH2} - J_{KGDH} + J_{AAT}  \\
\frac{d[SCoA]}{dt} &= J_{KGDH} - J_{SL}  \\
\frac{d [SUC]}{dt} &= J_{SL} - J_{SDH} \\
\frac{d [FUM]}{dt} &= J_{SDH} - J_{FH}  \\
\frac{d[MAL]}{dt} &= J_{FH} - J_{MDH}  \\
\frac{d [OAA]}{dt} & = J_{MDH} - J_{CS} - J_{AAT}  \\
\frac{d [NADH]_m}{dt} &= -J_{C1} + J_{IDH} + J_{KGDH} + J_{MDH} - J_{THD}  \\
\end{aligned}
$$

[^Wei2011]: Wei AC, Aon MA, O'Rourke B, Winslow RL, Cortassa S. Mitochondrial energetics, pH regulation, and ion dynamics: a computational-experimental approach. Biophys J. 2011;100(12):2894-903. [PMC3123977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3123977/)