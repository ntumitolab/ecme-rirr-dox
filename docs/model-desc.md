# RIRR-DOX model description

## TCA cycle rates

### Conservation relationship

$$
\begin{align}
\Sigma_{CAC} = [CIT] + [ISOC] + [\alpha KG] + [SCoA] + [SUC] + [FUM] + [MAL] + [OAA]
\end{align}
$$

| Parameter      | Value | Unit | Description                    |
| -------------- | ----- | ---- | ------------------------------ |
| $\Sigma_{CAC}$ | 1.300 | mM   | Sum of TCA cycle intermediates |

### Citrate synthase (CS)

$$
\begin{align}
J_{CS} = \frac{k_{cat}^{CS} E_T^{CS} ([AcCoA] / K_m^{AcCoA})([OAA] / K_m^{OAA})}{(1+[AcCoA] / K_m^{AcCoA})(1+[OAA] / K_m^{OAA})}
\end{align}
$$

| Parameter      | Value   | Unit | Description                  |
| :------------- | ------- | ---- | ---------------------------- |
| $k_{cat}^{CS}$ | 0.23523 | Hz   | Catalytic constant           |
| $E_T^{CS}$     | 0.4     | mM   | Enzyme concentration of CS   |
| $K_m^{AcCoA}$  | 12.6    | μM   | Michaelis constant for AcCoA |
| $K_m^{OAA}$    | 0.64    | μM   | Michaelis constant for OAA   |
| $[AcCoA]$      | 0.1     | mM   | Acetyl CoA concentration     |

### Aconitase (ACO)

$$
J_{ACO} = k_f^{ACO} ([CIT] - [ISOC] / K_{eq}^{ACO})
$$

| Parameter      | Value | Unit | Description                  |
| -------------- | ----- | ---- | ---------------------------- |
| $k_f^{ACO}$    | 0.1   | Hz   | Forward rate constant of ACO |
| $K_{eq}^{ACO}$ | 2.22  | -    | Equilibrium constant of ACO  |

### Isocitrate dehydrogenase, NADH-producing (IDH3)

$$
\begin{align}
J_{IDH3} &= \frac{k_{cat}^{IDH3} E_T^{IDH3} AB}{f_H AB + f_i B + f_a A + f_a f_i} \\
f_H & = 1 + \frac{[H^+]_m}{K_{H1}^{IDH3}} + \frac{K_{H2}^{IDH3}}{[H^+]_m}  \\
A &= [NAD] / K_{NAD}^{IDH3} \\
B &= ([ISOC] / K_{ISOC}^{IDH3})^n_{IDH3}  \\
f_a &= \frac{K_A^{IDH3}}{K_A^{IDH3} + [ADP]_m} \frac{K_{CA}^{IDH3}}{K_{CA}^{IDH3} + [Ca^{2+}]_m}  \\
f_i &= 1 + \frac{[NADH]}{K_{NADH}^{IDH3}}  \\
\end{align}
$$

| Parameter         | Value | Unit | Description                       |
| ----------------- | ----- | ---- | --------------------------------- |
| $k_{cat}^{IDH3}$  | 535   | Hz   | Rate constant of IDH3             |
| $E_T^{IDH3}$      | 0.109 | mM   | Concentration of IDH3             |
| $K_{H1}^{IDH3}$   | 1     | nM   | Ionization constant of IDH3       |
| $K_{H2}^{IDH3}$   | 900   | nM   | Ionization constant of IDH3       |
| $K_{NAD}^{IDH3}$  | 0.923 | mM   | Michaelis constant for NAD        |
| $K_{ISOC}^{IDH3}$ | 1.520 | mM   | Michaelis constant for isocitrate |
| $n_{IDH3}$        | 2     | -    | Cooperativity for isocitrate      |
| $K_A^{IDH3}$      | 0.62  | mM   | Activation constant by ADP        |
| $K_{CA}^{IDH3}$   | 0.5   | μM   | Activation constant for calcium   |
| $K_{NADH}^{IDH3}$ | 0.19  | mM   | Inhibition constant by NADH       |

### Alpha-ketoglutarate dehydrogenase (KGDH)

$$
\begin{align}
J_{KGDH} &= \frac{k_{cat}^{KGDH} E_T^{KGDH} AB}{f_H AB + f_a (A + B)} \\
f_H & = 1 + \frac{[H^+]_m}{K_{H1}^{KGDH}} + \frac{K_{H2}^{KGDH}}{[H^+]_m}  \\
A &= [NAD] / K_{NAD}^{KGDH} \\
B &= ([\alpha KG] / K_{AKG}^{KGDH})^n_{KGDH}  \\
f_a &= \frac{K_{MG}^{KGDH}}{K_{MG}^{KGDH} + [Mg^{2+}]_m} \frac{K_{CA}^{KGDH}}{K_{CA}^{KGDH} + [Ca^{2+}]_m}  \\
\end{align}
$$

| Parameter        | Value | Unit | Description                 |
| ---------------- | ----- | ---- | --------------------------- |
| $k_{cat}^{KGDH}$ | 17.9  | Hz   | Rate constant of KGDH       |
| $E_T^{KGDH}$     | 0.5   | mM   | Concentration of KGDH       |
| $K_{H1}^{KGDH}$  | 40    | nM   | Ionization constant of KGDH |
| $K_{H2}^{KGDH}$  | 70    | nM   | Ionization constant of KGDH |
| $K_{NAD}^{KGDH}$ | 38.7  | mM   | Michaelis constant for NAD  |
| $K_{AKG}^{KGDH}$ | 30    | mM   | Michaelis constant for αKG  |
| $n_{KGDH}$       | 1.2   | -    | Hill coefficient for αKG    |
| $K_{MG}^{KGDH}$  | 30.8  | μM   | Activation constant for Mg  |
| $K_{CA}^{KGDH}$  | 0.15  | μM   | Activation constant for Ca  |

### Succinate-CoA ligase (SL)

$$
\begin{align}
J_{SL} &= k_f^{SL} ([SCoA][ADP]_m[Pi]_m - [SUC][ATP]_m[CoA]/K_{eq}^{app}) \\
K_{eq}^{app} &= K_{eq}^{SL} \frac{P_{SUC}P_{ATP}}{P_{Pi}P_{ADP}} \\
\end{align}
$$

| Parameter     | Value | Unit    | Description                 |
| ------------- | ----- | ------- | --------------------------- |
| $k_f^{SL}$    | 28.4  | 1Hz/mM² | Forward rate constant of SL |
| $K_{eq}^{SL}$ | 3.11  | -       | Equilibrium constant of SL  |
| [CoA]         | 0.020 | mM      | Coenzyme A concentration    |

### Succinate dehydrogenase (SDH)

See OXPHOS part: complex II (Succinate dehydrogenase).

### Fumarate hydratase (FH)

$$
J_{FH} = k_f^{FH} ([FUM] - [MAL] / K_{eq}^{FH})
$$

| Parameter     | Value | Unit | Description           |
| ------------- | ----- | ---- | --------------------- |
| $k_f^{FH}$    | 8.3   | Hz   | Forward rate constant |
| $K_{eq}^{FH}$ | 1.0   | -    | Equilibrium constant  |

### Malate dehydrogenase (MDH)

$$
\begin{align}
J_{MDH} &= \frac{k_{cat}^{MDH} E_T^{MDH} AB f_a f_i}{(1+A)(1+B)}  \\
A &= \frac{[MAL]}{K_{MAL}^{MDH}}\frac{K_{OAA}^{MDH}}{K_{OAA}^{MDH} + [OAA]}  \\
B &= [NAD] / K_{NAD}^{MDH}  \\
f_a &= k_{offset}^{MDH} + \left( 1 + \frac{[H^+]_m}{K_{H1}^{MDH}} (1 + \frac{[H^+]_m}{K_{H2}^{MDH}})    \right)^{-1}  \\
f_i &= \left( 1 + \frac{K_{H3}^{MDH}}{[H^+]_m} (1 + \frac{K_{H4}^{MDH}}{[H^+]_m})    \right)^{2}  \\
\end{align}
$$

| Parameter          | Value  | Units | Description                          |
| ------------------ | ------ | ----- | ------------------------------------ |
| $k_{cat}^{MDH}$    | 125.9  | Hz    | Rate constant                        |
| $E_T^{MDH}$        | 154    | μM    |                                      |
| $K_{H1}^{MDH}$     | 11.31  | nM    | Ionization constant                  |
| $K_{H2}^{MDH}$     | 26.7   | mM    | Ionization constant                  |
| $K_{H3}^{MDH}$     | 6.68   | pM    | Ionization constant                  |
| $K_{H4}^{MDH}$     | 5.62   | nM    | Ionization constant                  |
| $k_{offset}^{MDH}$ | 0.0399 | -     | Offset of MDH pH activation factor   |
| $K_{NAD}^{MDH}$    | 224.4  | μM    | Michaelis constant for NAD           |
| $K_{MAL}^{MDH}$    | 1.493  | mM    | Michaelis constant for malate        |
| $K_{OAA}^{MDH}$    | 31     | μM    | Inhibition constant for oxaloacetate |

### Aspartate aminotransferase (AAT)

$$
\begin{align}
J_{AAT} = k_f^{AAT} [OAA][GLU] \frac{k_{ASP}^{AAT} K_{eq}^{AAT}}{k_{ASP}^{AAT} K_{eq}^{AAT} + k_f[\alpha KG]}
\end{align}
$$

| Parameter       | Value  | Units | Description                            |
| --------------- | ------ | ----- | -------------------------------------- |
| $k_f^{AAT}$     | 21.7   | Hz/mM | Forward rate constant                  |
| $k_{ASP}^{AAT}$ | 0.0015 | Hz    | Rate constant of aspartate consumption |
| $K_{eq}^{AAT}$  | 6.6    | -     | Equilibrium constant                   |
| [GLU]           | 30     | mM    | Glutamate concentration                |

### ODEs in the citric acid cycle

$$
\begin{align}
\frac{d [ISOC]}{dt} &= J_{ACO} -J_{IDH3} -J_{IDH2}  \\
\frac{d [\alpha KG]}{dt} &= J_{IDH3} + J_{IDH2} - J_{KGDH} + J_{AAT}  \\
\frac{d [SCoA]}{dt} &= J_{KGDH} - J_{SL}  \\
\frac{d [SUC]}{dt} &= J_{SL} - J_{SDH} \\
\frac{d [FUM]}{dt} &= J_{SDH} - J_{FH}  \\
\frac{d [MAL]}{dt} &= J_{FH} - J_{MDH}  \\
\frac{d [OAA]}{dt} & = J_{MDH} - J_{CS} - J_{AAT}  \\
\end{align}
$$

## Endoplasmic reticulum

### Ryanodine receptor (Jrel)

$$
\begin{align}
P_{C1} &= 1 - P_{O1} - P_{O2} - P_{C2}  \\
v_{o1c1} &= -k_a^-P_{O1} + k_a^+[Ca^{2+}]_{ss}^n P_{C1} \\
v_{o1o2} &= k_b^+ [Ca^{2+}]_{ss}^m_{RyR} P_{O1} - k_b^- P_{O2} \\
v_{o1c2} &= k_c^+ P_{O1} - k_c^- P_{C2} \\
\dot{P_{O1}}  &= -v_{o1c1} - v_{o1o2} - v_{o1c2}  \\
\dot{P_{O2}}  &= v_{o1o2}  \\
\dot{P_{C2}}  &= v_{o1c2}  \\
J_{rel} &= r_{ryr} (P_{O1} + P_{O2})([Ca^{2+}]_{JSR} - [Ca^{2+}]_{ss}) \\
\end{align}
$$

| Parameter | Value   | Units  | Description               |
| --------- | ------- | ------ | ------------------------- |
| $r_{RyR}$ | 3600    | Hz     | RyR flux channel constant |
| $n_{RyR}$ | 4       | -      | Cooperativity parameter   |
| $m_{RyR}$ | 3       | -      | Cooperativity parameter   |
| $k_a^+$   | 12.15   | Hz/μM⁴ | RyR rate constant         |
| $k_a^-$   | 576     | Hz     | RyR rate constant         |
| $k_b^+$   | 0.00405 | Hz/μM³ | RyR rate constant         |
| $k_b^-$   | 1930    | Hz     | RyR rate constant         |
| $k_c^+$   | 100     | Hz     | RyR rate constant         |
| $k_c^-$   | 0.8     | Hz     | RyR rate constant         |

### SERCA (Jup)

Michaelis-Menten dependence of enzyme activity with respect to ATP and mixed-type inhibition of the enzyme by ADP. Reversible during diastole with low cytoplasmic calcium levels.

$$
\begin{align}
J_{up} &= \frac{V_{f}^{up}f_b-V_{r}^{up}r_b}{(1 + f_b + r_b)f_{ATP}^{SERCA}} \\
f_b &= \left( \frac{[Ca^{2+}]_i}{K_{fb}} \right)^{N_{fb}} \\
r_b &= \left( \frac{[Ca^{2+}]_{NSR}}{K_{rb}} \right)^{N_{rb}} \\
f_{ATP}^{SERCA} &= \frac{K_{m,up}^{ATP}}{[ATP]_i} ( \frac{[ADP]_i}{K_{i1, up}} + 1) + \frac{[ADP]_i}{K_{i2, up}} + 1 \\
\end{align}
$$

| Parameter            | Value   | Units | Description                                    |
| -------------------- | ------- | ----- | ---------------------------------------------- |
| $V_{max, f}^{SERCA}$ | 0.2989  | Hz*mM | SERCA forward rate parameter                   |
| $V_{max, b}^{SERCA}$ | 0.3179  | Hz*mM | SERCA reverse  rate parameter                  |
| $K_{f}^{SERCA}$      | 0.24    | μM    | Forward Ca2+ half-saturation constant of SERCA |
| $K_{r}^{SERCA}$      | 1.64269 | mM    | Reverse Ca2+ half-saturation constant of SERCA |
| $N_{f}^{SERCA}$      | 1.4     | -     | Forward cooperativity constant of SERCA        |
| $N_{r}^{SERCA}$      | 1.0     | -     | Reverse  cooperativity constant of SERCA       |
| $K_{ATP}^{SERCA}$    | 10      | μM    | ATP half-saturation constant for SERCA         |
| $K_{ADP1}^{SERCA}$   | 140     | μM    | ADP first inhibition constant for SERCA        |
| $K_{ADP2}^{SERCA}$   | 5.1     | mM    | ADP second  inhibition constant for SERCA      |

## Sarcoplasmic ion currents

GHK current equation

$$
\begin{align}
\Phi_s(P_s, z_s, V_m, [S]_i, [S]_o) := P_sz^2_s\frac{V_mF^2}{RT}\frac{[S]_i - [S]_o\exp(-z_sV_mF/RT)}{1-\exp(-z_sV_mF/RT)}
\end{align}
$$

### Time-dependent delayed rectifier potassium current (IK)

$$
\begin{align}
I_K &= \bar G_K X_1 X_K^2 (V - E_K) \\
E_K &= \frac{RT}{F} \ln \frac{[K^+]_o + P_{Na,K}[Na^+]_o}{ [K^+]_i + P_{Na,K}[Na^+]_i} \\
\bar G_K &= 0.282 (mS/cm^2)\sqrt{[K^+]_o /5.4mM} \\
X_1 &= (1+ e^{(V_m-40)/40})^{-1} \\
\frac{dX_k}{dt} &= \alpha_X - X_k (\alpha_X + \beta_X) \\
\alpha_X &= \frac{V_m+30}{1 - e^{-0.148(V_m+30)}} * 0.0719Hz \\
\beta_X &= \frac{V_m+30}{e^{0.0687(V_m+30)} -1} * 0.131Hz \\
\end{align}
$$

### Time-independent potassium current (IK1)

$$
\begin{align}
\Delta V &= V_m - E_{K1} \\
I_{K1} &= \bar G_{K1}K_{1 \infty}\Delta V \\
E_{K1} &= \frac{RT}{F} \ln \frac{[K^+]_o}{[K^+]_i}  \\
\bar G_{K1} &= 0.748(mS/cm^2)\sqrt{[K^+]_o / 5.4mM}  \\
K_{1 \infty} &= \frac{\alpha_{K_1}}{\alpha_{K_1} + \beta_{K_1}} \\
\alpha_{K_1} &= \frac{1.02}{1 + e^{0.2385(\Delta V -59.215)}} * \text{kHz}  \\
\beta_{K_1} &= \frac{0.4912e^{0.28032(\Delta V + 5.476)} + e^{0.06175(\Delta V -594.31)}}{1 + e^{-0.5143(\Delta V + 4.753)}} * \text{kHz}
\end{align}
$$

### Plateau potassium current (IKp)

$$
\begin{align}
E_{Kp} &= \frac{RT}{F} \ln \frac{[K^+]_o}{[K^+]_i} \\
I_{Kp} &= \frac{\bar G_{Kp} (V - E_{Kp})}{1 + e^{(7.488-V_m) / 5.98}} \\
\end{align}
$$

### Fast Na current (INa)

$$
\begin{align}
I_{Na} &= \bar G_{Na} m_{Na}^{3} h_{Na} j_{Na} (V_m-E_{Na}) \\
E_{Na} &= \frac{RT}{F} \ln \frac{[Na^{+}]_o}{[Na^{+}]_i} \\
\frac{dm_{Na}}{dt} &= \alpha_{m} - m_{Na}(\alpha_{m} + \beta_{m}) \\
\frac{dh_{Na}}{dt} &= \alpha_{h} - h_{Na}(\alpha_{h} + \beta_{h}) \\
\frac{dj_{Na}}{dt} &= \alpha_{j} - m_{Na}(\alpha_{j} + \beta_{j}) \\
\alpha_{m} &= 0.32kHz \frac{V_m + 47.13}{1 - e^{-0.1(V_m+47.13)}} \\
\beta_{m} &= 0.08kHz \times e^{-V_m / 11} \\
\\
For \ V_m & \ge -40mV \\
\alpha_{h} &= \alpha_{j} = 0 \\
\beta_{h} &= (0.13 ms (1+e^{-(V_m+10.66)/11.1}))^{-1} \\
\beta_{j} &= 0.3kHz\frac{e^{-2.535 \times 10^{-7}V_m}}{1 + e^{-0.1(V_m + 32)}} \\
\\
For \ V_m & < -40mV \\
\alpha_{h} &= 0.135kHz * e^{-(V_m+80)/6.8} \\
\alpha_{j} &= (-127140e^{0.2444 V_m}-3.474 \times 10^{-5}e^{-0.04391 V_m})\frac{V_m + 37.78}{1+e^{0.311( V_m +79.23)}} \times \text{kHz} \\
\beta_{h} &= (3.56e^{0.079 V_m} + 3.1 \times 10^{5}e^{0.35 V_m}) \times \text{kHz} \\
\beta_{j} &= \frac{0.1212e^{-0.01052 V_m}}{1+e^{-0.1378(V_m + 40.14)}} \times \text{kHz} \\
\end{align}
$$

### Sodium-calcium exchanger current (INaCa)

$$
\begin{align}
I_{NaCa} &= k_{NaCa}  \cdot f_{Nao }  \cdot f_{Cao }\frac{exp(V_mF/RT)\phi_{Na}^3 - \phi_{Ca}}{exp((1 - \eta) V_mF/RT ) + k_{sat}} \\
f_{Nao} &= \frac{([Na^+]_o)^3}{([Na^+]_o)^3 + (K_{M,Na}^{NaCa})^3} \\
f_{Cao} &= \frac{[Ca^+]_o}{[Ca^+]_o + K_{M,Ca}^{NaCa}} \\
\phi_{Na} &= \frac{[Na^+]_i}{ [Na^+]_o} \\
\phi_{Ca} &= \frac{[Ca^{2+}]_i}{[Ca^{2+}]_o} \\
\end{align}
$$

### Background calcium ($I_{Ca,b}$) and sodium currents ($I_{Na,b}$)

$$
\begin{align}
I_{Ca,b} &= \bar G_{Ca,b} (V_m - \frac{RT}{2F} \ln \frac{[Ca^{2+}]_o}{[Ca^{2+}]_i})  \\
I_{Na,b} &= \bar G_{Na,b} (V_m - \frac{RT}{F} \ln \frac{[Na^{+}]_o}{[Na^{+}]_i}) \\
\end{align}
$$

### Non-specific calcium-activated current (InsCa)

$$
\begin{align}
f_{Ca} &= \frac{([Ca^{2+}]_i)^3}{([Ca^{2+}]_i)^3 + (K_{m}^{nsCa})^3}\\
I_{nsNa} &= 0.75  \cdot f_{Ca}  \cdot  \Phi_{Na}(P_{nsNa}, z_{Na}, V_m, [Na^+]_i, [Na^+]_o)  \\
I_{nsK} &= 0.75  \cdot f_{Ca}  \cdot  \Phi_{K}(P_{nsK}, z_{K}, V_m, [K^+]_i, [K^+]_o)  \\
\end{align}
$$

### Sodium-potassium ATPase current (INaK)

The Na+/K+ ATPase activity depends on the ATP concentration, as well as the competitive inhibition by ADP.

$$
\begin{align}
I_{NaK} &= \bar I_{NaK}  \cdot f_{ATP}  \cdot f_{Na}  \cdot f_{K} \cdot f_{NaK}  \\
\sigma &= \frac{e^{[Na^+]_o / 67.3mM}-1}{7}  \\
f_{NaK} &= (1 + 0.1245 \cdot \exp(-0.1V_m F / RT) + 0.0365 \sigma \cdot \exp(-V_m F / RT))^{-1}  \\
f_{Na} &= \frac{([Na^+]_i)^{1.5}}{([Na^+]_i)^{1.5} + (K_{m, Na_i})^{1.5}} \\
f_{K} &= \frac{[K^+]_o}{[K^+]_o + K_{m, K_o}} \\
f_{ATP} &= \frac{[ATP]_i}{[ATP]_i + K_{M,ATP}^{NaK} / f_{ADP}} \\
f_{ADP} &= \frac{K_{i,ADP}^{NaK}}{K_{i,ADP}^{NaK} + [ADP]_i} \\
\end{align}
$$

### L-type Ca current (ICa & ICaK)

"Common pool" subspace calcium model.

$$
\begin{align}
\alpha &= 0.4 e^{(V_m+2) / 10}  \\
\beta &= 0.4 e^{-(V_m+2) / 13}  \\
\alpha^\prime  &=  a \alpha \\
\beta^\prime  &=  \beta / b \\
\gamma &= \gamma_0 [Ca^{2+}]_{ss}  \\
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
\end{align}
$$

$$
\begin{align}
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
P_{K}  &= P_{K}^{max} \frac{I_{Ca}^{half}}{I_{Ca}^{half} + I_{Ca}^{max}}  \\
y_\infty &= \frac{1}{1 + e^{(V_m + 55) / 7.5}} + \frac{0.5}{1 + e^{(-V_m + 21) / 6}}  \\
\tau_y &= 20ms + \frac{600ms}{1 + e^{(V_m + 30) / 9.5}}  \\
\frac{dy_{Ca}}{dt}  &=  \frac{y_\infty - y_{Ca}}{\tau_y}  \\
\end{align}
$$

| Parameter      | Value                 | Units            | Description                                |
| -------------- | --------------------- | ---------------- | ------------------------------------------ |
| $A$            | 2                   |                  | Mode transition parameter                  |
| $B$            | 2                   |                  | Mode transition parameter                  |
| $\gamma_0$     | 187.5              | Hz/μM            | Mode transition parameter                  |
| $\omega$       | 10                  | Hz               | Mode transition parameter                  |
| $f$            | 300                 | Hz               | Transition rate into open state            |
| $g$            | 2000                | Hz               | Transition rate into open state            |
| $P_{Ca}^{LCC}$ | $8 \cdot 10^{-4}$  | cm/s             | L-type Ca2+ channel permeability to Ca2+ (*)  |
| $P_{K}^{LCC}$  | $1.11 \cdot 10^{-11}$ | cm/s             | L-type Ca2+ channel permeability to K+     |
| $I_{Ca, half}$ | $-0.4583$             | $\mu A / cm^{2}$ | ICa level that reduces equation Pk by half |

(*): adjusted for proper calcium transients.

### Plasma membrane calcium ATPase (PMCA) current (IpCa)

Modified rate expression incorporating the ATP-dependence of pump activity.
Plasma membrane calcium ATPase (PMCA) rate exhibits two different K0.5 values for ATP.

$$
\begin{align}
I_{pCa} &= I_{max}^{PMCA} \times \frac{[Ca^{2+}]_i}{[Ca^{2+}]_i +  K_{M, Ca}^{PMCA}} \times f_{ATP} \\
f_{ATP} &= \frac{[ATP]_i}{[ATP]_i + K_{M2,ATP}^{PMCA}} + \frac{[ATP]_i}{[ATP]_i + K_{M1,ATP}^{PMCA} / f{ADP}} \\
f_{ADP} &= \frac{K_{i,ADP}^{PMCA}}{K_{i,ADP}^{PMCA} + [ADP]_i} \\
\end{align}
$$

| Parameter         | Value   | Units        | Description                                                   |
| ----------------- | ------- | ------------ | ------------------------------------------------------------- |
| $I_{max}^{PMCA}$  | $0.575$ | $\mu A/cm^2$ | Maximum sarcolemmal Ca2+ pump current                         |
| $K_{Ca}^{PMCA}$   | $0.5$   | $uM$         | Ca2+ half-saturation constant for sarcolemmal Ca2+ pump       |
| $K_{ATP1}^{PMCA}$ | $0.012$ | $mM$         | First ATP half-saturation constant for sarcolemmal Ca2+ pump  |
| $K_{ATP2}^{PMCA}$ | $0.23$  | $mM$         | Second ATP half-saturation constant for sarcolemmal Ca2+ pump |
| $K_{ADP}^{PMCA}$  | $1.0$   | $mM$         | ADP inhibition constant for sarcolemmal Ca2+ pump             |

### Electrophysiology ODEs

$$
\begin{align}
\frac{d[Na^+]_i}{dt} &= -(I_{Na} + 3I_{NaCa} + 3I_{NaK})\frac{A_{cap}}{V_{myo}F} + (V_{NHE} - 3V_{NaCa}) \frac{V_{mito}}{V_{myo}} \\
\frac{d[K^+]_i}{dt} &= -(I_{Ks} + I_{Kr} + I_{K1} + I_{Kp} + I_{Ca,K}-2I_{NaK})\frac{A_{cap}}{V_{myo}F} \\
C_m\frac{dV_m}{dt} &= -(I_{Na} + I_{CaL} + I_{Kr} + I_{Ks} + I_{K1} + I_{Kp} + I_{NaCa} + I_{NaK} + I_{pCa} + I_{Ca, b} + I_{K_{ATP}} + I_{stim}) \\
β_i &= \frac{(K_m^{CMDN} + [Ca^{2+}]_i)^2}{ (K_m^{CMDN} + [Ca^{2+}]_i)^2 + K_m^{CMDN}  \cdot  [CMDN]_{tot}} \\
β_{SR} &= \frac{(K_m^{CSQN} + [Ca^{2+}]_{SR})^2}{(K_m^{CSQN} + [Ca^{2+}]_{SR})^2 + K_m^{CSQN}  \cdot  [CSQN]_{tot}} \\
\frac{d[Ca^{2+}]_i}{dt} &= \beta_i(J_{xfer}\frac{V_{ss}}{V_{myo}} - J_{up} - J_{trpn} - (I_{Ca,b} -2I_{NaCa} + I_{pCa})\frac{A_{cap}}{2V_{myo}F} + (V_{NaCa} - V_{uni})\frac{V_{mito}}{V_{myo}}) \\
\frac{d[Ca^{2+}]_{SR}}{dt} &= \beta_{SR}(J_{up}\frac{V_{myo}}{V_{SR}} - J_{rel}\frac{V_{ss}}{V_{SR}}) \\
\end{align}
$$

| Symbol          | Value                | Units        | Description                                           |
| --------------- | -------------------- | ------------ | ----------------------------------------------------- |
| $G_{Na}$        | $12.8$               | $mS/cm^2$    | Maximal Na channel conductance                        |
| $G_{Kp}$        | $0.00828$            | $mS/cm^2$    | Maximal plateau K channel conductance                 |
| $G_{K,0}$       | $0.282$              | $mS/cm^2$    | IK conductance                                        |
| $G_{K1,0}$      | $0.748$              | $mS/cm^2$    | IK1 conductance                                       |
| $P_{NaK}$       | $0.01833$            |              | Na+ permeability ratio of K+ channel                  |
| $K_{NaCa}$      | $9000$               | $\mu A/cm^2$ | NCX current                                           |
| $K_{Na}^{NCX}$  | $87.5$               | $mM$         | Dissociation constant of sodium for NCX               |
| $K_{Ca}^{NCX}$  | $1.38$               | $mM$         | Dissociation constant of calcium for NCX              |
| $K_{sat}^{NCX}$ | $0.1$                |              | NCX saturation factor at negative potentials          |
| $\eta^{NCX}$    | $0.35$               |              | Voltage dependence of NCX                             |
| $P_{ns,Na}$     | $1.75 \cdot 10^{-7}$ | $cm/s$       | Nonspecific channel current Na permeability           |
| $P_{ns,K}$      | $0$                  | $cm/s$       | Nonspecific channel current K permeability            |
| $K_{ca}^{ns}$   | $1.2$                | $\mu M$      | Ca2+ half-saturation constant for nonspecific current |
| $G_{Ca,b}$      | $0.003217$           | $mS/cm^2$    | Maximum background current Ca2+ conductance           |
| $G_{Na,b}$      | $0.003217$           | $mS/cm^2$    | Maximum background current Na+ conductance            |
| $\tau_{tr}$     | $574.7$              | Hz           | Time constant for transfer from subspace to myoplasm  |
| $\tau_{xfer}$   | $9090$               | Hz           | Time constant for transfer from NSR to JSR            |
| $K_{m}^{CMDN}$  | $2.38$               | $\mu M$      | Ca2+ half saturation constant for calmodulin          |
| $K_{m}^{CSQN}$  | $800$                | $\mu M$      | Ca2+ half saturation constant for calsequestrin       |
| $\Sigma[HTRPN]$ | $140$                | $\mu M$      | Total troponin high-affinity sites                    |
| $\Sigma[LTRPN]$ | $70$                 | $\mu M$      | Total troponin low-affinity sites                     |
| $\Sigma[CMDN]$  | $50$                 | $\mu M$      | Total myoplasmic calmodulin concentration             |
| $\Sigma[CQSN]$  | $15$                 | $mM$         | Total NSR calsequestrin concentration                 |

## Force generation

The rate of ATP hydrolysis associated with force generation through actomyosin ATPase depends explicitly on both ATP and ADP. [^Rice2000]

$$
\begin{align}
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
\end{align}
$$

$$
\begin{align}
v_{01} &= f_{01} [P_0] - g_{01(SL)} [P_1] \\
v_{12} &= f_{12} [P_1] - g_{21(SL)} [P_2]  \\
v_{23} &= f_{23} [P_2] - g_{23(SL)} [P_3] \\
v_{04} &= k_{pn}^{trop} [P_0] - k_{np}^{trop} [N_0] \\
 [N_0] &= 1 - [P_0] - [P_1] - [P_2] - [P_3] - [N_1]  \\
v_{15} &= k_{pn}^{trop} [P_1] - k_{np}^{trop} [N_1] \\
v_{54} &= g_{01,off} [N_1]  \\
 [HTRPN] &=  [HTRPN]_{tot} - [HTRPNCa]  \\
 [LTRPN] &=  [LTRPN]_{tot} - [LTRPNCa]  \\
f_{ATP}^{AM} &= \frac{[ATP]_i}{[ATP]_i + K_{m,AM}^{ATP}/f_{ADP}^{AM} } \\
f_{ADP}^{AM} &= \frac{K_{i,AM}^{ADP}}{[ADP]_i + K_{i,AM}^{ADP}} \\
V_{AM} &= V_{max}^{AM}  \cdot f_{ATP}^{AM}  \cdot  \frac{f_{01}[P_0] + f_{12}[P_1] + f_{23}[P_2]}{f_{01} + f_{12} + f_{23}} \\
J_{trpn} &= \frac{d[HTRPNCa]}{dt} + \frac{d[LTRPNCa]}{dt} \\
\frac{d[HTRPNCa]}{dt} &= k^{+}_{htrpn}[Ca^{2+}]_i[HTRPN] - k^{-}_{htrpn}[HTRPNCa]  \\
\frac{d[LTRPNCa]}{dt} &= k^{+}_{ltrpn}[Ca^{2+}]_i[LTRPN] - k^{-}_{ltrpn}(1-\frac{2}{3}Force_{norm})[LTRPNCa]  \\
\frac{d[P_0]}{dt} &= - v_{01} - v_{04}  \\
\frac{d[P_1]}{dt} &= v_{01} - v_{12} - v_{15}  \\
\frac{d[P_2]}{dt} &= v_{12} - v_{23}  \\
\frac{d[P_3]}{dt} &= v_{23}  \\
\frac{d[N_1]}{dt} &= v_{15} - v_{54} \\
\end{align}
$$


[^Rice2000]: Rice JJ, Jafri MS, Winslow RL. Modeling short-term interval-force relations in cardiac muscle. Am J Physiol Heart Circ Physiol. 2000 Mar;278(3):H913-31. [APS](https://www.physiology.org/doi/full/10.1152/ajpheart.2000.278.3.H913)

| Symbol          | Value    | Units           | Description                                                   |
| --------------- | -------- | --------------- | ------------------------------------------------------------- |
| $k_{pn}^{trop}$ | $40$     | $\text{Hz}$     | Transition rate from tropomyosin permissive to non-permissive |
| $\text{SL}$     | $2.15$   | $\mu \text{m}$  | Sarcomere length                                              |
| $f_{XB}$        | $50$     | $\text{Hz}$     | Transition rate from weak to strong crossbridge               |
| $g_{XB}^{min}$  | $100$    | $\text{Hz}$     | Minimum transition rate from strong to weak crossbridge       |
| $\zeta$         | $0.1$    | $\text{N/mm}^2$ | Conversion factor normalizing to physiological force          |
| $V_{AM}^{max}$  | $7.2$    | $\text{mM/s}$   | Conversion factor normalizing to physiological force          |
| $K_{ATP}^{AM}$  | $0.03$   | $\text{mM}$     | ATP half-saturation constant of AM ATPase                     |
| $K_{ADP}^{AM}$  | $0.26$   | $\text{mM}$     | ADP inhibition constant of AM ATPase                          |
| $h_{trpn}^{+}$  | $100000$ | $\text{Hz/mM}$  | Ca2+ on-rate for troponin high-affinity sites                 |
| $h_{trpn}^{-}$  | $0.33$   | $\text{Hz}$     | Ca2+ off-rate for troponin high-affinity sites                |
| $l_{trpn}^{+}$  | $100000$ | $\text{Hz/mM}$  | Ca2+ on-rate for troponin low-affinity sites                  |
| $l_{trpn}^{-}$  | $40$     | $\text{Hz}$     | Ca2+ off-rate for troponin low-affinity sites                 |

## OXPHOS

### Complex I

Assuming single electron transfer for each redox reaction.

$$
\begin{align}
\nu &= \exp((\Delta\Psi_m - \Delta\Psi_B) F/ RT) \\
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
K_{eq}^{ROS} &= \exp((E_{FMN} - E_{sox}) F / RT) \\
a_{24} &= a_{42} K_{eq}^{ROS} [O_2^{ \cdot  -}]_m  \\
a_{25} &= a_{52} = 0  \\
\end{align}
$$

$$
\begin{align}
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
\end{align}
$$

$$
\begin{align}
Δ &= e_{1} + e_{2} + e_{3} + e_{4} + e_{5} + e_{6} + e_{7}  \\
\rho_{C1}^\prime &= \rho_{C1}  \cdot mt_{prot} / \Delta  \\
J_{Hres}^{C1} &= 2\rho_{C1}^\prime (e_{6}a_{61} - e_{1}a_{16}) \\
J_{Q}^{C1} &= 0.5\rho_{C1}^\prime (e_{4}a_{47} - e_{7}a_{74}) \\
J_{NADH}^{C1} &= 0.5\rho_{C1}^\prime (e_{3}a_{34} - e_{4}a_{43}) \\
J_{ROS}^{C1} &= \rho_{C1}^\prime (e_{4}a_{42} - e_{2}a_{24})  \\
\end{align}
$$

| Parameter      | Value     | Units         | Desc.                                        |
| -------------- | --------- | ------------- | -------------------------------------------- |
| $\rho_{C1}$    | 5         | mM            | Concentration of complex I<br />(Adjustable) |
| $\Delta\Psi_B$ | 50        | mV            | Phase boundary potential                     |
| $k_{12}$       | 6.3396E11 | $Hz/mM^2$     |                                              |
| $k_{21}$       | 5         | Hz            |                                              |
| $k_{56}$       | 100       | Hz            |                                              |
| $k_{65}$       | 2.5119E13 | $Hz/mM^2$     |                                              |
| $k_{61}$       | 1E7       | Hz            |                                              |
| $k_{16}$       | 130       | Hz            |                                              |
| $k_{23}$       | 3886.7    | $Hz/mM^{1/2}$ |                                              |
| $k_{32}$       | 9.1295E6  | Hz            |                                              |
| $k_{34}$       | 639.1364  | Hz            |                                              |
| $k_{43}$       | 3.2882    | $Hz/mM^{1/2}$ |                                              |
| $k_{47}$       | 1.5962E7  | Hz/mM         |                                              |
| $k_{74}$       | 65.2227   | Hz            |                                              |
| $k_{75}$       | 24615     | Hz            |                                              |
| $k_{57}$       | 1166.7    | $Hz/mM^{1/2}$ |                                              |
| $k_{42}$       | 6.0318    | Hz/mM         |                                              |
| $E_{FMN}$      | -375      | mV            | Midpoint potential of flavin mononucleotide  |
| $E_{sox}$      | -150      | mV            | Midpoint potential of superoxide             |

### Complex II (Succinate dehydrogenase)

$$
\begin{align}
f_Q &= \frac{[Q]_n}{[Q]_n + [QH_2]_n}  \\
f_{OAA} &= \frac{K_{i, OAA}}{[OAA] + K_{i, OAA}} \\
f_{FUM} &= \frac{K_{i, FUM}}{[FUM] + K_{i, FUM}} \\
f_{SUC} &= \frac{[SUC]}{[SUC] + K_{m, SUC} / f_{OAA} / f_{FUM}} \\
J_{SDH} &= V_{SDH} C2_{inhib} f_{SUC} \frac{f_Q}{f_Q + K_{m, Q}}  \\
J_{c2} &= J_{SDH} \\
\end{align}
$$

| Parameter    | Value | Units       | Desc.                                |
| ------------ | ----- | ----------- | ------------------------------------ |
| $V_{SDH}$    | 250   | mM / minute | Maximum rate of SDH                  |
| $K_{i, OAA}$ | 0.150 | mM          | Inhibition constant for oxaloacetate |
| $K_{m, Q}$   | 0.6   | -           | Michaelis constant for CoQ           |
| $K_{i, FUC}$ | 0.150 | mM          | Inhibition constant for fumarate     |
| $K_{m, SUC}$ | 0.6   | -           | Michaelis constant for succinate     |

### Complex III

$$
\begin{align}
f_{hi} & = [H^+]_{i}  / 10^{-7}M   \\
v_{1} &= v_{Q}^{C1} + v_{Q}^{C2}   \\
v_2 &= k_d([QH_2]_{n} - [QH_2]_{p})  \\
k_{3} &= k_{03}K_{eq3}f_{hi} \\
k_{-3} &= k_{03} \\
v_{3} &= k_3[QH_2]_{p} [FeS]_{ox} - k_{-3}[Q^-]_p [FeS]_{rd}  \\
k_{4, ox} &= k_{04}K_{eq4, ox} \exp(-\alpha\delta_1\Delta\Psi_m F/ RT) \\
k_{4, rd} &= k_{04}K_{eq4, rd} \exp(-\alpha\delta_1\Delta\Psi_m F/ RT) \\
k_{-4, ox} &= k_{04} \exp(\alpha(1-\delta_1)\Delta\Psi_m F/ RT)  \\
k_{-4, rd} &= k_{04} \exp(\alpha(1-\delta_1)\Delta\Psi_m F/ RT)  \\
v_{4, ox} &= k_{4, ox}[Q^-]_p [cytb_1] - k_{-4, ox}[Q]_{p} [cytb_2]  \\
v_{4, rd} &= k_{4, rd}[Q^-]_p [cytb_3] - k_{-4, rd}[Q]_{p} [cytb_4]  \\
v_{5} &= k_d([Q]_{p} - [Q]_{n})  \\
k_{6} &= K_{06}K_{eq6} \exp( -\beta\delta_2\Delta\Psi_m / V_T) \\
k_{-6} &= k_{06}  \exp( \beta(1-\delta_2)\Delta\Psi_m / V_T)  \\
v_{6} &= k_{6} [cytb_2] - k_{-6} [cytb_3]  \\
k_{7, ox} &= k_{07, ox}K_{eq7, ox}\exp(-\gamma\delta_3\Delta\Psi_m F/ RT) \\
k_{7, rd} &= k_{07, rd}K_{eq7, rd}\exp(-\gamma\delta_3\Delta\Psi_m F/ RT) \\
k_{-7, ox} &= k_{07, ox} \exp(\gamma(1-\delta_3)\Delta\Psi_m F/ RT)  \\
k_{-7, rd} &= k_{07, rd} \exp(\gamma(1-\delta_3)\Delta\Psi_m F/ RT)  \\
v_{7, ox} &= (k_{7, ox}[Q]_{n}[cytb_3] - k_{-7, ox}[Q^-]_n [cytb_1])C3_{inhib} \\
v_{7, rd} &= (k_{7, rd}[Q]_{n}[cytb_4] - k_{-7, rd}[Q^-]_n [cytb_2])C3_{inhib} \\
\end{align}
$$

$$
\begin{align}
f_{hm} & = [H^+]_{m}   / 10^{-7}M  \\
k_{8, ox} &= k_{08, ox}K_{eq8, ox}\exp(-\gamma\delta_3\Delta\Psi_m F/ T)(f_{hm})^2  \\
k_{8, rd} &= k_{08,rd}K_{eq8, rd}\exp(-\gamma\delta_3\Delta\Psi_m F/ T)(f_{hm})^2  \\
k_{-8, ox} &= k_{08, ox} \exp(\gamma(1-\delta_3)\Delta\Psi_m F/ RT)  \\
k_{-8, rd} &= k_{08, rd} \exp(\gamma(1-\delta_3)\Delta\Psi_m F/ RT)  \\
v_{8, ox} &= (k_{8, ox}[Q^-]_{n}[cytb_3] - k_{-8, ox}[QH_2]_{n}[cytb_1])C3_{inhib} \\
v_{8, rd} &= (k_{8, rd}[Q^-]_{n}[cytb_4] - k_{-8, rd}[QH_2]_{n}[cytb_2])C3_{inhib} \\
k_9 &= k_{09}K_{eq9}  \\
k_{-9} &= k_{09} \\
v_{9} &= k_{9}[FeS]_{rd}[cytc1]_{ox} - k_{-9}[FeS]_{ox}[cytc1]_{rd}\\
k_{10} &= k_{010}K_{eq10}  \\
k_{-10} &= k_{010} \\
v_{10} &= k_{10}[Q^-]_p[O_2] - k_{-10}[Q]_p[O_2^-]  \\
v_{10b} &= v_{10}  \\
v_{33} &= k_{33}(K_{eq}[cytc1]_{rd}[cytc]_{ox} - [cytc]_{rd}[cytc1]_{ox})  \\
\rho_{C3}^{\prime} &= \rho_{C3} \cdot mt_{prot}  \\
\rho_{C4}^{\prime} &= \rho_{C4} \cdot mt_{prot}  \\
FeS_{rd} &= \rho_{C3}^{\prime} - FeS_{ox} \\
cytc1_{rd} &= \rho_{C3}^{\prime} - cytc1_{ox}  \\
cytc_{rd} &= \rho_{C4}^{\prime} - cytc_{ox}  \\
 [cytb_4] &= \rho_{C3}^{\prime} - [cytb_1] - [cytb_2] - [cytb_3]  \\
 [QH_2]_p &= \Sigma Q - [Q]_n - [Q]_p - [QH_2]_n - [Q^-]_p - [Q^-]_n  \\
J_{hRes}^{C3} &= 2v_{3}     \\
J_{ROS, m}^{C3} &= v_{10}   \\
J_{ROS, i}^{C3} &= v_{10b}  \\
\end{align}
$$

$$
\begin{align}
\frac{d[Q]_n}{dt} &= v_5 - v_{7,ox}- v_{7,rd} - v_1  \\
\frac{d[Q^-]_n}{dt} &= v_{7,ox} + v_{7,rd} - v_{8,ox}- v_{8,rd}  \\
\frac{d[QH_2]_n}{dt} &= v_{8,ox} + v_{8,rd} + v_1 - v_2   \\
\frac{d[QH_2]_p}{dt} &= v_2 -v_3 \\
\frac{d[Q^-]_p}{dt} &= v_3 - v_{10} - v_{10b} - v_{4,ox} - v_{4,rd}   \\
\frac{d[Q]_p}{dt} &= v_{10} + v_{10b} + v_{4,ox} + v_{4,rd} - v_5   \\
\frac{d[cytb_1]}{dt} &= v_{7,ox} + v_{8,ox} - v_{4,ox}    \\
\frac{d[cytb_2]}{dt} &= v_{4,ox} + v_{7,rd} - v_{8,rd} - v_6   \\
\frac{d[cytb_3]}{dt} &= v_6 - v_{4,rd} + v_{7,ox} - v_{8,ox}    \\
\frac{d[cytb_4]}{dt} &= v_{4,rd} - v_{7,rd} - v_{8,rd}   \\
\frac{d[FeS]_{ox}}{dt} &= v_9 - v_3      \\
\frac{d[cytc1]_{ox}}{dt} &= v_{33} - v_9   \\
\frac{d[cytc]_{ox}}{dt} &= V_e - v_{33}   \\
\end{align}
$$

| Parameter    | Value    | Unit  | Desc.                                                    |
| ------------ | -------- | ----- | -------------------------------------------------------- |
| $k_{03}$     | 1,666.63 | Hz/mM | Reverse rate constant for reaction 3                     |
| $K_{eq3}$    | 0.6877   | -     | Equilibrium constant for reaction 3                      |
| $k_{04}$     | 60.67    | Hz/mM | Reverse rate constant for reaction 4                     |
| $K_{eq4,ox}$ | 129.9853 | -     | Equilibrium constant for reaction 4 <br />(bH oxidized)  |
| $K_{eq4,rd}$ | 13.7484  | -     | Equilibrium constant for reaction 4 <br />(bH reduced)   |
| $\delta_1$   | 0.5      | -     |                                                          |
| $\alpha$     | 0.2497   | -     |                                                          |
| $k_d$        | 22000    | Hz    | Diffusion rate of ubiquinone across the membrane         |
| $k_{06}$     | 166.67   | Hz/mM | Reverse rate constant for reaction 6                     |
| $K_{eq6}$    | 9.4596   | -     | Equilibrium constant for reaction 6                      |
| $\delta_2$   | 0.5      | -     |                                                          |
| $\beta$      | 0.5006   | -     |                                                          |
| $k_{07,ox}$  | 13.33    | Hz/mM | Reverse rate constant for reaction 7 <br />(bL oxidized) |
| $K_{eq7,ox}$ | 3.0748   | -     | Equilibrium constant for reaction 7 <br />(bL oxidized)  |
| $k_{07,rd}$  | 1.667    | Hz/mM | Reverse rate constant for reaction 7 <br />(bL reduced)  |
| $K_{eq7,rd}$ | 29.0714  | -     | Equilibrium constant for reaction 7 <br />(bL reduced)   |
| $\delta_3$   | 0.5      | -     |                                                          |
| $\gamma$     | 0.2497   | -     | $\alpha + \beta + \gamma = 1$                            |
| $k_{08,ox}$  | 83.33    | Hz/mM | Reverse rate constant for reaction 8 <br />(bL oxidized) |
| $K_{eq8,ox}$ | 129.9853 | -     | Equilibrium constant for reaction 8 <br />(bL oxidized)  |
| $k_{08,rd}$  | 8.333    | Hz/mM | Reverse rate constant for reaction 8 <br />(bL reduced)  |
| $K_{eq8,rd}$ | 9.4596   | -     | Equilibrium constant for reaction 8 <br />(bL reduced)   |
| $k_{09}$     | 833      | Hz/mM | Reverse rate constant for reaction 9                     |
| $K_{eq9}$    | 0.2697   | -     | Equilibrium constant for reaction 9                      |
| $k_{010}$    | 0.8333   | Hz/mM | Reverse rate constant for reaction 10                    |
| $K_{eq10}$   | 1.4541   | -     | Equilibrium constant for reaction 10                     |
| $k_{33}$     | 2469.13  | Hz/mM | Reverse rate constant for reaction 33                    |
| $K_{eq33}$   | 2.1145   | -     | Equilibrium constant for reaction 33                     |
| $\rho_{C3}$  | 0.325    | mM    | Complex III content                                      |


### Complex IV

$$
\begin{align}
f_{H_{m}} &= \exp(-\delta_5\Delta\Psi_m F/ RT) ([H^+]_m /10^{-7}M) \\
f_{H_{i}} &= \exp((1-\delta_5)\Delta\Psi_m F/ RT) ([H^+]_i /10^{-7}M)  \\
f_{C_{rd}} &= [cytc]_{rd} \\
f_{C_{ox}} &= \text{exp}((1-\delta_5)\Delta\Psi_m F/ RT) [cytc]_{ox} \\
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
\end{align}
$$

| Parameter     | Value     | Unit    | Desc.                  |
| ------------- | --------- | ------- | ---------------------- |
| $\Sigma cytc$ | 0.325     | mM      | Cytochrome c pool      |
| $\rho_{C4}$   | 0.325     | mM      | Complex IV content     |
| $k_{34}$      | 2.9445E10 | Hz/mM^3 | Rate constant @ pH = 7 |
| $k_{-34}$     | 290.03    | Hz/mM^3 | Rate constant @ pH = 7 |
| $k_{35}$      | 45000     | Hz/mM   |                        |
| $k_{36}$      | 4.826E11  | Hz/mM   | Rate constant @ pH = 7 |
| $k_{-36}$     | 4.826     | Hz/mM   | Rate constant @ pH = 7 |
| $k_{37}$      | 1.7245E8  | Hz      | Rate constant @ pH = 7 |
| $k_{-37}$     | 17.542    | Hz      | Rate constant @ pH = 7 |

### Complex V (ATP synthase)

$$
\begin{align}
J_{F1Fo} &= -\rho^{F1} ((100 p_a + p_{c1} v_B) v_a - (p_a + p_{c2} v_a) v_h)  / \Delta \\
J_H^{F1Fo} &= -3\rho^{F1} (100p_a(1 + v_a) - (p_a + p_b)v_h) / \Delta \\
\Delta &= (1 + p_1 v_a)v_B + (p_2 + p_3 v_a)v_h \\
v_B &= \text{exp}(3\Delta\Psi_B / V_T)   \\
v_h &= \text{exp}(3\Delta p / V_T)  \\
v_a &= \frac{K_{eq}^{'} \cdot \Sigma[ATP]_m}{ \Sigma[Pi]_m \cdot \Sigma[ADP]_m } \\
\end{align}
$$

| Parameter      | Value     | Unit | Desc.                                                          |
| -------------- | --------- | ---- | -------------------------------------------------------------- |
| $\rho_{F1}$    | 5         | mM   | Concentration of F1-Fo ATPase                                  |
| $K_{eq}^{'}$   | 6.47E5    | M    | Apparent equilibrium constant for ATP hydrolysis[^Golding1995] |
| $\Delta\Psi_B$ | 50        | mV   | Phase boundary potential                                       |
| $p_{a}$        | 1.656E-5  | Hz   | Sum of products of rate constants                              |
| $p_{b}$        | 3.373E-7  | Hz   | Sum of products of rate constants                              |
| $p_{c1}$       | 9.651E-14 | Hz   | Sum of products of rate constants                              |
| $p_{c2}$       | 4.585E-14 | Hz   | Sum of products of rate constants                              |
| $p_{1}$        | 1.346E-4  | -    | Sum of products of rate constants                              |
| $p_{2}$        | 7.739E-7  | -    | Sum of products of rate constants                              |
| $p_{3}$        | 6.65E-15  | -    | Sum of products of rate constants                              |

[^Golding1995]: Golding, E. M., Teague, W. E., & Dobson, G. P. (1995). Adjustment of K’ to varying pH and pMg for the creatine kinase, adenylate kinase and ATP hydrolysis equilibria permitting quantitative bioenergetic assessment. The Journal of Experimental Biology, 198(Pt 8), 1775–1782.

## Reactive oxygen species (ROS) scavenging and transport

### Catalase (CAT)

Includes inhibition by high levels of hydrogen peroxide

$$
\begin{align}
V_{CAT} = 2k_1E_T[H_2O_2]_i \cdot e^{-fr[H_2O_2]_i} \\
\end{align}
$$

| Parameter | Value | Unit      | Desc.                                  |
| --------- | ----- | --------- | -------------------------------------- |
| $k_1$     | 17    | 1/(mM*ms) | Rate constant of catalase              |
| $E_T$     | 0.01  | mM        | Extra-matrix concentration of catalase |
| $fr$      | 0.05  | 1/mM      | Hydrogen peroxide inhibition factor    |

### Superoxide dismutase (SOD)

Based on (McADAM, 1976) model.

$$
\begin{align}
J_{SOD} &= \frac{2 k_5 E_T f_{sox} (k_1 + k_3^\prime)}{ k_5 (2 k_1 + k_3^\prime) + k_3^\prime f_{sox}} \\
k_3^\prime &= k_3 (1 + \frac{[H_2O_2]}{K_{H_2O_2}})  \\
f_{sox} &= k_1^{SOD} [O_2^-]
\end{align}
$$

| Parameter | Value | Unit      | Desc.                                 |
| --------- | ----- | --------- | ------------------------------------- |
| $k_1$     | 1200  | 1/(mM*ms) | Rate constant for EA -> EB            |
| $k_3$     | 24    | 1/(mM*ms) | Rate constant for EB -> EC            |
| $k_5$     | 0.24  | 1/s       | Rate constant for EC -> EA            |
| $K_{i}$   | 500   | μM        | Inhibition constant for H2O2          |
| $E_{T}$   | 3     | μM        | Concentration of Cu,ZnSOD (cytosolic) |

### Glutathione peroxidase (GPX)

Dalziel type Ping-pong mechanism.

$$
\begin{align}
J_{GPX} &= \frac{E_T}{A + B}      \\
A &= \frac{\Phi_1}{[H_2O_2] }  \\
B & = \frac{\Phi_2}{[GSH] }  \\
\end{align}
$$

| Parameter | Value | Unit | Desc.               |
| --------- | ----- | ---- | ------------------- |
| $E_T$     | 10    | μM   | GPX content         |
| $\Phi_1$  | 5     | mM/s | Dalziel coefficient |
| $\Phi_2$  | 75    | mM/s | Dalziel coefficient |

### Glutathione reductase (GR)

Michaelis-Menten kinetics.

$$
\begin{align}
J_{GR} &= k_1^{GR} E_T \frac{[GSSG]}{[GSSG] + K_{GSSG}} \frac{[NADPH]}{[NADPH] + K_{NADPH}} \\
\Sigma [GSH] &= [GSH] + 2 [GSSG]
\end{align}
$$

| Parameter      | Value | Unit | Desc.                        |
| -------------- | ----- | ---- | ---------------------------- |
| $E_T$          | 10    | μM   | GR content (cytosolic)       |
| $k_1^{GR}$     | 5     | Hz   | Catalytic constant of GR     |
| $K_{GSSG}$     | 60    | μM   | Michaelis constant for GSSG  |
| $K_{NADPH}$    | 15    | μM   | Michaelis constant for NADPH |
| $\Sigma [GSH]$ | 1     | mM   | Cytosolic GSH pool           |

### Inner mitochondrial anion channel (IMAC)

$$
\begin{align}
g_{IMAC} &= \left( a + b \frac{[O_2^-]_i}{[O_2^-]_i + K_{CC}} \right) \left( G_L + \frac{G_{max}}{1 + e^{κ(\Delta\Psi_m^b + \Delta\Psi_m)}} \right) \\
V_{IMAC} &= g_{IMAC}\Delta\Psi_m \\
V_{tr}^{ROS} &= j \cdot g_{IMAC} \left( \Delta\Psi_m + V_T ln \left( \frac{[O_2^-]_m}{[O_2^-]_i} \right) \right) \\
\end{align}
$$

| Parameter        | Value  | Unit         | Desc.                             |
| ---------------- | ------ | ------------ | --------------------------------- |
| a                | 0.001  | -            | Basal IMAC conductance            |
| b                | 10000  | -            | Activation factor by superoxide   |
| $K_{CC}$         | 10     | μM           | Activation constant by superoxide |
| $G_L$            | 0.035  | μM * Hz / mV | Integral conductance for IMAC     |
| $G_{max}$        | 3.9085 | μM * Hz / mV | Leak conductance of IMAC          |
| $\kappa$         | 0.07   | 1/mV         | Steepness factor                  |
| $\Delta\Psi_m^b$ | 4      | mV           | Potential at half saturation      |
| j                | 0.1    | -            | Fraction of IMAC conductance      |

### ODEs for ROS transport and scavenging

$$
\begin{align}
\frac{d [ O_{2}^{-}]_{m}}{dt} &= J_{ROS,m} - J^{Tr}_{ROS}  \\
\frac{d [ O_{2}^{-}]_{i}}{dt} &= \frac{V_{mito}}{V_{cyto}} J^{Tr}_{ROS} -J_{SOD,i}  \\
\frac{d[H_2O_2]_i}{dt} &= 0.5J_{SOD,i}  -J_{GPX,i} - J_{CAT}  \\
\frac{d[GSH]_i}{dt} &= J_{GR,i} - J_{GPX,i} \\
\end{align}
$$


[^Cortassa2004]:Cortassa S, Aon MA, Winslow RL, O'Rourke B. A mitochondrial oscillator dependent on reactive oxygen species. Biophys J. 2004;87(3):2060-73. [PMC1304608](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1304608/)

[^Kembro2013]:Kembro JM, Aon MA, Winslow RL, O'Rourke B, Cortassa S. Integrating mitochondrial energetics, redox and ROS metabolic networks: a two-compartment model. Biophys J. 2013;104(2):332-43. [PMC3552263](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3552263/)

## Mitochondrial ion transport

### Adenine Nucleotide translocator (ANT)

$$
\begin{align}
J_{ANT} &= V_{max}^{ANT}\frac{AB - \delta PQ}{(B + \delta^{h_{ANT}} P)(A + Q)}  \\
A &= [ATP^{4-}]_m = 0.025 [ATP]_m  \\
B &= [ADP^{3-}]_i = 0.45 [ADP]_i \\
P &= [ATP^{4-}]_i = 0.25 [ATP]_i \\
Q &= [ADP^{3-}]_m = 0.17 [ADP]_m \\
\delta &= \text{exp}(-\Delta\Psi_m F / RT)
\end{align}
$$


| Parameter       | Value | Unit | Desc.               |
| --------------- | ----- | ---- | ------------------- |
| $V_{max}^{ANT}$ | 5     | mM/s | Maximal rate of ANT |
| $h_{ANT}$       | 0.5   | -    | Fraction of MMP     |

### Mitochondrial calcium uniporter (MCU)

$$
\begin{align}
J_{uni} &= V_{max}^{Uni} \frac{S (1+S)^3}{(1+S)^4 + L(1 + A)^n} \frac{\delta}{e^\delta-1}  \\
S &= [Ca^{2+}]_i / K_{trans}  \\
A &= [Ca^{2+}]_i / K_{act}    \\
\delta &= -2 (\Delta\Psi_m - \Delta\Psi_0) F/RT \\
\end{align}
$$

| Parameter       | Value | Unit | Desc.                              |
| --------------- | ----- | ---- | ---------------------------------- |
| $V_{max}^{Uni}$ | 4.46  | mM/s | Maximal rate                       |
| $\Delta\Psi_0$  | 91    | mV   | Offset potential                   |
| $K_{act}$       | 0.38  | μM   | Activation constant for calcium    |
| $K_{trans}$     | 19    | μM   | Dissociation constant for calcium  |
| n               | -2.8  | -    | Activation cooperativity           |
| L               | 110   | -    | Keq for conformational transitions |

### Mitochondrial sodium-calcium exchanger (NCLX)

$$
\begin{align}
J_{NCLX} = V_{max}^{NCLX} \exp(b\Delta\Psi_m F/RT) \frac{[Ca^{2+}]_m}{[Ca^{2+}]_i} \left( \frac{[Na^+]_i}{[Na^+]_i + K_{Na}^{NCLX}} \right)^n \frac{[Ca^{2+}]_m}{[Ca^{2+}]_m + K_{Ca}^{NCLX}}
\end{align}
$$

| Parameter        | Value   | Unit | Desc.                             |
| ---------------- | ------- | ---- | --------------------------------- |
| $V_{max}^{NCLX}$ | 0.04665 | mM/s | Maximal rate of NCLX              |
| b                | 0.5     | -    | Fraction of MMP                   |
| $K_{Na}^{NCLX}$  | 9.4     | mM   | Dissociation constant for sodium  |
| $K_{Ca}^{NCLX}$  | 0.375   | μM   | Dissociation constant for calcium |
| $n$              | 3       | -    | Cooperativity                     |

### Mitochondrial proton leak

$$
J_{hleak} = g_H\Delta\Psi_m
$$

| Parameter | Value | Unit            | Desc.                                                 |
| --------- | ----- | --------------- | ----------------------------------------------------- |
| $g_{H}$   | 2     | mM / (Volt * s) | Ionic conductance of the inner mitochondrial membrane |

### ODEs for mitochondrial ion transport

$$
\begin{align}
\frac{d [Ca^{2+}]_m}{dt} &=\delta_{Ca}( J_{uni} - J_{NCLX}) \\
\frac{d [Na^+]_m}{dt} &= J_{NCLX} - J_{NaH} \\
C_{m}\frac{d \Delta \Psi_m}{dt} &= J_{Hres} - J_{Hu} - J_{ANT} - J_{Hleak} -J_{NCLX} - J_{uni} - J_{IMAC} \\
\end{align}
$$

## General parameters

| Parameter       | Value                  | Unit                  | Desc.                                    |
| --------------- | ---------------------- | --------------------- | ---------------------------------------- |
| F               | 96485                  | C/mol                 | Faraday constant                         |
| T               | 310                    | K                     | Absolute temperature                     |
| R               | 8.314                  | J/molK                | Universal gas constant                   |
| $V_T$           | 26.71                  | mV                    | Thermal voltage (=${RT}/{F}$)            |
| $C_m$           | 1.0                    | $\text{μF/cm}^2$      | Plasma membrane capacitance              |
| $C_{mito}$      | 1.812                  | mM/V                  | Mitochondrial inner membrane capacitance |
| $\delta_{Ca}$   | 0.0003                 | -                     | Mitochondrial free calcium fraction      |
| $\delta_H$      | 1E-5                   | -                     | Mitochondrial proton buffering factor    |
| $V_{myo}$       | $25.84$                | $pL$                  | Cytosolic volume                         |
| $V_{mito}$      | $15.89$                | $pL$                  | Mitochondrial volume                     |
| $V_{NSR}$       | $1.4$                  | $pL$                  | Network SR volume                        |
| $V_{JSR}$       | $0.16$                 | $pL$                  | Junctional SR volume                     |
| $V_{SS}$        | $0.000495$             | $pL$                  | Subspace volume                          |
| $A_{cap}$       | $1.534 \cdot 10^{-4} $ | $cm^{2}$              | Capacitance area                         |
| $C_{m}$         | $1.0$                  | $\mu F \cdot cm^{-2}$ | Plasma membrane capacitance              |
| $[K^+]_{o}$     | $5.4$                  | $mM$                  | Extracellualr potassium                  |
| $[Na^+]_{o}$    | $140$                  | $mM$                  | Extracellualr sodium                     |
| $[Ca^{2+}]_{o}$ | $2$                    | $mM$                  | Extracellualr calcium                    |
| $C_{mito}$      | $1.812 \cdot 10^{-3}$  | $mM/mV$               | Inner membrane capacitance               |
| $g_{H}$         | $1 \cdot 10^{-8}$      | $mM/msmV$             | Inner membrane conductance               |

### Fixed concentrations

| Parameter          | Value   | Unit | Desc.                                    |
| ------------------ | ------- | ---- | ---------------------------------------- |
| $pH_i$             | 7       |      | CytosoliWc pH                            |
| $pH_m$             | 7.3-7.8 |      | Mitochondrial pH                         |
| $[O_2]$            | 0.006   | mM   | Tissue oxygen concentration              |
| $[Mg^{2+}]_i$      | 1.0     | mM   | Cytosolic magnesium concentration        |
| $[Mg^{2+}]_m$      | 0.4     | mM   | Mitochondrial magnesium concentration    |
| $\Sigma[Pi]_m$     | 8.6512  | mM   | Sum of mitochondrial inorganic phosphate |
| $\Sigma{[N]}$      | 1       | mM   | Sum of mitochondrial NAD and NADH        |
| $\Sigma[A]_m$      | 1.5     | mM   | Sum of mitochondrial ATP and ADP         |
| $\Sigma{[NADP]_m}$ | 0.1     | mM   | Sum of mitochondrial NADPH plus NADP     |
| $[Ca^{2+}]_i$      | 1E-4    | mM   | Cytosolic calcium concentration          |
