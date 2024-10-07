# Sarcoplasmic ion currents

## Time-dependent delayed rectifier potassium current (IK)

$$
\begin{aligned}
I_K &= \bar G_K X_1 X_K^2 (V - E_K) \\
E_K &= \frac{RT}{F} \ln \frac{[K^+]_o + P_{Na,K}[Na^+]_o}{ [K^+]_i + P_{Na,K}[Na^+]_i} \\
\bar G_K &= 0.282\sqrt{[K^+]_o /5.4mM} * (mS/cm^2) \\
X_1 &= (1+ e^{(V_m-40)/40})^{-1} \\
\frac{dX_k}{dt} &= \alpha_X - X_k (\alpha_X + \beta_X) \\
\alpha_X &= \frac{V_m+30}{1 - e^{-0.148(V_m+30)}} * 0.0719Hz \\
\beta_X &= \frac{V_m+30}{e^{0.0687(V_m+30)} -1} * 0.131Hz \\
\end{aligned}
$$

## Time-independent potassium current (IK1)

$$
\begin{aligned}
\Delta V &= V_m - E_{K1} \\
I_{K1} &= \bar G_{K1}K_{1 \infty}\Delta V \\
E_{K1} &= \frac{RT}{F} \ln \frac{[K^+]_o}{[K^+]_i}  \\
\bar G_{K1} &= 0.748\sqrt{[K^+]_o / 5.4mM} * (mS/cm^2) \\
K_{1 \infty} &= \frac{\alpha_{K_1}}{\alpha_{K_1} + \beta_{K_1}} \\
\alpha_{K_1} &= \frac{1.02}{1 + e^{0.2385(\Delta V -59.215)}}  * kHz  \\
\beta_{K_1} &= \frac{0.4912e^{0.28032(\Delta V + 5.476)} + e^{0.06175(\Delta V -594.31)}}{1 + e^{-0.5143(\Delta V + 4.753)}} * kHz
\end{aligned}
$$

## Plateau potassium current (IKp)

$$
\begin{aligned}
E_{Kp} &= \frac{RT}{F} \ln \frac{[K^+]_o}{[K^+]_i} \\
I_{Kp} &= \frac{\bar G_{Kp} (V - E_{Kp})}{1 + e^{(7.488-V_m) / 5.98}} \\
\end{aligned}
$$

## Fast Na current (INa)

$$
\begin{aligned}
I_{Na} &= \bar G_{Na} m_{Na}^{3} h_{Na} j_{Na} (V_m-E_{Na}) \\
E_{Na} &= \frac{RT}{F} \ln \frac{[Na^{+}]_o}{[Na^{+}]_i} \\
\frac{dm_{Na}}{dt} &= \alpha_{m} - m_{Na}(\alpha_{m} + \beta_{m}) \\
\frac{dh_{Na}}{dt} &= \alpha_{h} - h_{Na}(\alpha_{h} + \beta_{h}) \\
\frac{dj_{Na}}{dt} &= \alpha_{j} - m_{Na}(\alpha_{j} + \beta_{j}) \\
\alpha_{m} &= 0.32kHz \frac{V + 47.13}{1 - e^{-0.1(V_m+47.13)}} \\
\beta_{m} &= 0.08kHz * e^{-V_m / 11} \\
\\
For \ V & \ge -40mV \\
\alpha_{h} &= \alpha_{j} = 0 \\
\beta_{h} &= (0.13 ms (1+e^{-(V_m+10.66)/11.1}))^{-1} \\
\beta_{j} &= 0.3kHz\frac{e^{-2.535  \cdot 10^{-7}V_m}}{1 + e^{-0.1(V_m + 32)}} \\
\\
For \ V & < -40mV \\
\alpha_{h} &= 0.135kHz * e^{-(V_m+80)/6.8} \\
\alpha_{j} &= (-127140e^{0.2444 V_m}-3.474 \cdot 10^{-5}e^{-0.04391 V_m})\frac{V_m + 37.78}{1+e^{0.311( V_m +79.23)}} * kHz \\
\beta_{h} &= (3.56e^{0.079 V_m} + 3.1  \cdot 10^{5}e^{0.35 V_m}) * kHz \\
\beta_{j} &= \frac{0.1212e^{-0.01052 V_m}}{1+e^{-0.1378(V_m + 40.14)}} * kHz \\
\end{aligned}
$$

## Sodium-calcium exchanger current (INaCa)

$$
\begin{aligned}
I_{NaCa} &= k_{NaCa}  \cdot f_{Nao }  \cdot f_{Cao }\frac{exp(V_mF/RT)\phi_{Na}^3 - \phi_{Ca}}{exp((1 - \eta) V_mF/RT ) + k_{sat}} \\
f_{Nao} &= \frac{([Na^+]_o)^3}{([Na^+]_o)^3 + (K_{M,Na}^{NaCa})^3} \\
f_{Cao} &= \frac{[Ca^+]_o}{[Ca^+]_o + K_{M,Ca}^{NaCa}} \\
\phi_{Na} &= \frac{[Na^+]_i}{ [Na^+]_o} \\
\phi_{Ca} &= \frac{[Ca^{2+}]_i}{[Ca^{2+}]_o} \\
\end{aligned}
$$

## Background calcium ($I_{Ca,b}$) and sodium currents ($I_{Na,b}$)

$$
\begin{aligned}
I_{Ca,b} &= \bar G_{Ca,b} (V_m - \frac{RT}{2F} \ln \frac{[Ca^{2+}]_o}{[Ca^{2+}]_i})  \\
I_{Na,b} &= \bar G_{Na,b} (V_m - \frac{RT}{F} \ln \frac{[Na^{+}]_o}{[Na^{+}]_i}) \\
\end{aligned}
$$

## Non-specific calcium-activated current (InsCa)

$$
\begin{aligned}
f_{Ca} &= \frac{([Ca^{2+}]_i)^3}{([Ca^{2+}]_i)^3 + (K_{m}^{nsCa})^3}\\
I_{nsNa} &= 0.75  \cdot f_{Ca}  \cdot  \Phi_{Na}(P_{nsNa}, z_{Na}, V_m, [Na^+]_i, [Na^+]_o)  \\
I_{nsK} &= 0.75  \cdot f_{Ca}  \cdot  \Phi_{K}(P_{nsK}, z_{K}, V_m, [K^+]_i, [K^+]_o)  \\
\end{aligned}
$$

## Sodium-potassium ATPase current (INaK)

The Na+/K+ ATPase activity depends on the ATP concentration, as well as the competitive inhibition by ADP.

$$
\begin{aligned}
I_{NaK} &= \bar I_{NaK}  \cdot f_{ATP}  \cdot f_{Na}  \cdot f_{K} \cdot f_{NaK}  \\
\sigma &= \frac{e^{[Na^+]_o / 67.3mM}-1}{7}  \\
f_{NaK} &= (1 + 0.1245 \cdot \exp(-0.1V_m F / RT) + 0.0365 \sigma \cdot \exp(-V_m F / RT))^{-1}  \\
f_{Na} &= \frac{([Na^+]_i)^{1.5}}{([Na^+]_i)^{1.5} + (K_{m, Na_i})^{1.5}} \\
f_{K} &= \frac{[K^+]_o}{[K^+]_o + K_{m, K_o}} \\
f_{ATP} &= \frac{[ATP]_i}{[ATP]_i + K_{M,ATP}^{NaK} / f_{ADP}} \\
f_{ADP} &= \frac{K_{i,ADP}^{NaK}}{K_{i,ADP}^{NaK} + [ADP]_i} \\
\end{aligned}
$$

## Electrophysiology ODEs

$$
\begin{aligned}
\frac{d[Na^+]_i}{dt} &= -(I_{Na} + 3I_{NaCa} + 3I_{NaK})\frac{A_{cap}}{V_{myo}F} + (V_{NHE} - 3V_{NaCa}) \frac{V_{mito}}{V_{myo}} \\
\frac{d[K^+]_i}{dt} &= -(I_{Ks} + I_{Kr} + I_{K1} + I_{Kp} + I_{Ca,K}-2I_{NaK})\frac{A_{cap}}{V_{myo}F} \\
C_m\frac{dV_m}{dt} &= -(I_{Na} + I_{CaL} + I_{Kr} + I_{Ks} + I_{K1} + I_{Kp} + I_{NaCa} + I_{NaK} + I_{pCa} + I_{Ca, b} + I_{K_{ATP}} + I_{stim}) \\
\end{aligned}
$$

## GHK current equation

$$
\Phi_s(P_s, z_s, V_m, [S]_i, [S]_o) := P_sz^2_s\frac{V_mF^2}{RT}\frac{[S]_i - [S]_o\exp(-z_sV_mF/RT)}{1-\exp(-z_sV_mF/RT)}
$$

## Parameters

| Symbol          | Value                | Units                 | Description                                           |
| --------------- | -------------------- | --------------------- | ----------------------------------------------------- |
| $G_{Na}$        | $12.8$               | $mS/cm^2$    | Maximal Na channel conductance                        |
| $G_{Kp}$        | $0.00828$            | $mS/cm^2$    | Maximal plateau K channel conductance                 |
| $G_{K,0}$       | $0.282$              | $mS/cm^2$    | IK conductance                                        |
| $G_{K1,0}$      | $0.748$              | $mS/cm^2$    | IK1 conductance                                       |
| $P_{NaK}$       | $0.01833$            |                       | Na+ permeability ratio of K+ channel                  |
| $K_{NaCa}$      | $9000$               | $\mu A/cm^2$ | NCX current                                           |
| $K_{Na}^{NCX}$  | $87.5$               | $mM$                  | Dissociation constant of sodium for NCX               |
| $K_{Ca}^{NCX}$  | $1.38$               | $mM$                  | Dissociation constant of calcium for NCX              |
| $K_{sat}^{NCX}$ | $0.1$                |                       | NCX saturation factor at negative potentials          |
| $\eta^{NCX}$    | $0.35$               |                       | Voltage dependence of NCX                             |
| $P_{ns,Na}$     | $1.75 \cdot 10^{-7}$ | $cm/s$         | Nonspecific channel current Na permeability           |
| $P_{ns,K}$      | $0$                  | $cm/s$         | Nonspecific channel current K permeability            |
| $K_{ca}^{ns}$   | $1.2$                | $\mu M$               | Ca2+ half-saturation constant for nonspecific current |
| $G_{Ca,b}$      | $0.003217$           | $mS/cm^2$    | Maximum background current Ca2+ conductance           |
| $G_{Na,b}$      | $0.003217$           | $mS/cm^2$    | Maximum background current Na+ conductance            |
