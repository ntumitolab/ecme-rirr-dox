# Reactive oxygen species (ROS) scavenging and transport

## Catalase (CAT)

Includes inhibition by high levels of hydrogen peroxide

$$
V_{CAT} = 2k_1E_T[H_2O_2]_i \cdot e^{-fr[H_2O_2]_i} \\
$$

| Parameter | Value | Unit  | Desc.                                  |
| --------- | ----- | ----- | -------------------------------------- |
| $k_1$     | 17    | 1/(mM*ms) | Rate constant of catalase          |
| $E_T$     | 0.01  | mM    | Extra-matrix concentration of catalase |
| $fr$      | 0.05  | 1/mM  | Hydrogen peroxide inhibition factor    |

## Superoxide dismutase (SOD)

Based on (McADAM, 1976) model.

$$
\begin{aligned}
J_{SOD} &= {2k_5E_Tf_{sox}(k_1 + k_3^\prime) \over k_5(2 k_1 + k_3^\prime) + k_3^\prime f_{sox}}   \\
k_3^\prime &= k_3 (1 + \frac{[H_2O_2]}{K_{H_2O_2}})  \\
f_{sox} &= k_1^{SOD} [O_2^-]
\end{aligned}
$$

| Parameter | Value   | Unit  | Desc.                                  |
| --------- | ------- | ----- | -------------------------------------- |
| $k_1$     | 1200    | 1/(mM*ms) | Rate constant for EA -> EB         |
| $k_3$     | 24      | 1/(mM*ms) | Rate constant for EB -> EC         |
| $k_5$     | 0.24    | 1/s    | Rate constant for EC -> EA            |
| $K_{i}$   | 500     | μM    | Inhibition constant for H2O2           |
| $E_{T}$   | 3       | μM    | Concentration of Cu,ZnSOD (cytosolic)  |

## Glutathione (GSH) systems

### Glutathione peroxidase (GPX)

Dalziel type Ping-pong mechanism.

$$
\begin{aligned}
J_{GPX} &= \frac{E_T}{A + B}      \\
A &= \frac{\Phi_1}{[H_2O_2] }  \\
B & = \frac{\Phi_2}{[GSH] }  \\
\end{aligned}
$$

| Parameter | Value  | Unit | Desc.                       |
| --------- | ------ | ---- | --------------------------- |
| $E_T$     | 10     | μM   | GPX content                 |
| $\Phi_1$  | 5      | mM/s | Dalziel coefficient         |
| $\Phi_2$  | 75     | mM/s | Dalziel coefficient         |

### Glutathione reductase (GR)

Michaelis-Menten kinetics.

$$
J_{GR} = k_1^{GR} E_T \frac{[GSSG]}{[GSSG] + K_{GSSG}} \frac{[NADPH]}{[NADPH] + K_{NADPH}}
$$

| Parameter   | Value | Unit | Desc.                        |
| ----------- | ----- | ---- | ---------------------------- |
| $E_T$       | 10    | μM   | GR content (cytosolic)       |
| $k_1^{GR}$  | 5     | Hz   | Catalytic constant of GR     |
| $K_{GSSG}$  | 60    | μM   | Michaelis constant for GSSG  |
| $K_{NADPH}$ | 15    | μM   | Michaelis constant for NADPH |

### Conservation relationship of glutathione

$$
\Sigma [GSH] = [GSH] + 2 [GSSG]
$$

| Parameter | Value | Unit | Desc.                  |
| --------- | ----- | ---- | ---------------------- |
| $\Sigma [GSH]$ |   1    | mM   | Cytosolic GSH pool     |

## Inner mitochondrial anion channel (IMAC)

$$
\begin{aligned}
g_{IMAC} &= \left( a + b \frac{[O_2^-]_i}{[O_2^-]_i + K_{CC}} \right) \left( G_L + \frac{G_{max}}{1 + e^{κ(\Delta\Psi_m^b + \Delta\Psi_m)}} \right) \\
V_{IMAC} &= g_{IMAC}\Delta\Psi_m \\
V_{tr}^{ROS} &= j \cdot g_{IMAC} \left( \Delta\Psi_m + V_T ln \left( \frac{[O_2^-]_m}{[O_2^-]_i} \right) \right) \\
\end{aligned}
$$

| Parameter        | Value  | Unit         | Desc.                              |
| ---------------- | ------ | ------------ | ---------------------------------- |
| a                | 0.001  | -            | Basal IMAC conductance             |
| b                | 10000  | -            | Activation factor by superoxide   |
| $K_{CC}$         | 10     | μM           | Activation constant by superoxide |
| $G_L$            | 0.035  | μM * Hz / mV | Integral conductance for IMAC      |
| $G_{max}$        | 3.9085 | μM * Hz / mV | Leak conductance of IMAC           |
| $\kappa$         | 0.07   | 1/mV         | Steepness factor                   |
| $\Delta\Psi_m^b$ | 4      | mV           | Potential at half saturation       |
| j                | 0.1    | -            | Fraction of IMAC conductance       |

## ODE system for ROS transport and scavenging

$$
\begin{aligned}
\frac{d [ O_{2}^{-}]_{m}}{dt} &= J_{ROS,m} - J^{Tr}_{ROS}  \\
\frac{d [ O_{2}^{-}]_{i}}{dt} &= \frac{V_{mito}}{V_{cyto}} J^{Tr}_{ROS} -J_{SOD,i}  \\
\frac{d[H_2O_2]_i}{dt} &= 0.5J_{SOD,i}  -J_{GPX,i} - J_{CAT}  \\
\frac{d[GSH]_i}{dt} &= J_{GR,i} - J_{GPX,i} \\
\end{aligned}
$$


[^Cortassa2004]:Cortassa S, Aon MA, Winslow RL, O'Rourke B. A mitochondrial oscillator dependent on reactive oxygen species. Biophys J. 2004;87(3):2060-73. [PMC1304608](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1304608/)

[^Kembro2013]:Kembro JM, Aon MA, Winslow RL, O'Rourke B, Cortassa S. Integrating mitochondrial energetics, redox and ROS metabolic networks: a two-compartment model. Biophys J. 2013;104(2):332-43. [PMC3552263](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3552263/)
