# General things

## Functions

Nernst potential

$$
E_N(X_o, X_i, z_X) := \frac{RT}{Fz_X} ln(\frac{X_o}{X_i}) \approx \frac{26.7mV}{z_X} ln(\frac{X_o}{X_i}) 
$$

Hill function

$$
Hill(x, k, n) := \frac{x^n}{x^n + k^n}
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


## Fixed concentrations

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

## Initial conditions of state variables

| State variable | Value   | Unit |
| -------------- | ------- | ---- |
| $[Ca^{2+}]_m$  | 0.02738 | μM   |
| $[ADP]_m$      | 15.8    | μM   |
| $\Delta\Psi_m$ | 193     | mV   |
| $[NADH]$       | 965     | μM   |
| $[H^+]_m$      | 0.0697  | μM   |
| $[Pi]_m$       | 8280    | μM   |
| $[ISOC]$       | 121     | μM   |
| $[\alpha KG]$  | 130     | μM   |
| $[SCoA]$       | 16.1    | μM   |
| $[SUC]$        | 37      | μM   |
| $[FUM]$        | 235     | μM   |
| $[MAL]$        | 228     | μM   |
| $[OAA]$        | 1.28    | μM   |
| $[Na^+]_m$     | 98.5    | μM   |
| $[O_2^-]_m$    | 6.39E-3 | μM   |
| $[O_2^-]_i$    | 4.83E-5 | μM   |
| $[H_2O_2]_m$   | 0.0823  | μM   |
| $[H_2O_2]_i$   | 2.83E-4 | μM   |
| $[GSH]_m$      | 1650    | μM   |
| $[GSH]_i$      | 1650    | μM   |
| $[TrxSH_2]_m$  | 24.3    | μM   |
| $[TrxSH_2]_i$  | 49.9    | μM   |
| $[PSSG]_m$     | 0.676   | μM   |
| $[PSSG]_i$     | 0.0264  | μM   |

## Acid-base equilibria and binding polynomials[^Wei2011]

For both cytoplasmic and mitochondrial compartments.

$$
\begin{aligned}
P_{ATP} &= 1 + \frac{[H^+]}{K_{a}^{ATP}} + \frac{[Mg^{2+}]}{K_{Mg}^{ATP}}  \\
P_{ADP} &= 1 + \frac{[H^+]}{K_{a}^{ADP}} + \frac{[Mg^{2+}]}{K_{Mg}^{ADP}}  \\
P_{Pi} &= 1 + \frac{[H^+]}{K_{a}^{Pi}}      \\
P_{SUC} &= 1 + \frac{[H^+]}{K_{a}^{SUC}}    \\
P_{H_2O} &= 1 + \frac{[H^+]}{K_w}     \\
[OH^-] &= K_w / [H^+]     \\
[ATP^{4-}] &= \Sigma ATP / P_{ATP}   \\
[ADP^{3-}] &= \Sigma ADP / P_{ADP}   \\
[HPO_4^{2-}] &= \Sigma Pi / P_{Pi}   \\
[HATP^{3-}] &= [ATP^{4-}] \frac{[H^+]}{K_{a}^{ATP}}  \\
[HADP^{2-}] &= [ADP^{3-}] \frac{[H^+]}{K_{a}^{ADP}}  \\
[H2PO_4^-] &= [HPO_4^{2-}] \frac{[H^+]}{K_{a}^{Pi}} \\
[MgATP^{2-}] &= [ATP^{4-}] \frac{[Mg^{2+}]}{K_{Mg}^{ATP}}  \\
[MgADP^-] &= [ADP^{3-}] \frac{[Mg^{2+}]}{K_{Mg}^{ADP}}  \\
\end{aligned}
$$

| Parameter       | Value | Unit | Description                                 |
| --------------- | ----- | ---- | ------------------------------------------- |
| $\delta_H$      | 1E-5  | -    | mitochondria [H^+] buffering capacity       |
| $pK_{a}^{ATP}$  | 6.48  | -    | pK of ATP acid dissociation constant        |
| $pK_{a}^{ADP}$  | 6.38  | -    | pK of ADP acid dissociation constant        |
| $pK_{a}^{Pi}$   | 6.75  | -    | pKa of phosphate acid dissociation constant |
| $pK_{Mg}^{ATP}$ | 4.19  | -    | pK of ATP magnesium dissociation constant   |
| $pK_{Mg}^{ADP}$ | 3.25  | -    | pK of ADP magnesium dissociation constant   |
| $pK_{a}^{SUC}$  | 5.2   | -    | pK of succinic acid dissociation constant   |
| $pKw$           | 14    |      | pK of water acid dissociation constant      |


[^Wei2011]: Wei AC, Aon MA, O'Rourke B, Winslow RL, Cortassa S. Mitochondrial energetics, pH regulation, and ion dynamics: a computational-experimental approach. Biophys J. 2011;100(12):2894-903. [PMC3123977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3123977/)
