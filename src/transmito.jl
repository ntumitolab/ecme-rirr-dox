function get_mitoh_sys(na_m, na_i, h_m, h_i, pi_m, pi_i; name=:mitohsys)
    @parameters begin
        ### Na-H exchanger (NHE)
        K1_P_NHE = 0.0252 / ms    # NHE forward rate constant
        K1_M_NHE = 0.0429 / ms    # NHE backward rate constant
        K4_P_NHE = 0.16 / ms      # NHE forward rate constant
        K4_M_NHE = 0.0939 / ms    # NHE backward rate constant
        KM_NA_NHE = 24.0mM      # Na dissociation constant
        KM_H_NHE = 1.584E-4mM   # H+ dissociation constant
        H_KI_NHE = 3.02E-6mM
        NI_NHE = 3              # Hill coefficient for proton binding
        ET_NHE = 0.00785mM      # NHE concentration
        ### PiC (phosphate carrier)
        ET_PIC = 4.9mM              # PIC concentration (mM)
        KM_PI_IN_PIC = 11.06mM      # Cytosolic Pi binding constant
        KM_PI_MT_PIC = 11.06mM      # Mitochondrial Pi binding constant
        KM_OH_IN_PIC = 4.08E-5mM    # Cytosolic OH binding constant
        KM_OH_MT_PIC = 4.08E-5mM    # Mitochondrial OH binding constant
        KF_PIC = 90 / minute          # Max forward rate of PIC (mM/ms)
        KB_PIC = 90 / minute          # Max backward rate of PIC (mM/ms)
    end
    @variables vNHE(t) vPIC(t)
    # NHE transition rates
    β₁⁺ = K1_P_NHE * hil(na_m * hil(KM_H_NHE, h_m), KM_NA_NHE)
    β₁⁻ = K1_M_NHE * hil(na_i * hil(KM_H_NHE, h_i), KM_NA_NHE)
    β₂⁺ = K4_P_NHE * hil(h_i * hil(KM_NA_NHE, na_i), KM_H_NHE)
    β₂⁻ = K4_M_NHE * hil(h_m * hil(KM_NA_NHE, na_m), KM_H_NHE)

    # PiC transition rates
    dhpi⁻_i = pi_i * hil(h_i, KA_PI)
    dhpi⁻_m = pi_m * hil(h_m, KA_PI)
    f_pi = dhpi⁻_i / KM_PI_IN_PIC
    f_pm = dhpi⁻_m / KM_PI_MT_PIC
    f_ohi = KWATER / h_i / KM_OH_IN_PIC
    f_ohm = KWATER / h_m / KM_OH_MT_PIC
    ϕf = f_pi * f_ohm
    ϕb = f_pm * f_ohi

    eqs = [
        vNHE ~ ET_NHE / (1 + (H_KI_NHE / h_i)^NI_NHE) * (β₁⁺ * β₂⁺ - β₁⁻ * β₂⁻) / (β₁⁺ + β₂⁺ + β₁⁻ + β₂⁻),
        vPIC ~ ET_PIC * (KF_PIC * ϕf - KB_PIC * ϕb) / (1 + f_pi + f_pm + f_ohi + f_ohm + ϕf + ϕb)
    ]
    return System(eqs, t; name)
end

"Mitochondrial calcium handling"
function get_mitoca_sys(; na_i, ca_m, ca_i, dpsi, name=:mitocasys)
    @parameters begin
        # Mitochondrial NCE (NCLX)
        VMAX_NCLX = 0.1μM / ms
        B_NCLX = 0.5
        KM_NA_NCLX = 9.4mM
        KM_CA_NCLX = 0.375μM
        # Mitochondrial Ca uniporter (MCU)
        VMAX_MCU = 27.5μM / ms
        DPSI_OFFSET_MCU = 91mV
        KACT_MCU = 0.38μM
        iKTRANS_MCU = inv(19μM)
        L_MCU = 110
        N_MCU = 2.8
    end
    @variables vUni(t) vNaCa(t)
    v2frt = 2 * iVT * (dpsi - DPSI_OFFSET_MCU)
    f_trans = ca_i * iKTRANS_MCU
    f_act = L_MCU * (NaNMath.pow(hil(KACT_MCU, ca_i), N_MCU))
    v_mcu = VMAX_MCU * f_trans * (f_trans + 1)^3 * exprel(-v2frt) / ((f_trans + 1)^4 + f_act)
    f_na = hil(na_i, KM_NA_NCLX)^3
    f_ca = hil(ca_m, KM_CA_NCLX)
    r_ca = ca_m / ca_i
    v_naca = VMAX_NCLX * exp(iVT * B_NCLX * (dpsi - DPSI_OFFSET_MCU)) * r_ca * f_na * f_ca
    eqs = [vUni ~ v_mcu, vNaCa ~ v_naca]
    return System(eqs, t; name)
end
