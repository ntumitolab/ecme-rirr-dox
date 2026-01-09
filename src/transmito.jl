# Mitochondrial transport systems
get_default_NHE_params() = ComponentArray(
    K1_P_NHE = 0.0252 / ms,    # NHE forward rate constant
    K1_M_NHE = 0.0429 / ms,    # NHE backward rate constant
    K4_P_NHE = 0.16 / ms,      # NHE forward rate constant
    K4_M_NHE = 0.0939 / ms,    # NHE backward rate constant
    KM_NA_NHE = 24.0mM,      # Na dissociation constant
    KM_H_NHE = 1.584E-4mM,   # H+ dissociation constant
    H_KI_NHE = 3.02E-6mM,
    NI_NHE = 3,              # Hill coefficient for proton binding
    ET_NHE = 0.00785mM       # NHE concentration
)

function vNHE(; na_m, na_i, h_m, h_i, p=get_default_NHE_params())
    @unpack K1_P_NHE, K1_M_NHE, K4_P_NHE, K4_M_NHE, KM_NA_NHE, KM_H_NHE, H_KI_NHE, NI_NHE, ET_NHE= p
    # NHE transition rates
    β₁⁺ = K1_P_NHE * hil(na_m * hil(KM_H_NHE, h_m), KM_NA_NHE)
    β₁⁻ = K1_M_NHE * hil(na_i * hil(KM_H_NHE, h_i), KM_NA_NHE)
    β₂⁺ = K4_P_NHE * hil(h_i * hil(KM_NA_NHE, na_i), KM_H_NHE)
    β₂⁻ = K4_M_NHE * hil(h_m * hil(KM_NA_NHE, na_m), KM_H_NHE)
    return ET_NHE / (1 + (H_KI_NHE / h_i)^NI_NHE) * (β₁⁺ * β₂⁺ - β₁⁻ * β₂⁻) / (β₁⁺ + β₂⁺ + β₁⁻ + β₂⁻)
end

get_default_PiC_params() = ComponentArray(
    ET_PIC = 4.9mM,              # PIC concentration
    KM_PI_IN_PIC = 11.06mM,      # Cytosolic Pi binding constant
    KM_PI_MT_PIC = 11.06mM,      # Mitochondrial Pi binding constant
    KM_OH_IN_PIC = 4.08E-5mM,    # Cytosolic OH binding constant
    KM_OH_MT_PIC = 4.08E-5mM,    # Mitochondrial OH binding constant
    KF_PIC = 90 / minute,        # Max forward rate of PIC
    KB_PIC = 90 / minute,        # Max backward rate of PIC
    KWATER = 1E-14 * Molar^2     # Water ion product
)

function vPiC(; pi_m, pi_i, h_m, h_i, p=get_default_PiC_params())
    @unpack ET_PIC, KM_PI_IN_PIC, KM_PI_MT_PIC, KM_OH_IN_PIC, KM_OH_MT_PIC, KF_PIC, KB_PIC, KWATER = p
    # PiC transition rates
    dhpi⁻_i = pi_i * hil(h_i, KA_PI)
    dhpi⁻_m = pi_m * hil(h_m, KA_PI)
    f_pi = dhpi⁻_i / KM_PI_IN_PIC
    f_pm = dhpi⁻_m / KM_PI_MT_PIC
    f_ohi = KWATER / h_i / KM_OH_IN_PIC
    f_ohm = KWATER / h_m / KM_OH_MT_PIC
    ϕf = f_pi * f_ohm
    ϕb = f_pm * f_ohi
    return ET_PIC * (KF_PIC * ϕf - KB_PIC * ϕb) / (1 + f_pi + f_pm + f_ohi + f_ohm + ϕf + ϕb)
end

get_default_nclx_params() = ComponentArray(
    VMAX_NCLX = 0.1μM / ms,        # Mitochondrial NCE (NCLX)
    B_NCLX = 0.5,
    KM_NA_NCLX = 9.4mM,
    KM_CA_NCLX = 0.375μM,
    DPSI_OFFSET_MCU = 91mV,
)

get_default_mcu_params() = ComponentArray(
    VMAX_MCU = 27.5μM / ms,       # Mitochondrial Ca uniporter (MCU)
    DPSI_OFFSET_MCU = 91mV,
    KACT_MCU = 0.38μM,
    iKTRANS_MCU = inv(19μM),
    L_MCU = 110,
    N_MCU = 2.8
)

function vNCLX(; na_i, ca_m, ca_i, dpsi, p=get_default_nclx_params())
    @unpack VMAX_NCLX, B_NCLX, KM_NA_NCLX, KM_CA_NCLX, DPSI_OFFSET_MCU = p
    f_na = hil(na_i, KM_NA_NCLX)^3
    f_ca = hil(ca_m, KM_CA_NCLX)
    r_ca = ca_m / ca_i
    return VMAX_NCLX * exp(iVT * B_NCLX * (dpsi - DPSI_OFFSET_MCU)) * r_ca * f_na * f_ca
end

function vMCU(; ca_i, dpsi, p=get_default_mcu_params())
    @unpack VMAX_MCU, DPSI_OFFSET_MCU, KACT_MCU, iKTRANS_MCU, L_MCU, N_MCU = p
    v2frt = 2 * iVT * (dpsi - DPSI_OFFSET_MCU)
    f_trans = ca_i * iKTRANS_MCU
    f_act = L_MCU * (NaNMath.pow(hil(KACT_MCU, ca_i), N_MCU))
    return VMAX_MCU * f_trans * (f_trans + 1)^3 * exprel(-v2frt) / ((f_trans + 1)^4 + f_act)
end
