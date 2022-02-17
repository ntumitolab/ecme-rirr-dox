using Parameters
import .Utils: mM, μm, kHz, mm, mmr, μM

"Myofibril force generation, troponin and crossbridge (XB) parameters"
@with_kw struct MyoFibril
    ΣLTRPN = 0.07mM             # Total pool of low affinity troponin ca binding sites
    ΣHTRPN = 0.14mM             # Total pool of high affinity troponin ca binding sites
    K_P_LTRPN = 100.0kHz / mM   # Forward constant of low-affinity troponin sites
    K_M_LTRPN = 4e-2kHz         # Backward constant of low-affinity troponin sites
    K_P_HTRPN = 100.0kHz / mM   # Forward constant of high-affinity troponin sites
    K_M_HTRPN = 3.3e-4kHz       # Backward constant of high-affinity troponin sites
    K_PN_TROP = 0.04kHz         # Troponin rate constant (permissive -> nonpermissive)
    SL = 2.15                   # Sarcomere length (μm)
    F_XB = 0.05kHz              # Transition rate from weak to strong cross bridge
    G_MIN = 0.1kHz              # Minimum transition rate from strong to weak cross bridge
    G_OFF = 0.01kHz
    V_MAX_AM = 7.2E-3mM / kHz   # Maximal rate of ATP hydrolysis by myofibrils (AM ATPase) 
    KM_ATP_AM = 0.03mM          # ATP half saturation constant of AM ATPase
    KI_ADP_AM = 0.26mM          # ADP inhibition constant of AM ATPase
    ζ = 0.1                     # Conversion factor normalizing to  physiological force (N/mm²)
    G01 = 1 * G_MIN
    G12 = 2 * G_MIN
    G23 = 3 * G_MIN
    F01 = 3 * F_XB
    F12 = 10 * F_XB
    F23 = 7 * F_XB
    ΣPATHS = G01 * G12 * G23 + F01 * G12 * G23 + F01 * F12 * G23 + F01 * F12 * F23
    P1_MAX = F01 * G12 * G23 / ΣPATHS
    P2_MAX = F01 * F12 * G23 / ΣPATHS
    P3_MAX = F01 * F12 * F23 / ΣPATHS
    Φ = 1 + (2.3 - SL) / ((2.3 - 1.7)^1.6)
    G01_SL = Φ * G01
    G12_SL = Φ * G12
    G23_SL = Φ * G23
    G01_OFF = Φ * G_OFF
    N_TROP = 3.5 * SL - 2.0
    K_CA_TRPN = K_M_LTRPN / K_P_LTRPN  # Activation factor for calcium of low affinity troponin sites
    K½_TRPN = mm(1.7μM - 0.8μM * (SL - 1.7) / 0.6, K_CA_TRPN)
end

# Contractile force (N/mm^2)
force(p1, p2, p3, n1, p::MyoFibril) = p.ζ * (p1 + n1 + 2 * p2 + 3 * p3) / (p.P1_MAX + 2 * p.P2_MAX + 3 * p.P3_MAX)
# Normalized force
force_norm(p1, p2, p3, n1, p::MyoFibril) = (p1 + p2 + p3 + n1) / (p.P1_MAX + p.P2_MAX + p.P3_MAX)

# ATP hydrolysis rate by myofibrils
function v_am(p0, p1, p2, atp_i, adp_i, p::MyoFibril)
    @unpack V_MAX_AM, F01, F12, F23, KM_ATP_AM, KI_ADP_AM = p
    f_atp = mm(atp_i * mmr(adp_i, KI_ADP_AM), KM_ATP_AM)
    return V_MAX_AM / (F01 + F12 + F23) * (F01 * p0 + F12 * p1 + F23 * p2) * f_atp
end

# Calculates the rates of the troponin and crossbridge (xb) system 
function xb_system(p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i, p::MyoFibril)
    @unpack K_PN_TROP, ΣLTRPN, ΣHTRPN, K½_TRPN, N_TROP = p
    # Cross bridge
    k_np_trop = K_PN_TROP * (ltr_ca / (ΣLTRPN * K½_TRPN))^N_TROP
    n0 = one_m(p0 + p1 + p2 + p3 + n1)
    v01 = F01 * p0
    v10 = G01_SL * p1
    v12 = F12 * p1
    v21 = G12_SL * p2
    v23 = F23 * p2
    v32 = G23_SL * p3
    v04 = K_PN_TROP * p0
    v40 = k_np_trop * n0
    v15 = K_PN_TROP * p1
    v51 = k_np_trop * n1
    v54 = G01_OFF * n1
    # Troponin system
    d_ltr_ca = K_P_LTRPN * ca_i * (ΣLTRPN - ltr_ca) - K_M_LTRPN * ltr_ca * (1 - 2 / 3 * force_norm(p1, p2, p3, n1, p))
    d_htr_ca = K_P_HTRPN * ca_i * (ΣHTRPN - htr_ca) - K_M_HTRPN * htr_ca
    return (;
        d_p0 = v10 + v40 - v01 - v04,
        d_p1 = v01 + v21 + v51 - v10 - v12 - v15,
        d_p2 = v12 + v32 - v21 - v23,
        d_p3 = v23 - v32,
        d_n1 = v15 - v51 - v54,
        d_ltr_ca = d_ltr_ca,
        d_htr_ca = d_htr_ca,
        j_trpn = d_ltr_ca + d_htr_ca)
end
