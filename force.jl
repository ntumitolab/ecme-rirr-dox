#=
Force generation and troponin model by Rice et al. Use in Li et al (2015)
=#
using Parameters
include("common.jl")
# Force generation, troponin and crossbridge (XB) parameters
@with_kw struct XBParams
    TOT_LTRPN = 0.07  # Total pool of low affinity troponin ca binding sites (mM)
    TOT_HTRPN = 0.14  # Total pool of high affinity troponin ca binding sites (mM)
    K_P_LTRPN = 100.0  # Forward constant of low-affinity troponin sites
    K_M_LTRPN = 4e-2  # Backward constant of low-affinity troponin sites
    K_P_HTRPN = 100.  # Forward constant of high-affinity troponin sites
    K_M_HTRPN = 3.3e-4  # Backward constant of high-affinity troponin sites
    K_PN_TROP = 0.04  # Troponin rate constant (permissive -> nonpermissive)
    SL = 2.15  # Sarcomere length (μm)
    F_XB = 0.05
    G_MIN = 0.1
    G_OFF = 0.01
    V_MAX_AM = 7.2E-3
    KM_ATP_AM = 0.03
    KI_ADP_AM = 0.26
    ζ = 0.1  # Force constant (N/mm²)
    G01 = G_MIN
    G12 = 2 * G_MIN
    G23 = 3 * G_MIN
    F01 = 3 * F_XB
    F12 = 10 * F_XB
    F23 = 7 * F_XB
    ΣPATHS = G01 * G12 * G23 + F01 * G12 * G23 + F01 * F12 * G23 + F01 * F12 * F23
    P1_MAX = F01 * G12 * G23 / ΣPATHS
    P2_MAX = F01 * F12 * G23 / ΣPATHS
    P3_MAX = F01 * F12 * F23 / ΣPATHS
    Φ = 1 + (2.3 - SL) / ((2.3 - 1.7) ^ 1.6)
    G01_SL = Φ * G01
    G12_SL = Φ * G12
    G23_SL = Φ * G23
    G01_OFF = Φ * G_OFF
    N_TROP = 3.5 * SL - 2.0
    K_CA_TRPN = K_M_LTRPN / K_P_LTRPN  # Activation factor for calcium of low affinity troponin sites
    K½_TRPN = _mm(1.7e-3 - 0.8e-3 * (SL - 1.7) / 0.6, K_CA_TRPN)
    FORCE_NORM_DENOM = P1_MAX + P2_MAX + P3_MAX
    FORCE_DENOM = P1_MAX + 2 * P2_MAX + 3 * P3_MAX
    V_MAX_AM_SCALE = V_MAX_AM / (F01 + F12 + F23)
end

# Contractile force (N/mm^2)
force(p1, p2, p3, n1, ζ, FORCE_DENOM) = ζ / FORCE_DENOM * (p1 + n1 + 2 * p2 + 3 * p3)

# ATP hydrolysis rate by myofibrils (mM/ms)
function vam(p0, p1, p2, atp_i, adp_i, V_MAX_AM_SCALE, F01, F12, F23, KM_ATP_AM, KI_ADP_AM)
    vAm = V_MAX_AM_SCALE * (F01 * p0 + F12 * p1 + F23 * p2) * _mm(atp_i * _mm(KI_ADP_AM, adp_i), KM_ATP_AM)
end
function vam(p0, p1, p2, atp_i, adp_i, pXB::XBParams)
    @unpack V_MAX_AM_SCALE, F01, F12, F23, KM_ATP_AM, KI_ADP_AM = pXB
    vam(p0, p1, p2, atp_i, adp_i, V_MAX_AM_SCALE, F01, F12, F23, KM_ATP_AM, KI_ADP_AM)
end

function _force_denom(p1, p2, p3, n1, DENOM)
    return (p1 + p2 + p3 + n1) / DENOM
end

# Insider rates of crossbridge(xb) system
function _d_ltr_ca(ltr_ca, ca_i, p1, p2, p3, n1, K_P_LTRPN, K_M_LTRPN, TOT_LTRPN, FORCE_NORM_DENOM)
    K_P_LTRPN * ca_i * (TOT_LTRPN - ltr_ca) - K_M_LTRPN * ltr_ca * (1 - 2 / 3 * _force_denom(p1, p2, p3, n1, FORCE_NORM_DENOM))
end
function _d_htr_ca(htr_ca, ca_i, K_P_HTRPN, K_M_HTRPN, TOT_HTRPN)
    K_P_HTRPN * ca_i * (TOT_HTRPN - htr_ca) - K_M_HTRPN * htr_ca
end

# Calculates the rates of the troponin and crossbridge(xb) system (scalar version)
function xb_system(p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i, TOT_LTRPN, TOT_HTRPN,
                  K_PN_TROP, K½_TRPN, N_TROP, F01, F12, F23, G01_SL, G12_SL,
                  G23_SL, G01_OFF, FORCE_NORM_DENOM, K_P_LTRPN, K_P_HTRPN, K_M_HTRPN, K_M_LTRPN)
    k_np_trop = K_PN_TROP * (ltr_ca / (TOT_LTRPN * K½_TRPN)) ^ N_TROP
    n0 = one(p0) - (p0 + p1 + p2 + p3 + n1)
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
    d_p0 = v10 + v40 - v01 - v04
    d_p1 = v01 + v21 + v51 - v10 - v12 - v15
    d_p2 = v12 + v32 - v21 - v23
    d_p3 = v23 - v32
    d_n1 = v15 - v51 - v54
    d_ltr_ca = _d_ltr_ca(ltr_ca, ca_i, p1, p2, p3, n1, K_P_LTRPN, K_M_LTRPN, TOT_LTRPN, FORCE_NORM_DENOM)
    d_htr_ca = _d_htr_ca(htr_ca, ca_i, K_P_HTRPN, K_M_HTRPN, TOT_HTRPN)
    return d_p0, d_p1, d_p2, d_p3, d_n1, d_ltr_ca, d_htr_ca
end

function xb_system(p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i, pXB::XBParams)
    @unpack (TOT_LTRPN, TOT_HTRPN, K_PN_TROP, K½_TRPN, N_TROP, F01, F12, F23, G01_SL, G12_SL,
             G23_SL, G01_OFF, FORCE_NORM_DENOM, K_P_LTRPN, K_P_HTRPN, K_M_HTRPN, K_M_LTRPN) = pXB
    xb_system(p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i, TOT_LTRPN, TOT_HTRPN,
             K_PN_TROP, K½_TRPN, N_TROP, F01, F12, F23, G01_SL, G12_SL,
             G23_SL, G01_OFF, FORCE_NORM_DENOM, K_P_LTRPN, K_P_HTRPN, K_M_HTRPN, K_M_LTRPN)
end

# Calcium flux to troponin (mM/ms)
jtrpn(d_ltr_ca, d_htr_ca) = d_ltr_ca + d_htr_ca
