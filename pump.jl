using Parameters
include("common.jl")

# Na-K ATPase (NKA) parameters
@with_kw struct NKAParams
    k_o = 5.4
    na_o = 140.0
    # Max Na-K ATPase current (μA/cm²)
    IMAX = 3.147
    # Na half-saturate constant of Na-K ATPase (mM)
    KM_NA = 10.0
    # K half-saturate constant of Na-K ATPase (mM)
    KM_K = 1.5
    # ATP half-saturate constant of Na-K ATPase (mM)
    KM_ATP = 8E-3
    # ADP half-saturate constant of Na-K ATPase (mM)
    KI_ADP = 0.1
    # Factor of extracellular potassium
    FAC_K = _mm(k_o, KM_K)
    # Factor of extracellular sodium
    FAC_NA = expm1(na_o / 67.3) / 7
end

# Na/K ATPase current (iNaK) (uA/cm^2), scalar version
function inak(na_i, atp_i, adp_i, vm, evfrtm1, KM_ATP, KI_ADP, KM_NA, IMAX, FAC_K, FAC_NA)
    evfrt = evfrtm1 + 1
    f_na = _hills(na_i, KM_NA, 1.5)
    f_nak_inv = 1.0 + 0.1245 * _ra(-0.1 * vm) + 0.0365 * FAC_NA / evfrt
    f_atp = _mm(atp_i * _mm(KI_ADP, adp_i), KM_ATP)
    iNak = IMAX * FAC_K * f_na * f_atp / f_nak_inv
end

function inak(na_i, atp_i, adp_i, vm, evfrtm1, pNKA::NKAParams)
    @unpack KM_ATP, KI_ADP, KM_NA, IMAX, FAC_K, FAC_NA = pNKA
    inak(na_i, atp_i, adp_i, vm, evfrtm1, KM_ATP, KI_ADP, KM_NA, IMAX, FAC_K, FAC_NA)
end

# Plasma membrane calcium ATPase (PMCA) parameters
@with_kw struct PMCAParams
    # Max PMCA current (μA/cm²)
    IMAX = 0.575
    # Ca half-saturation constant (mM)
    KM_CA = 5E-4
    # ATP 1st half-saturation constant (mM)
    KM1_ATP = 0.012
    # ATP 2nd half-saturation constant (mM)
    KM2_ATP = 0.23
    # ADP half-inhibition constant (mM)
    KI_ADP = 1.0
end

# PMCA current (iPCa) (uA/cm^2), scalar version
function ipca(ca_i, atp_i, adp_i, KM1_ATP, KI_ADP, KM2_ATP, KM_CA, IMAX)
    f_atp = _mm(atp_i * _mm(KI_ADP, adp_i), KM1_ATP) + _mm(atp_i, KM2_ATP)
    f_ca = _mm(ca_i, KM_CA)
    iPca = IMAX * f_atp * f_ca
end

function ipca(ca_i, atp_i, adp_i, pPMCA::PMCAParams)
    @unpack KM1_ATP, KI_ADP, KM2_ATP, KM_CA, IMAX = pPMCA
    ipca(ca_i, atp_i, adp_i, KM1_ATP, KI_ADP, KM2_ATP, KM_CA, IMAX)
end

# SER ca ATPase (SERCA) parameters
# Params from Li et al. (2015)
@with_kw struct SERCAParams
    # Max forward rate (mM/ms)
    VF = 2.989E-4
    # Max reverse rate (mM/ms)
    VR = 3.179E-4
    # Michaelis constant for Ca of forward reaction (mM)
    KMF_CA = 1.5E-4
    # Michaelis constant for Ca of reverse reaction (mM)
    KMR_CA = 3.3
    # Cooperativity of forward reaction
    NFB = 0.55
    # Cooperativity of reverse reaction
    NRB = 0.5
    # Michaelis constant for ATP (mM)
    KM_ATP = 0.01
    # 1st michaelis constant for ADP (mM)
    KI1_ADP = 0.14
    # 2nd michaelis constant for ADP (mM)
    KI2_ADP = 5.1
end

# SERCA calcium flux (jUp) (mM/ms), scalar version
function jup(ca_i, ca_nsr, atp_i, adp_i,
             NFB, NRB, KMF_CA, KMR_CA, KM_ATP, KI1_ADP, KI2_ADP, VF, VR)
    fb = (ca_i / KMF_CA) ^ NFB
    rb = (ca_nsr / KMR_CA) ^ NRB
    f_atp = inv(KM_ATP / atp_i * _mm_reci(KI1_ADP, adp_i) + _mm_reci(KI2_ADP, adp_i))
    jUp = f_atp * (VF * fb - VR * rb) / (1.0 + fb + rb)
end

function jup(ca_i, ca_nsr, atp_i, adp_i, pSERCA::SERCAParams)
    @unpack NFB, NRB, KMF_CA, KMR_CA, KM_ATP, KI1_ADP, KI2_ADP, VF, VR = pSERCA
    jup(ca_i, ca_nsr, atp_i, adp_i,
        NFB, NRB, KMF_CA, KMR_CA, KM_ATP, KI1_ADP, KI2_ADP, VF, VR)
end
