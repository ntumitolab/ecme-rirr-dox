using Parameters
import .Utils: mM, kHz, μA, μM, cm², hill, hillr, iVT, pow_s

"Creatine kinase and creatine shuttle parameters"
@with_kw struct CK{R}
    KCK_IN::R = 1.4E-4kHz              # Rate constant of creatine kinase (cytosolic)
    KCK_MT::R = 1.33E-6kHz             # Rate constant of creatine kinase (juxta-mitochondrial)
    KTR_CR::R = 2E-3kHz                # Diffusion rate of creatine phosphate
    KEQ::R = 0.0095                    # Eqilibrium constant of creatine kinase (PCr-forming)
    V_ATPASE_CYTO::R = 1E-5mM * kHz    # Basal consumption of cytosolic ATP
    ΣCR::R = 25.0mM                    # Pool of creatine and creatine phosphate
    ΣATPc::R = 8.0mM                   # Pool of ATP and ADP
end

"CK shuttle rates"
function ck_system(atp_i, adp_i, atp_ic, crp_i, crp_ic, p::CK)
    @unpack V_ATPASE_CYTO, KTR_CR, KCK_IN, KCK_MT, KEQ, ΣCR, ΣATP = p
    cr_ic = ΣCR - crp_ic
    cr_i = ΣCR - crp_i
    adp_ic = ΣATP - atp_ic
    vck_cyto = KCK_IN * (atp_ic * cr_ic - adp_ic * crp_ic / KEQ)
    vck_mito = KCK_MT * (atp_i * cr_i - adp_i * crp_i / KEQ)
    vtr_crp = KTR_CR * (crp_i - crp_ic)
    v_aptase_cyto = V_ATPASE_CYTO * hill(atp_ic, 10μM)
    return (;
        d_atp_in_ck = -vck_mito,
        d_atp_ic = -vck_cyto - v_aptase_cyto,
        d_crp_i = vck_mito - vtr_crp,
        d_crp_ic = vtr_crp + vck_cyto)
end

"Na-K ATPase (NKA) parameters"
@with_kw struct NKA{R}
    k_o::R = 5.4mM
    na_o::R = 140.0mM
    IMAX0::R = 3.147μA / cm²   # Max Na-K ATPase current
    KM_NA::R = 10.0mM          # Na half-saturate constant of Na-K ATPase
    KM_K::R = 1.5mM            # K half-saturate constant of Na-K ATPase
    KM_ATP::R = 8μM            # ATP half-saturate constant of Na-K ATPase
    KI_ADP::R = 100μM          # ADP half-saturate constant of Na-K ATPase
    F_NA::R = expm1(na_o / 67.3mM) / 7  # Factor of extracellular sodium
    IMAX::R = hill(k_o, KM_K) * IMAX0
end

"Na-K ATPase (NKA) current"
function i_nak(na_i, atp_i, adp_i, vm, p::NKA, evfrt = exp(vm * iVT))
    @unpack KM_ATP, KI_ADP, KM_NA, IMAX, F_NA = p
    f_na = hill(na_i, KM_NA, 1.5)
    f_nak_inv = 1.0 + 0.1245 * exp(-0.1iVT * vm) + 0.0365 * F_NA / evfrt
    f_atp = hill(atp_i * hillr(adp_i, KI_ADP), KM_ATP)
    return IMAX * f_na * f_atp / f_nak_inv
end

(p::NKA)(na_i, atp_i, adp_i, vm, evfrt = exp(vm * iVT)) = i_nak(na_i, atp_i, adp_i, vm, p, evfrt)

"Plasma membrane calcium ATPase (PMCA) parameters"
@with_kw struct PMCA{R}
    IMAX::R = 0.575μA / cm²  # Max PMCA current
    KM_CA::R = 0.5μM         # Ca half-saturation constant
    KM1_ATP::R = 12μM        # ATP 1st half-saturation constant
    KM2_ATP::R = 230μM       # ATP 2nd half-saturation constant (mM)
    KI_ADP::R = 1.0mM        # ADP half-inhibition constant (mM)
end

"Plasma membrane calcium pump current"
function i_pca(ca_i, atp_i, adp_i, p::PMCA)
    @unpack KM1_ATP, KI_ADP, KM2_ATP, KM_CA, IMAX = p
    f_atp = hill(atp_i * hillr(adp_i, KI_ADP), KM1_ATP) + hill(atp_i, KM2_ATP)
    f_ca = hill(ca_i, KM_CA)
    return IMAX * f_atp * f_ca
end

(p::PMCA)(ca_i, atp_i, adp_i) = i_pca(ca_i, atp_i, adp_i, p)

"SER Ca ATPase (SERCA) parameters"
@with_kw struct SERCA{R}
    VF::R = 2.989E-4mM * kHz   # Max forward rate
    VR::R = 3.179E-4mM * kHz   # Max reverse rate
    KMF_CA::R = 0.24μM         # Michaelis constant for Ca of forward reaction
    KMR_CA::R = 1.64269mM      # Michaelis constant for Ca of reverse reaction
    NFB::R = 1.4               # Cooperativity of forward reaction
    NRB::R = 1.0               # Cooperativity of reverse reaction
    KM_ATP::R = 10μM           # Michaelis constant for ATP (mM)
    KI1_ADP::R = 140μM         # 1st michaelis constant for ADP (mM)
    KI2_ADP::R = 5.1mM         # 2nd michaelis constant for ADP (mM)
end

"SERCA calcium flux (jUp)"
function j_up(ca_i, ca_nsr, atp_i, adp_i, p::SERCA)
    @unpack NFB, NRB, KMF_CA, KMR_CA, KM_ATP, KI1_ADP, KI2_ADP, VF, VR = p
    fb = pow_s(ca_i / KMF_CA, NFB)
    rb = pow_s(ca_nsr / KMR_CA, NRB)
    f_atp_inv = KM_ATP / atp_i * p_one(adp_i / KI1_ADP) + p_one(adp_i / KI2_ADP)
    return f_atp * (VF * fb - VR * rb) / (p_one(fb + rb) * f_atp_inv)
end

(p::SERCA)(ca_i, ca_nsr, atp_i, adp_i) = j_up(ca_i, ca_nsr, atp_i, adp_i, p)