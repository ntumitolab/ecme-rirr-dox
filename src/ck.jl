using Parameters
include("common.jl")

# Creatine kinase and creatine shuttle parameters
@with_kw struct CKParams
    KCK_IN= 1.4E-4
    KCK_MT = 1.33E-6
    KTR_CR = 2E-3
    KEQ_CR = 0.0095
    # Basal consumption of cytosolic ATP (mM/ms)
    V_ATPASE_CYTO = 1E-5
    # Pool of creatine and creatine phosphate (mM)
    CR_T = 25.0
end

# CK shuttle rates
function ck_system(atp_i, atp_ic, crp_i, crp_ic, V_ATPASE_CYTO, KTR_CR, KCK_IN,
                  KCK_MT, KEQ, CR_T, AXP_I_T)
    cr_ic = CR_T - crp_ic
    cr_i = CR_T - crp_i
    adp_ic = AXP_I_T - atp_ic
    adp_i = AXP_I_T - atp_i

    vck_cyto = KCK_IN * (atp_ic * cr_ic - adp_ic * crp_ic / KEQ)
    vck_mito = KCK_MT * (atp_i * cr_i - adp_i * crp_i / KEQ)
    vtr_crp = KTR_CR * (crp_i - crp_ic)

    d_atp_in_ck = -vck_mito
    d_atp_ic = -vck_cyto - V_ATPASE_CYTO
    d_crp_i = vck_mito - vtr_crp
    d_crp_ic = vtr_crp + vck_cyto
    return (d_atp_in_ck, d_atp_ic, d_crp_i, d_crp_ic)
end

function ck_system(atp_i, atp_ic, crp_i, crp_ic, pCK::CKParams, AXP_I_T)
    @unpack V_ATPASE_CYTO, KTR_CR, KCK_IN, KCK_MT, KEQ_CR, CR_T = pCK
    ck_system(atp_i, atp_ic, crp_i, crp_ic, V_ATPASE_CYTO, KTR_CR, KCK_IN, KCK_MT, KEQ_CR, CR_T, AXP_I_T)
end
