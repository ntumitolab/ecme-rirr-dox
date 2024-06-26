"Creatine kinase (CK)"
function get_ck_sys(atp_i, adp_i; name=:cksys)
    @parameters KCK_IN = 0.14Hz [description = "Rate constant of creatine kinase (cytosolic)"]
    @parameters KCK_MT = 1.33E-3Hz [description = "Rate constant of creatine kinase (juxta-mitochondrial)"]
    @parameters KTR_CR = 2Hz [description = "Diffusion rate of creatine phosphate"]
    @parameters KEQ_CK = 0.0095 [description = "Eqilibrium constant of creatine kinase (PCr-forming)"]
    @parameters V_ATPASE_CYTO = 0.13E-5kHz [description = "Basal consumption of cytosolic ATP"]
    @parameters ΣCR = 25.0mM [description = "Pool of creatine and creatine phosphate"]
    @parameters ΣA_ic = 8mM
    @variables begin
        cr_ic(t)
        crp_ic(t) = 5.1291mM
        cr_i(t)
        crp_i(t) = 5.0297mM
        adp_ic(t) = 0.29mM
        atp_ic(t)
        vck_mito(t)
    end

    vck_cyto = KCK_IN * (atp_ic * cr_ic - adp_ic * crp_ic / KEQ_CK)
    vtr_crp = KTR_CR * (crp_i - crp_ic)
    v_aptase_cyto = V_ATPASE_CYTO * atp_ic

    eqs = [
        vck_mito ~ KCK_MT * (atp_i * cr_i - adp_i * crp_i / KEQ_CK),
        ΣCR ~ cr_ic + crp_ic,
        ΣCR ~ cr_i + crp_i,
        ΣA_ic ~ adp_ic + atp_ic,
        D(adp_ic) ~ vck_cyto + v_aptase_cyto,
        D(crp_i) ~ vck_mito - vtr_crp,
        D(crp_ic) ~ vtr_crp + vck_cyto
    ]
    return ODESystem(eqs, t; name)
end
