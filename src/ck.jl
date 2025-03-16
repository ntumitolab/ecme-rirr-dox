"Creatine kinase (CK)"
function get_ck_sys(atp_i, adp_i; name=:cksys)
    @parameters begin
        KCK_IN = 0.14Hz     # Rate constant of creatine kinase (cytosolic)
        KCK_MT = 1.33E-3Hz  # Rate constant of creatine kinase (juxta-mitochondrial)
        KTR_CR = 2Hz        # Diffusion rate of creatine phosphate
        KEQ_CK = 0.0095     # Eqilibrium constant of creatine kinase (PCr-forming)
        K_ATPASE_CYTO = 0.0013Hz # Basal consumption of cytosolic ATP
        ΣCR = 25.0mM        # Pool of creatine and creatine phosphate
        ΣA_ic = 8mM         # Pool of cytosolic ATP/ADP
    end

    @variables begin
        cr_ic(t)  # Conserved
        crp_ic(t) = 5.1291mM
        cr_i(t)   # Conserved
        crp_i(t) = 5.0297mM
        adp_ic(t) = 0.29mM
        atp_ic(t)
        vCK_mito(t)
        vCK_cyto(t)
        vTR_crp(t)
        vATPase_cyto(t)
    end
    eqs = [
        ΣCR ~ cr_ic + crp_ic,
        ΣCR ~ cr_i + crp_i,
        ΣA_ic ~ adp_ic + atp_ic,
        vATPase_cyto ~ K_ATPASE_CYTO * atp_ic,
        vCK_mito ~ KCK_MT * (atp_i * cr_i - adp_i * crp_i / KEQ_CK),
        vCK_cyto ~ KCK_IN * (atp_ic * cr_ic - adp_ic * crp_ic / KEQ_CK),
        vTR_crp ~ KTR_CR * (crp_i - crp_ic),
        D(adp_ic) ~ vCK_cyto + vATPase_cyto,
        D(crp_i) ~ vCK_mito - vTR_crp,
        D(crp_ic) ~ vTR_crp + vCK_cyto
    ]
    return ODESystem(eqs, t; name)
end
