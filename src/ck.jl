function get_ck_eqs(; atp_i=7.9mM, adp_i=0.1mM)
    @parameters begin
        KCK_IN = 0.14Hz     # Rate constant of creatine kinase (cytosolic)
        KCK_MT = 1.33E-3Hz  # Rate constant of creatine kinase (juxta-mitochondrial)
        KTR_CR = 2Hz        # Diffusion rate of creatine phosphate
        iKEQ_CK = inv(0.0095) # Eqilibrium constant of creatine kinase (Cr-forming)
        KATPASE_CYTO = 1.3E-3Hz # Basal consumption rate of cytosolic ATP
        ΣCR = 25mM          # Pool of creatine and creatine phosphate
        ΣA_ic = 8mM         # Pool of cytosolic ATP + ADP
    end
    @variables begin
        cr_ic(t)  ## Conserved
        crp_ic(t) = 5.1291mM
        cr_i(t)   ## Conserved
        crp_i(t) = 5.0297mM
        adp_ic(t) = 0.29mM
        atp_ic(t)  ## Conserved
        vCK_mito(t)
        vCK_cyto(t)
        vTR_crp(t)
    end

    eqs_ck = [
        ΣCR ~ cr_ic + crp_ic,
        ΣCR ~ cr_i + crp_i,
        ΣA_ic ~ adp_ic + atp_ic,
        vCK_mito ~ KCK_MT * (atp_i * cr_i - adp_i * crp_i * iKEQ_CK),
        vCK_cyto ~ KCK_IN * (atp_ic * cr_ic - adp_ic * crp_ic * iKEQ_CK),
        vTR_crp ~ KTR_CR * (crp_i - crp_ic),
        D(adp_ic) ~ vCK_cyto + KATPASE_CYTO * atp_ic,
        D(crp_i) ~ vCK_mito - vTR_crp,
        D(crp_ic) ~ vTR_crp + vCK_cyto
    ]

    return (; eqs_ck, vCK_mito)
end

"Creatine kinase (CK)"
function get_ck_sys(; atp_i, adp_i, name=:cksys)
    eqs_ck, _ = get_ck_eqs(; atp_i, adp_i)
    return System(eqs_ck, t; name)
end
