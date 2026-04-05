# Citrate synthase
_vcs(kcat, et, oaa, ACCOA, km_oaa, km_accoa) = kcat * et * oaa * ACCOA / ((km_oaa + oaa) * (km_accoa + ACCOA))

# Aconitase
_vaco(kf, rkeq, citrate, isocitrate) = kf * (citrate - isocitrate * rkeq)

# IDH3 (Isocitrate dehydrogenase, NADH-producing)
function _vidh3(kcat, et, isocitrate, km_isoc, adp, km_adp, ca, km_ca, nad, km_nad, nadh, ki_nadh)
    vmax = kcat * et
    A = (isocitrate / km_isoc)^2 * (1 + adp / km_adp) * (1 + ca / km_ca)
    B = nad * ki_nadh / (km_nad * (ki_nadh + nadh))
    H = 1 + h_m / KH1_IDH + KH2_IDH / h_m
    return vmax * A * B / (H * A * B + A + B + 1)
end

# KGDH (alpha-ketoglutarate dehydrogenase)
function _vkgdh(kcat, et, akg, km_akg, nad_m, km_nad, mg_m, km_mg, ca_m, km_ca, h_m, kh1, kh2, ni)
    vmax = et * kcat
    f_h = 1 + h_m / kh1 + kh2 / h_m
    f_akg = NaNMath.pow(akg / km_akg, ni)
    f_mgca = (1 + mg_m / km_mg) * (1 + ca_m / km_ca)
    f_nad = nad_m / km_nad
    return vmax * f_akg * f_nad * f_mgca / (f_h * f_mgca * f_akg * f_nad + f_akg + f_nad)
end

# SL (Succinyl-coA lyase)
_vsl(kf, rkeq, scoa, adp_m, pi_m, succinate, atp_m, COA) = kf * (scoa * adp_m * pi_m - succinate * atp_m * COA * rkeq)

# SL (Succinyl-coA lyase) with binding polynomials
function _vsl_poly(kf, rkeq, scoa, adp_m, pi_m, succinate, atp_m, COA, h_m, mg_m)
    atp4, hatp, mgatp, atp_poly = breakdown_atp(atp_m, h_m, mg_m)
    adp3, hadp, mgadp, adp_poly = breakdown_adp(adp_m, h_m, mg_m)
    pi_poly = pipoly(h_m)
    suc_poly = sucpoly(h_m)
    rkeq = rkeq * (adp_poly * pi_poly) / (atp_poly * suc_poly)
    return kf * (scoa * adp_m * pi_m - succinate * (atp4 + hatp) * COA * rkeq)
end

# Fumarate hydratase
_vfh(kf, req, fumarate, malate) = kf * (fumarate - malate * req)

# Malate dehydrogenase (reversible)
function _vmdh_reverisble(vf, vb, oaa, km_oaa, malate, km_mal, nad_m, km_nad, nadh_m, km_nadh, h_m, koffset, kh1, kh2, kh3, kh4)
    f_ha = koffset + (kh1 * kh2 / (kh1 * kh2 + kh2 * h_m + h_m^2))
    f_hi = (h_m^2 / (h_m^2 + h_m * kh3 + kh3 * kh4))^2
    f_h = f_ha * f_hi
    f_nad = nad_m / km_nad
    f_mal = malate / km_mal
    f_oaa = oaa / km_oaa
    f_nadh = nadh_m / km_nadh
    return f_h * (vf * f_nad * f_mal - vb * f_oaa * f_nadh) / (1 + f_nad + f_nad * f_mal + f_oaa * f_nadh + f_nadh)
end

_vast_reverisble(kf, oaa, glu, akg, asp, keq) = kf * (oaa * glu - akg * asp / keq)

# TCA cycle model
function get_tca_eqs(; atp_m, adp_m, nad_m, nadh_m, ca_m, h_m=exp10(-7.6) * Molar, pi_m=8mM, mg_m=0.4mM, use_mg=false)
    @parameters begin
        ## Total TCA metabolite pool
        TCA_T = 1.3mM
        ## Citrate synthase
        KCAT_CS = 50Hz
        ET_CS = 400μM
        KM_ACCOA_CS = 12.6μM
        KM_OAA_CS = 0.64μM
        ACCOA = 1000μM
        ## ACO (aconitase)
        KF_ACO = 12.5Hz # 12.5Hz (Zhou, 2009)
        rKEQ_ACO = inv(2.22)
        ## IDH3 (Isocitrate dehydrogenase, NADH-producing)
        KCAT_IDH = 43Hz ## 50Hz
        ET_IDH = 109μM
        KI_NADH_IDH = 190μM
        KH1_IDH = 10nM
        KH2_IDH = 900nM
        KM_ISOC_IDH = 1520μM
        NI_ISOC_IDH = 2
        KM_NAD_IDH = 923μM
        KM_ADP_IDH = 620μM
        KM_CA_IDH = 0.5μM
        ## KGDH (alpha-ketoglutarate dehydrogenase)
        ET_KGDH = 500μM
        KCAT_KGDH = 50Hz
        KM_AKG_KGDH = 1940μM
        KM_NAD_KGDH = 38.7mM
        KH1_KGDH = 40nM
        KH2_KGDH = 70nM
        KM_MG_KGDH = 30.8μM
        KM_CA_KGDH = 0.15μM
        NI_AKG_KGDH = 1.2
        ## SL (Succinyl-coA lyase)
        COA = 20μM
        KF_SL = 28Hz / mM^2  ## Gauthier, 2013 vs 0.644 in Zhou, 2009
        rKEQ_SL = inv(3.115)
        ## FH (Fumarate hydrase) parameters
        KF_FH = 8.3Hz
        rKEQ_FH = 1.0
        ## MDH (Malate dehydrogenase), reversible
        KCAT_MDH = 126Hz
        ET_MDH = 154μM
        VF_MDH = KCAT_MDH * ET_MDH
        KEQ_MDH = 3.08e-4
        KH1_MDH = 11.31nM
        KH2_MDH = 26.7mM
        KH3_MDH = 6.68E-3nM
        KH4_MDH = 5.62nM
        K_OFFSET_MDH = 0.0399
        KM_NAD_MDH = 110μM
        KM_MAL_MDH = 450μM
        KM_OAA_MDH = 7μM
        KM_NADH_MDH = 17μM
        VR_MDH = VF_MDH * (KM_OAA_MDH * KM_NADH_MDH) / (KEQ_MDH * KM_NAD_MDH * KM_MAL_MDH)
        ### AAT (alanine aminotransferase, aka AST)
        KF_AAT = 21.4Hz / mM
        KEQ_AAT = 6.6
        GLU = 10mM   ## Glutamate
        ASP = 100μM  ## Aspartate
    end

    sts = @variables begin
        oaa(t) = 11.6μM     ## oxaloacetate
        isocitrate(t) = 51.6μM    ## isocitrate
        akg(t) = 51μM       ## alpha-ketoglutarate
        scoa(t) = 35μM      ## succinyl-CoA
        succinate(t) = 1.9μM      ## succinate
        fumarate(t) = 175μM      ## fumarate
        malate(t) = 160μM      ## malate
    end

    @variables begin
        vCS(t)
        vACO(t)
        vIDH(t)
        vKGDH(t)
        vSL(t)
        vFH(t)
        vMDH(t)
        vAAT(t)
        vSDH(t)
        citrate(t) ## Conserved = TCA_T - isocitrate - oaa - akg - scoa - succinate - fumarate - malate
    end

    if use_mg
        v_sl = _vsl_poly(KF_SL, rKEQ_SL, scoa, adp_m, pi_m, succinate, atp_m, COA, h_m, mg_m)
    else
        v_sl = _vsl(KF_SL, rKEQ_SL, scoa, adp_m, pi_m, succinate, atp_m, COA)
    end

    eqs_tca = [
        TCA_T ~ citrate + isocitrate + oaa + akg + scoa + succinate + fumarate + malate,
        vCS ~ _vcs(KCAT_CS, ET_CS, oaa, ACCOA, KM_OAA_CS, KM_ACCOA_CS),
        vACO ~ _vaco(KF_ACO, rKEQ_ACO, citrate, isocitrate),
        vIDH ~ _vidh3(KCAT_IDH, ET_IDH, isocitrate, KM_ISOC_IDH, adp_m, KM_ADP_IDH, ca_m, KM_CA_IDH, nad_m, KM_NAD_IDH, nadh_m, KI_NADH_IDH),
        vKGDH ~ _vkgdh(KCAT_KGDH, ET_KGDH, akg, KM_AKG_KGDH, nad_m, KM_NAD_KGDH, mg_m, KM_MG_KGDH, ca_m, KM_CA_KGDH, h_m, KH1_KGDH, KH2_KGDH, NI_AKG_KGDH),
        vSL ~ v_sl,
        vFH ~ _vfh(KF_FH, rKEQ_FH, fumarate, malate),
        vMDH ~ _vmdh_reverisble(VF_MDH, VR_MDH, oaa, KM_OAA_MDH, malate, KM_MAL_MDH, nad_m, KM_NAD_MDH, nadh_m, KM_NADH_MDH, h_m, K_OFFSET_MDH, KH1_MDH, KH2_MDH, KH3_MDH, KH4_MDH),
        vAAT ~ _vast_reverisble(KF_AAT, oaa, GLU, akg, ASP, KEQ_AAT),
        D(isocitrate) ~ vACO - vIDH,
        D(akg) ~ vIDH - vKGDH + vAAT,
        D(scoa) ~ vKGDH - vSL,
        D(succinate) ~ vSL - vSDH,
        D(fumarate) ~ vSDH - vFH,
        D(malate) ~ vFH - vMDH,
        D(oaa) ~ vMDH - vCS - vAAT,
    ]
    return (; eqs_tca, vIDH, vKGDH, vMDH, vSL, succinate, fumarate, oaa)
end

function get_tca_sys(; atp_m, adp_m, nad_m, nadh_m, ca_m, h_m=exp10(-7.6) * Molar, pi_m=8mM, mg_m=0.4mM, use_mg=false, name=:tcasys)
    @unpack eqs_tca = get_tca_eqs(; atp_m, adp_m, nad_m, nadh_m, ca_m, h_m, pi_m, mg_m, use_mg)
    return System(eqs_tca, t; name)
end
