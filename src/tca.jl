"Citrate synthase rate"
function v_cs(oaa; accoa=10μM, KM_ACCOA_CS=12.6μM, KM_OAA_CS=0.64μM, KCAT_CS=50Hz, ET_CS=400μM)
    return KCAT_CS * ET_CS * hil(oaa, KM_OAA_CS) * hil(accoa, KM_ACCOA_CS)
end

"Aconitase rate"
function v_aco(cit, isoc; rKeq=inv(2.22), KF_ACO=12.5Hz)
    return KF_ACO * (cit - isoc * rKeq)
end

"IDH3 (Isocitrate dehydrogenase, NADH-producing) rate (Wei, 2011)"
function v_idh3(isoc, nad_m, nadh_m; ca_m, adp_m, h_m, KI_NADH_IDH=190μM, KCAT_IDH=60Hz, ET_IDH=109μM, KM_ISOC_IDH=1520μM, KM_NAD_IDH=923μM, KM_ADP_IDH=620μM, KM_CA_IDH=0.5μM, KH1_IDH=10nM, KH2_IDH=900nM)
    vmax = KCAT_IDH * ET_IDH
    a = (isoc / KM_ISOC_IDH)^2 * (1 + adp_m / KM_ADP_IDH) * (1 + ca_m / KM_CA_IDH)
    b = nad_m / KM_NAD_IDH * hil(KI_NADH_IDH, nadh_m)
    h = 1 + h_m / KH1_IDH + KH2_IDH / h_m
    return vmax * a * b / (h * a * b + a + b + 1)
end

"Oxoglutarate dehydrogenase complex (OGDC / KGDH) rate (Wei, 2011)"
function v_kgdh(akg, nad_m; ca_m, mg_m=400μM, h_m=exp10(-7.6) * Molar, ET_KGDH=500μM, KCAT_KGDH=50Hz, KM_AKG_KGDH=1940μM, KM_NAD_KGDH=38.7mM, KH1_KGDH=0.04μM, KH2_KGDH=0.07μM, KM_MG_KGDH=30.8μM, KM_CA_KGDH=0.15μM, NI_AKG_KGDH=1.2)
    vmax = ET_KGDH * KCAT_KGDH
    f_h = 1 + h_m / KH1_KGDH + KH2_KGDH / h_m
    f_akg = NaNMath.pow(akg / KM_AKG_KGDH, NI_AKG_KGDH)
    f_mgca = (1 + mg_m / KM_MG_KGDH) * (1 + ca_m / KM_CA_KGDH)
    f_nad = nad_m / KM_NAD_KGDH
    vmax = ET_KGDH * KCAT_KGDH
    return vmax * f_akg * f_nad * f_mgca / (f_h * f_mgca * f_akg * f_nad + f_akg + f_nad)
end

"Succinyl-CoA lyase rate"
function v_sl(scoa, adp_m, suc, atp_m; pi_m=8mM, COA=20μM, KF_SL=28 / (μM * ms), rKEQ_SL=inv(3.115))
    return KF_SL * (scoa * adp_m * pi_m - suc * atp_m * COA * rKEQ_SL)
end

function v_sl2(scoa, adp_m, suc, atp_m; pi_m=8mM, COA=20μM, KF_SL=28 / (μM * ms), rKEQ_SL=inv(3.115), h_m=exp10(-7.6) * Molar, mg_m=400μM)
    atp4, hatp, _, atp_poly = breakdown_atp(atp_m, h_m, mg_m)
    _, _, _, adp_poly = breakdown_adp(adp_m, h_m, mg_m)
    pi_poly = pipoly(h_m)
    suc_poly = sucpoly(h_m)
    rkeq = rKEQ_SL * (adp_poly * pi_poly) / (atp_poly * suc_poly)
    atp = atp4 + hatp
    return v_sl(scoa, adp_m, suc, atp; pi_m, COA, KF_SL, rKEQ_SL=rkeq)
end

"Fumarate hydrase rate"
function v_fh(fum, mal; KF_FH=8.3Hz, rKEQ_FH=1.0)
    KF_FH * (fum - mal * rKEQ_FH)
end

"Malate dehydrogenase rate"
function v_mdh(mal, oaa, nad_m; h_m=exp10(-7.6) * Molar, KCAT_MDH=126Hz, ET_MDH=154μM, KH1_MDH=1.131E-5mM, KH2_MDH=26.7mM, KH3_MDH=6.68E-9mM, KH4_MDH=5.62E-6mM, K_OFFSET_MDH=3.99E-2, KM_MAL_MDH=1493μM, KI_OAA_MDH=3.1μM, KM_NAD_MDH=224.4μM)
    vmax = KCAT_MDH * ET_MDH
    f_ha = K_OFFSET_MDH + hil(KH1_MDH * hil(KH2_MDH, h_m), h_m)
    f_hi = 1 + KH3_MDH / h_m * (1 + KH4_MDH / h_m)
    f_oaa = hil(KI_OAA_MDH, oaa)
    f_mal = hil(mal * f_oaa, KM_MAL_MDH)
    f_nad = hil(nad_m, KM_NAD_MDH)
    return vmax * ET_MDH * f_ha * f_hi * f_nad * f_mal
end

# TCA cycle model, default parameters values from Gauthier et al. (2013)
function get_tca_sys(; atp_m, adp_m, nad_m, nadh_m, h_m=exp10(-7.6) * Molar, ca_m, pi_m=8mM, mg_m=0.4mM, use_mg=false, name=:tcasys)
    @parameters begin
        # Total TCA metabolite pool
        TCA_T = 1.3mM

        ## Citrate synthase
        KCAT_CS = 50Hz
        ET_CS = 400μM
        KM_ACCOA_CS = 12.6μM
        KM_OAA_CS = 0.64μM
        ACCOA = 1000μM

        ## ACO (aconitase)
        KF_ACO = 12.5Hz # Zhou, 2009
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

        ### SL (Succinyl-coA lyase)
        COA = 20μM
        KF_SL = 28 / (μM * ms)
        rKEQ_SL = inv(3.115)

        ### FH (Fumarate hydrase) parameters
        KF_FH = 8.3Hz
        rKEQ_FH = 1.0

        ### MDH (Malate dehydrogenase)
        KCAT_MDH = 126Hz
        ET_MDH = 154μM
        KH1_MDH = 11.31nM
        KH2_MDH = 26.7mM
        KH3_MDH = 6.68E-3nM
        KH4_MDH = 5.62nM
        K_OFFSET_MDH = 0.0399
        KM_MAL_MDH = 1493μM
        KI_OAA_MDH = 3.1μM
        KM_NAD_MDH = 224.4μM

        ### AAT (alanine aminotransferase)
        KF_AAT = 0.644Hz / mM
        KEQ_AAT = 6.6
        K_ASP = 1.5e-3Hz
        GLU = 10mM        # Glutamate
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
        oaa(t) = 11.6μM
        cit(t) # Conserved
        isoc(t) = 51.6μM    # isocitrate
        akg(t) = 51μM       # alpha-ketoglutarate
        scoa(t) = 35μM      # succinyl-CoA
        suc(t) = 1.9μM      # succinate
        fum(t) = 175μM      # fumarate
        mal(t) = 160μM      # malate
    end

    eqs = [
        TCA_T ~ cit + isoc + oaa + akg + scoa + suc + fum + mal,
        vCS ~ v_cs(oaa; accoa=ACCOA, KCAT_CS, ET_CS, KM_ACCOA_CS, KM_OAA_CS),
        vACO ~ v_aco(cit, isoc; rKeq=rKEQ_ACO, KF_ACO),
        vIDH ~ v_idh3(isoc, nad_m, nadh_m; ca_m, adp_m, h_m, KI_NADH_IDH, KCAT_IDH, ET_IDH, KH1_IDH, KH2_IDH, KM_ISOC_IDH, KM_NAD_IDH, KM_ADP_IDH, KM_CA_IDH),
        vKGDH ~ v_kgdh(akg, nad_m; ca_m, mg_m, h_m, ET_KGDH, KCAT_KGDH, KM_AKG_KGDH, KM_NAD_KGDH, KH1_KGDH, KH2_KGDH, KM_MG_KGDH, KM_CA_KGDH, NI_AKG_KGDH),
        vSL ~ use_mg ? v_sl2(scoa, adp_m, suc, atp_m; pi_m, COA, KF_SL, rKEQ_SL, h_m, mg_m) : v_sl(scoa, adp_m, suc, atp_m; pi_m, COA, KF_SL, rKEQ_SL),
        vFH ~ v_fh(fum, mal; KF_FH, rKEQ_FH),
        vMDH ~ v_mdh(mal, oaa, nad_m; h_m, KCAT_MDH, ET_MDH, KH1_MDH, KH2_MDH, KH3_MDH, KH4_MDH, K_OFFSET_MDH, KM_MAL_MDH, KI_OAA_MDH, KM_NAD_MDH),
        vAAT ~ KF_AAT * oaa * GLU * hil(K_ASP * KEQ_AAT, akg * KF_AAT),
        D(isoc) ~ vACO - vIDH,
        D(akg) ~ vIDH - vKGDH + vAAT,
        D(scoa) ~ vKGDH - vSL,
        D(suc) ~ vSL - vSDH,
        D(fum) ~ vSDH - vFH,
        D(mal) ~ vFH - vMDH,
        D(oaa) ~ vMDH - vCS - vAAT,
    ]
    return System(eqs, t; name)
end
