# TCA cycle model, default parameters values from Gauthier et al. (2013)
function get_tca_sys(; atp_m, adp_m, nad_m, nadh_m, h_m, ca_m, pi_m=8mM, mg_m=0.4mM, use_mg=false, name=:tcasys)
    @parameters begin
        TCA_T = 1300μM  # Total TCA metabolite pool

        ## Citrate synthase
        KCAT_CS = 0.23523Hz # Gauthier (2013), mitochondrial model
        ET_CS = 400μM
        KM_ACCOA_CS = 12.6μM
        KM_OAA_CS = 0.64μM
        ACCOA = 1000μM  # 100μM

        ## ACO (aconitase)
        KF_ACO = 12.5Hz # Zhou, 2009
        KEQ_ACO = 2.22

        ## IDH3 (Isocitrate dehydrogenase, NADH-producing)
        KI_NADH_IDH = 190μM
        KCAT_IDH = 41Hz
        ET_IDH = 109μM
        KH1_IDH = 1E-5mM
        KH2_IDH = 9E-4mM
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
        KH1_KGDH = 0.04μM
        KH2_KGDH = 0.07μM
        KM_MG_KGDH = 30.8μM
        KM_CA_KGDH = 0.15μM
        NI_AKG_KGDH = 1.2

        ### SL (Succinyl-coA lyase)
        COA = 20μM
        KF_SL = 28 / (μM * ms)
        KEQ_SL = 3.115

        ### FH (Fumarate hydrase) parameters
        KF_FH = 8.3Hz
        KEQ_FH = 1.0

        ### MDH (Malate dehydrogenase)
        KCAT_MDH = 126Hz
        ET_MDH = 154μM
        KH1_MDH = 1.131E-5mM
        KH2_MDH = 26.7mM
        KH3_MDH = 6.68E-9mM
        KH4_MDH = 5.62E-6mM
        K_OFFSET_MDH = 3.99E-2
        KM_MAL_MDH = 1493μM
        KI_OAA_MDH = 3.1μM
        KM_NAD_MDH = 224.4μM

        ### AAT (alanine aminotransferase)
        KF_AAT = 0.644Hz / mM
        KEQ_AAT = 6.6
        GLU = 10mM        # Glutamate
        ASP = GLU         # Aspartate
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

    v_kgdh = let
        f_mgca = (1 + mg_m / KM_MG_KGDH) * (1 + ca_m / KM_CA_KGDH)
        f_h = 1 + h_m / KH1_KGDH + KH2_KGDH / h_m
        f_akg = NaNMath.pow(akg / KM_AKG_KGDH, NI_AKG_KGDH)
        f_nad = nad_m / KM_NAD_KGDH
        vmax = ET_KGDH * KCAT_KGDH
        vmax * f_akg * f_nad * f_mgca / (f_h * f_mgca * f_akg * f_nad + f_akg + f_nad)
    end

    v_idh = let
        vmax = KCAT_IDH * ET_IDH
        a = NaNMath.pow(isoc / KM_ISOC_IDH, NI_ISOC_IDH) * (1 + adp_m / KM_ADP_IDH) * (1 + ca_m / KM_CA_IDH)
        b = nad_m / KM_NAD_IDH * hil(KI_NADH_IDH, nadh_m)
        h = 1 + h_m / KH1_IDH + KH2_IDH / h_m
        vmax * a * b / (h * a * b + a + b + 1)
    end

    v_mdh = let
        vmax = KCAT_MDH * ET_MDH
        f_ha = K_OFFSET_MDH + hil(KH1_MDH * hil(KH2_MDH, h_m), h_m)
        f_hi = hil(h_m * hil(h_m, KH4_MDH), KH3_MDH)^2
        f_oaa = hil(KI_OAA_MDH, oaa)
        f_mal = hil(mal * f_oaa, KM_MAL_MDH)
        f_nad = hil(nad_m, KM_NAD_MDH)
        vmax * ET_MDH * f_ha * f_hi * f_nad * f_mal
    end

    v_sl = let
        if use_mg
            atp4, hatp, _, atp_poly = breakdown_atp(atp_m, h_m, mg_m)
            _, _, _, adp_poly = breakdown_adp(adp_m, h_m, mg_m)
            pi_poly = pipoly(h_m)
            suc_poly = sucpoly(h_m)
            keq = KEQ_SL * (atp_poly * suc_poly) / (adp_poly * pi_poly)
            atp = atp4 + hatp
        else
            keq = KEQ_SL
            atp = atp_m
        end
        KF_SL * (scoa * adp_m * pi_m - suc * atp * COA / keq)
    end

    eqs = [
        TCA_T ~ cit + isoc + oaa + akg + scoa + suc + fum + mal,
        vCS ~ KCAT_CS * ET_CS * hil(ACCOA, KM_ACCOA_CS) * hil(oaa, KM_OAA_CS),
        vACO ~ KF_ACO * (cit - isoc / KEQ_ACO),
        vIDH ~ v_idh,
        vKGDH ~ v_kgdh,
        vSL ~ v_sl,
        vFH ~ KF_FH * (fum - mal / KEQ_FH),
        vMDH ~ v_mdh,
        vAAT ~ KF_AAT * (oaa * GLU - akg * ASP / KEQ_AAT),
        D(isoc) ~ vACO - vIDH,
        D(akg) ~ vIDH - vKGDH + vAAT,
        D(scoa) ~ vKGDH - vSL,
        D(suc) ~ vSL - vSDH,
        D(fum) ~ vSDH - vFH,
        D(mal) ~ vFH - vMDH,
        D(oaa) ~ vMDH - vCS - vAAT,
    ]
    return ODESystem(eqs, t; name)
end
