# TCA cycle model, default parameters values from Gauthier et al. (2013)
function get_tca_sys(; atp_m, adp_m, nad_m, nadh_m, ca_m, h_m=exp10(-7.6) * Molar, pi_m=8mM, mg_m=0.4mM, use_mg=false, name=:tcasys)
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
        VR_MDH = VF_MDH * KM_OAA_MDH * KM_NADH_MDH / (KEQ_MDH * KM_NAD_MDH * KM_MAL_MDH)
        ### AAT (alanine aminotransferase, aka AST)
        KF_AAT = 21.4Hz / mM
        KEQ_AAT = 6.6
        GLU = 10mM   ## Glutamate
        ASP = 100μM  ## Aspartate
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

    ## Citrate synthase
    v_cs = let
        vmax = KCAT_CS * ET_CS
        A = oaa / KM_OAA_CS
        B = ACCOA / KM_ACCOA_CS
        v_cs = vmax * A * B / (1 + A + B + A * B)
    end
    ## Aconitase
    v_aco = KF_ACO * (cit - isoc * rKEQ_ACO)
    ## IDH3 (Isocitrate dehydrogenase, NADH-producing)
    v_idh3 = let
        vmax = KCAT_IDH * ET_IDH
        A = (isoc / KM_ISOC_IDH)^2 * (1 + adp_m / KM_ADP_IDH) * (1 + ca_m / KM_CA_IDH)
        B = nad_m * KI_NADH_IDH / (KM_NAD_IDH * (KI_NADH_IDH + nadh_m))
        H = 1 + h_m / KH1_IDH + KH2_IDH / h_m
        v_idh3 = vmax * A * B / (H * A * B + A + B + 1)
    end
    ## Oxoglutarate dehydrogenase complex (OGDC / KGDH)
    v_kgdh = let
        vmax = ET_KGDH * KCAT_KGDH
        f_h = 1 + h_m / KH1_KGDH + KH2_KGDH / h_m
        f_akg = NaNMath.pow(akg / KM_AKG_KGDH, NI_AKG_KGDH)
        f_mgca = (1 + mg_m / KM_MG_KGDH) * (1 + ca_m / KM_CA_KGDH)
        f_nad = nad_m / KM_NAD_KGDH
        vmax = ET_KGDH * KCAT_KGDH
        v_kgdh = vmax * f_akg * f_nad * f_mgca / (f_h * f_mgca * f_akg * f_nad + f_akg + f_nad)
    end
    ## Succinyl-CoA lyase
    v_sl = let
        if use_mg
            atp4, hatp, _, atp_poly = breakdown_atp(atp_m, h_m, mg_m)
            _, _, _, adp_poly = breakdown_adp(adp_m, h_m, mg_m)
            pi_poly = pipoly(h_m)
            suc_poly = sucpoly(h_m)
            rkeq = rKEQ_SL * (adp_poly * pi_poly) / (atp_poly * suc_poly)
            atp = atp4 + hatp
        else
            rkeq = rKEQ_SL
            atp = atp_m
        end
        v_sl = KF_SL * (scoa * adp_m * pi_m - suc * atp * COA * rKEQ_SL)
    end
    ## Fumarate hydrase
    v_fh = KF_FH * (fum - mal * rKEQ_FH)
    ## Malate dehydrogenase
    v_mdh = let
        vmax = KCAT_MDH * ET_MDH
        f_ha = K_OFFSET_MDH + (KH1_MDH * KH2_MDH / (KH1_MDH * KH2_MDH + KH2_MDH * h_m + h_m^2))
        f_hi = (h_m^2 / (h_m^2 + h_m * KH3_MDH + KH3_MDH * KH4_MDH))^2
        f_h = f_ha * f_hi
        f_nad = nad_m / KM_NAD_MDH
        f_mal = mal / KM_MAL_MDH
        f_oaa = oaa / KM_OAA_MDH
        f_nadh = nadh_m / KM_NADH_MDH
        v_mdh = f_h * (VF_MDH * f_nad * f_mal - VR_MDH * f_oaa * f_nadh) / (1 + f_nad + f_nad * f_mal + f_oaa * f_nadh + f_nadh)
    end
    ## AST
    v_aat = KF_AAT * (oaa * GLU - akg * ASP / KEQ_AAT)

    eqs = [
        TCA_T ~ cit + isoc + oaa + akg + scoa + suc + fum + mal,
        vCS ~ v_cs,
        vACO ~ v_aco,
        vIDH ~ v_idh3,
        vKGDH ~ v_kgdh,
        vSL ~ v_sl,
        vFH ~ v_fh,
        vMDH ~ v_mdh,
        vAAT ~ v_aat,
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
