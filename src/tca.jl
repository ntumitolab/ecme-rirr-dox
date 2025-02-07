# TCA cycle model, default parameters values from Gauthier et al. (2013)
function get_tca_sys(atp_m, adp_m, nad_m, nadh_m, h_m, ca_m, mg_m; use_mg=false, name=:tcasys)
    @parameters begin
        TCA_T = 1.3mM  # Total TCA metabolite pool (mM)
        ### Citrate synthase
        # KCAT = 1.5891E-4 # Gauthier (2013), cellular model
        KCAT_CS = 0.23523Hz # Gauthier (2013), mitochondrial model
        ET_CS = 0.4mM
        KM_ACCOA_CS = 0.0126mM
        KM_OAA_CS = 6.4E-4mM
        # ACCOA = 1.0  # Li, 2015
        ACCOA = 0.1mM  # De Oliveira (2016)

        ### ACO (aconitase)
        KF_ACO = 0.1Hz # Adjusted
        KEQ_ACO = 2.22

        ### IDH3 (Isocitrate dehydrogenase, NADH-producing)
        KI_NADH_IDH = 0.19mM
        KCAT_IDH = 535Hz
        ET_IDH = 0.109mM
        KH1_IDH = 1E-5mM
        KH2_IDH = 9E-4mM
        KM_ISOC_IDH = 1.52mM
        NI_ISOC_IDH = 2
        KM_NAD_IDH = 0.923mM
        KM_ADP_IDH = 0.62mM
        KM_CA_IDH = 5E-4mM

        ### KGDH (alpha-ketoglutarate dehydrogenase)
        ET_KGDH = 0.5mM
        KCAT_KGDH = 17.9Hz
        KM_AKG_KGDH = 30mM
        KM_NAD_KGDH = 38.7mM
        KH1_KGDH = 4E-5mM
        KH2_KGDH = 7E-5mM
        KM_MG_KGDH = 0.0308mM
        KM_CA_KGDH = 1.5E-4mM
        NI_AKG_KGDH = 1.2

        ### SL (Succinyl-coA lyase)
        COA = 0.02mM
        KF_SL = 0.0284 / (mM * ms)
        KEQ_SL = 3.115mM

        ### FH (Fumarate hydrase) parameters
        KF_FH = 8.3Hz
        KEQ_FH = 1.0

        ### MDH (Malate dehydrogenase)
        KH1_MDH = 1.131E-5mM
        KH2_MDH = 26.7mM
        KH3_MDH = 6.68E-9mM
        KH4_MDH = 5.62E-6mM
        K_OFFSET_MDH = 3.99E-2
        KCAT_MDH = 0.1259 / ms
        ET_MDH = 0.154mM
        KM_MAL_MDH = 1.493mM
        KI_OAA_MDH = 3.1E-3mM
        KM_NAD_MDH = 0.2244mM

        ### AAT (alanine aminotransferase)
        KF_AAT = 0.0217 / (ms * mM)
        KEQ_AAT = 6.6
        KASP_AAT = 1.5E-6 / ms
        GLU = 30.0mM
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
        oaa(t) = 0.011576938766421891mM
        cit(t) # Conserved
        isoc(t) = 0.05159054318098895mM# isocitrate
        akg(t) = 0.051145197718677655mM# alpha-ketoglutarate
        scoa(t) = 0.03508849487000582mM# succinyl-CoA
        suc(t) = 0.0019107469302081612mM    # succinate
        fum(t) = 0.1751906841603877mM       # fumarate
        mal(t) = 0.15856757152954906mM# malate
    end

    v_cs = KCAT_CS * ET_CS * hil(ACCOA, KM_ACCOA_CS) * hil(oaa, KM_OAA_CS)
    v_aco = KF_ACO * (cit - isoc / KEQ_ACO)

    v_kgdh = let
        f_mgca = (1 + mg_m / KM_MG_KGDH) * (1 + ca_m / KM_CA_KGDH)
        f_h = 1 + h_m / KH1_KGDH + KH2_KGDH / h_m
        f_akg = NaNMath.pow(akg / KM_AKG_KGDH, NI_AKG_KGDH)
        f_nad = nad_m / KM_NAD_KGDH
        vmax = ET_KGDH * KCAT_KGDH
        v = vmax * f_akg * f_nad * f_mgca / (f_h * f_mgca * f_akg * f_nad + f_akg + f_nad)
    end

    v_idh = let
        vmax = KCAT_IDH * ET_IDH
        a = NaNMath.pow(isoc / KM_ISOC_IDH , NI_ISOC_IDH) * (1 + adp_m / KM_ADP_IDH) * (1 + ca_m / KM_CA_IDH)
        b = nad_m / KM_NAD_IDH * hil(KI_NADH_IDH, nadh_m)
        h = 1 + h_m / KH1_IDH + KH2_IDH / h_m
        vidh = vmax * a * b / (h * a * b + a + b + 1)
    end

    v_mdh = let
        f_ha = K_OFFSET_MDH + hil(KH1_MDH * hil(KH2_MDH, h_m), h_m)
        f_hi = hil(h_m * hil(h_m, KH4_MDH), KH3_MDH)^2
        f_oaa = hil(KI_OAA_MDH, oaa)
        f_mal = hil(mal * f_oaa, KM_MAL_MDH)
        f_nad = hil(nad_m, KM_NAD_MDH)
        vmdh = KCAT_MDH * ET_MDH * f_ha * f_hi * f_nad * f_mal
    end

    v_sl = let
        if use_mg
            atp4, hatp, _, poly_atp = breakdown_atp(atp_m, h_m, mg_m)
            _, _, _, poly_adp = breakdown_adp(adp_m, h_m, mg_m)
            pi_poly = pipoly(h_m)
            suc_poly = sucpoly(h_m)
            keapp = KEQ * (poly_atp * suc_poly) / (poly_adp * pi_poly)
            KF_SL * (scoa * adp_m - suc * (atp4 + hatp) * COA / keapp)
        else
            KF_SL * (scoa * adp_m - suc * atp_m * COA / KEQ_SL)
        end
    end

    v_fh = KF_FH * (fum - mal / KEQ_FH)
    v_aat = KF_AAT * oaa * GLU * hil(KASP_AAT * KEQ_AAT, akg * KF_AAT)

    eqs = [
        TCA_T ~ cit + isoc + oaa + akg + scoa + suc + fum + mal,
        vCS ~ v_cs,
        vACO ~ v_aco ,
        vIDH ~ v_idh,
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
    return ODESystem(eqs, t; name)
end
