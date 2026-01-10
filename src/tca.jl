# TCA cycle model
"Citrate concentration from total TCA pool"
_cit(isoc, oaa, akg, scoa, suc, fum, mal, TCA_T=1.3mM) = TCA_T - isoc - oaa - akg - scoa - suc - fum - mal
_cit(u, p, t) = _cit(u.isoc, u.oaa, u.akg, u.scoa, u.suc, u.fum, u.mal, p.TCA_T)

get_default_CS_params() = ComponentArray(
    KCAT_CS = 50Hz,
    ET_CS = 400μM,
    iKM_ACCOA_CS = inv(12.6μM),
    iKM_OAA_CS = inv(0.64μM),
    ACCOA = 1000μM,
)

"Citrate synthase rate"
function vCS(; oaa, p=get_default_CS_params())
    vmax = p.KCAT_CS * p.ET_CS
    A = oaa * p.iKM_OAA_CS
    B = p.ACCOA * p.iKM_ACCOA_CS
    return vmax * A * B / ((1 + A) * (1 + B))
end

get_default_ACO_params() = ComponentArray(
    TCA_T = 1.3mM,
    KF_ACO = 12.5Hz, # Zhou, 2009
    rKEQ_ACO = inv(2.22),
)

"Aconitase rate"
vACO(; cit, isoc, p=get_default_ACO_params()) = p.KF_ACO * (cit - isoc * p.rKEQ_ACO)

get_default_IDH3_params() = ComponentArray(
    ET_IDH = 109μM,
    KCAT_IDH = 43Hz,
    KI_NADH_IDH = 190μM,
    KH1_IDH = 10nM,
    KH2_IDH = 900nM,
    KM_ISOC_IDH = 1520μM,
    ## NI_ISOC_IDH = 2,
    iKM_NAD_IDH = inv(923μM),
    iKM_ADP_IDH = inv(620μM),
    iKM_CA_IDH = inv(0.5μM),
)

"Isocitrate dehydrogenase (NADH-producing) rate"
function vIDH3(; isoc, nad_m, nadh_m, adp_m, ca_m, h_m, p=get_default_IDH3_params())
    @unpack ET_IDH, KCAT_IDH, KI_NADH_IDH, KH1_IDH, KH2_IDH, iKM_ISOC_IDH, iKM_CA_IDH, iKM_ADP_IDH, NI_ISOC_IDH = p
    vmax = KCAT_IDH * ET_IDH
    A = (isoc * iKM_ISOC_IDH)^2 * (1 + adp_m * iKM_ADP_IDH) * (1 + ca_m * iKM_CA_IDH)
    B = nad_m * KI_NADH_IDH / (KM_NAD_IDH * (KI_NADH_IDH + nadh_m))
    H = 1 + h_m / KH1_IDH + KH2_IDH / h_m
    return vmax * A * B / (H * A * B + A + B + 1)
end

get_default_KGDH_params() = ComponentArray(
    ET_KGDH = 500μM,
    KCAT_KGDH = 50Hz,
    iKM_AKG_KGDH = inv(1940μM),
    iKM_NAD_KGDH = inv(38.7mM),
    KH1_KGDH = 40nM,
    KH2_KGDH = 70nM,
    iKM_MG_KGDH = inv(30.8μM),
    iKM_CA_KGDH = inv(0.15μM),
    NI_AKG_KGDH = 1.2,
)

"Alpha-ketoglutarate dehydrogenase complex (OGDC / KGDH) rate"
function vKGDH(; akg, nad_m, h_m, mg_m, ca_m, p=get_default_KGDH_params())
    @unpack ET_KGDH, KCAT_KGDH, iKM_AKG_KGDH, iKM_NAD_KGDH, KH1_KGDH, KH2_KGDH, iKM_MG_KGDH, iKM_CA_KGDH, NI_AKG_KGDH = p
    vmax = ET_KGDH * KCAT_KGDH
    f_h = 1 + h_m / KH1_KGDH + KH2_KGDH / h_m
    f_akg = NaNMath.pow(akg * iKM_AKG_KGDH, NI_AKG_KGDH)
    f_mgca = (1 + mg_m * iKM_MG_KGDH) * (1 + ca_m * iKM_CA_KGDH)
    f_nad = nad_m * iKM_NAD_KGDH
    vmax = ET_KGDH * KCAT_KGDH
    return vmax * f_akg * f_nad * f_mgca / (f_h * f_mgca * f_akg * f_nad + f_akg + f_nad)
end

get_default_vSL_params() = ComponentArray(
    COA = 20μM,
    KF_SL = 28Hz / mM^2,  ## Gauthier, 2013 vs 0.644 in Zhou, 2009
    rKEQ_SL = inv(3.115),
)

"Succinyl-CoA lyase rate"
function vSL(; scoa, suc, atp_m, adp_m, pi_m, p)
    @unpack COA, KF_SL, rKEQ_SL = p
    return KF_SL * (scoa * adp_m * pi_m - suc * atp_m * COA * rKEQ_SL)
end

"Succinyl-CoA lyase rate with magnesium"
function vSL(; scoa, suc, atp_m, adp_m, pi_m, h_m, mg_m, p)
    @unpack COA, KF_SL, rKEQ_SL = p
    atp4, hatp, _, atp_poly = breakdown_atp(atp_m, h_m, mg_m)
    _, _, _, adp_poly = breakdown_adp(adp_m, h_m, mg_m)
    pi_poly = pipoly(h_m)
    suc_poly = sucpoly(h_m)
    rkeq = rKEQ_SL * (adp_poly * pi_poly) / (atp_poly * suc_poly)
    atp = atp4 + hatp
    return KF_SL * (scoa * adp_m * pi_m - suc * atp * COA * rkeq)
end

get_default_FH_params() = ComponentArray(
    KF_FH = 8.3Hz,
    rKEQ_FH = 1.0,
)

"Fumarate hydrase rate"
vFH(; fum, mal, p=get_default_FH_params()) = p.KF_FH * (fum - mal * p.rKEQ_FH)

get_default_MDH_params() = ComponentArray(
    KCAT_MDH = 126Hz,
    ET_MDH = 154μM,
    KEQ_MDH = 3.08e-4,
    KH1_MDH = 11.31nM,
    KH2_MDH = 26.7mM,
    KH3_MDH = 6.68E-3nM,
    KH4_MDH = 5.62nM,
    K_OFFSET_MDH = 0.0399,
    iKM_NAD_MDH = inv(110μM),
    iKM_MAL_MDH = inv(450μM),
    iKM_OAA_MDH = inv(7μM),
    iKM_NADH_MDH = inv(17μM),
)

"Malate dehydrogenase (MDH) rate. Reversible."
function vMDH(; mal, oaa, nad_m, nadh_m, h_m, p=get_default_MDH_params())
    @unpack KCAT_MDH, ET_MDH, KEQ_MDH, KH1_MDH, KH2_MDH, KH3_MDH, KH4_MDH, K_OFFSET_MDH, iKM_NAD_MDH, iKM_MAL_MDH, iKM_OAA_MDH, iKM_NADH_MDH = p
    VF_MDH = KCAT_MDH * ET_MDH
    VR_MDH = VF_MDH * (iKM_NAD_MDH * iKM_MAL_MDH) / (iKM_OAA_MDH * iKM_NADH_MDH * KEQ_MDH)
    f_ha = K_OFFSET_MDH + (KH1_MDH * KH2_MDH / (KH1_MDH * KH2_MDH + KH2_MDH * h_m + h_m^2))
    f_hi = (h_m^2 / (h_m^2 + h_m * KH3_MDH + KH3_MDH * KH4_MDH))^2
    f_h = f_ha * f_hi
    f_nad = nad_m * iKM_NAD_MDH
    f_mal = mal * iKM_MAL_MDH
    f_oaa = oaa * iKM_OAA_MDH
    f_nadh = nadh_m * iKM_NADH_MDH
    return f_h * (VF_MDH * f_nad * f_mal - VR_MDH * f_oaa * f_nadh) / (1 + f_nad + f_nad * f_mal + f_oaa * f_nadh + f_nadh)
end

get_default_AAT_params() = ComponentArray(
    KF_AAT = 21.4Hz / mM,
    KEQ_AAT = 6.6,
    GLU = 10mM,
    ASP = 100μM,
)

function vAAT(; oaa, akg, p=get_default_AAT_params())
    @unpack KF_AAT, KEQ_AAT, GLU, ASP = p
    return KF_AAT * (oaa * GLU - akg * ASP / KEQ_AAT)
end

function tca_rates(u, p, t)
    @unpack isoc, akg, scoa, suc, fum, oaa, nad_m, nadh_m, adp_m, atp_m, ca_m, h_m, mg_m, pi_m = u
    cit = _cit(u, p, t)
    v_cs = vCS(; oaa, p)
    v_aco = vACO(; cit, isoc, p)
    v_idh3 = vIDH3(; isoc, nad_m, nadh_m, adp_m, ca_m, h_m, p)
    v_kgdh = vKGDH(; akg, nad_m, h_m, mg_m, ca_m, p)
    v_sl = vSL(; scoa, suc, atp_m, adp_m, pi_m, h_m, mg_m, p)
    v_fh = vFH(; fum, mal, p)
    v_mdh = vMDH(; mal, oaa, nad_m, nadh_m, h_m, p)
    v_aat = vAAT(; oaa, akg, p)
    return (; v_cs, v_aco, v_idh3, v_kgdh, v_sl, v_fh, v_mdh, v_aat)
end
