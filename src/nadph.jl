# Mitochondrial NADPH generation (Kembro 2013, Gauthier 2013)
get_default_thd_params() = ComponentArray(
    iKM_NADH_THD = 10μM,  # Michealis constant for NADH of THD
    iKM_NAD_THD = 125μM,  # Michealis constant for NAD of THD
    iKM_NADP_THD = 20μM,  # Michealis constant for NADP of THD
    iKM_NADPH_THD = 20μM,  # Michealis constant for nadph_m of THD
    ET_THD = 1μM,  # THD concentration
    KF_THD = 1.17474Hz,  # Max forward rate of THD
    KB_THD = 17.2756Hz,  # Max backward rate of THD
    D_THD = 0.5,  # Voltage assymetry factor of THD
    X_THD = 0.1,  # Voltage dependence factor of THD
)

"NADPH-NADH transhydrogenase (THD) rate"
function vTHD(; nadh_m, nad_m, nadp_m, nadph_m, ΔμH, p=get_default_thd_params())
    @unpack iKM_NADH_THD, iKM_NAD_THD, iKM_NADP_THD, iKM_NADPH_THD, ET_THD, KF_THD, KB_THD, D_THD, X_THD = p
    # Corrected an error from both the paper and code from Gauthier et al.
    vNAD = exp(iVT * X_THD * (1.0 - D_THD) * ΔμH)
    fNAD = nad_m * iKM_NAD_THD
    fvNAD = fNAD * vNAD
    fNADH = nadh_m * iKM_NADH_THD
    fNADP = nadp_m * iKM_NADP_THD
    vNADP = exp(iVT * X_THD * D_THD * ΔμH)
    fvNADP = vNADP * fNADP
    fNADPH = nadph_m * iKM_NADPH_THD
    denom = 1 + fNAD + fNADH + fNADP + fNADPH + (fNADPH + fvNAD) * (fNADPH + fvNADP)
    num = KF_THD * fNADH * fvNADP - KB_THD * fvNAD * fNADPH
    return ET_THD * num / denom
end

get_default_idh2_params() = ComponentArray(
    KM_H_IDH2 = 5μM,        # Dissociation constant for H+ of IDH2
    KM_ISOC_IDH2 = 45μM,    # Dissociation constant for isocitrate of IDH2
    KM_NADP_IDH2 = 46μM,    # Michealis constant for NADP of IDH2
    KI_NADP_IDH2 = 2E-6mM,  # Inhibition constant for NADP of IDH2
    KM_NADPH_IDH2 = 12μM,   # Michealis constant for nadph_m of IDH2
    KM_AKG_IDH2 = 80μM,     # Michealis constant for AKG of IDH2
    VF_IDH2 = 87.2mM/Hz,    # Max forward rate of IDH2
    VB_IDH2 = 5.45mM/Hz,    # Max backward rate of IDH2
)

function vIDH2(; h_m, isoc, akg, nadp_m, nadph_m, p=get_default_idh2_params())
    @unpack KM_H_IDH2, KM_ISOC_IDH2, KM_NADP_IDH2, KI_NADP_IDH2, KM_NADPH_IDH2, KM_AKG_IDH2, VF_IDH2, VB_IDH2 = p
    f_h = hil(KM_H_IDH2, h_m)
    f_isoc = isoc / KM_ISOC_IDH2
    f_nadp = nadp_m / (KM_NADP_IDH2 * hil(nadp_m, KI_NADP_IDH2))
    f_akg = akg / KM_AKG_IDH2
    f_nadph = nadph_m / KM_NADPH_IDH2
    denom = (1 + f_isoc + f_akg) * (1 + f_nadp + f_nadph)
    num = (VF_IDH2 * f_isoc * f_nadp - VB_IDH2 * f_akg * f_nadph)
    return num * f_h / denom
end
