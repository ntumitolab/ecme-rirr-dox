"NADPH generation (Gauthier 2013)"
function get_nadph_sys(h_m, isoc, akg, ΔμH; name=:nadphsys)
    @parameters begin
		KM_H_IDH2 = 0.5mM  # Dissociation constant for H+ of IDH2
		KM_ISOC_IDH2 = 0.045mM  # Dissociation constant for isocitrate of IDH2
		KM_NADP_IDH2 = 0.046mM  # Michealis constant for NADP of IDH2
		KI_NADP_IDH2 = 2E-6mM  # Inhibition constant for NADP of IDH2
		KM_NADPH_IDH2 = 1.2E-2mM  # Michealis constant for nadph_m of IDH2
		KM_AKG_IDH2 = 0.08mM  # Michealis constant for AKG of IDH2
		VF_IDH2 = 8.72E-2/ms  # Max forward rate of IDH2
		VB_IDH2 = 5.45E-3/ms  # Max backward rate of IDH2
		KM_NADH_THD = 0.01  # Michealis constant for NADH of THD
		KM_NAD_THD = 0.125  # Michealis constant for NAD of THD
		KM_NADP_THD = 0.02  # Michealis constant for NADP of THD
		KM_NADPH_THD = 0.02  # Michealis constant for nadph_m of THD
		ET_THD = 1E-3  # THD concnetration (mM)
		KF_THD = 1.17474  # Max forward rate of THD
		KB_THD = 17.2756  # Max backward rate of THD
		D_THD = 0.5  # Voltage assymetry factor of THD
		X_THD = 0.1  # Voltage dependence factor of THD
	end

    @variables begin
        nadph_m(t)
        nadp_m(t)
        vIDH2(t)
        vTHD(t)
        nadh_m(t)
        nad_m(t)
    end

    vidh2 = let
        f_h = hil(KM_H_IDH2, h_m)
        f_isoc = isoc / KM_ISOC_IDH2
        f_nadp = nadp_m / (KM_NADP_IDH2 * hil(nadp_m, KI_NADP_IDH2))
        f_akg = akg / KM_AKG_IDH2
        f_nadph = nadph_m / KM_NADPH_IDH2
        denom = (1 + f_isoc + f_akg) * (1 + f_nadp + f_nadph)
        num = (VF_IDH2 * f_isoc * f_nadp - VB_IDH2 * f_akg * f_nadph)
        vidh2 = num * f_h / denom
    end

    vthd = let
        # Corrected an error from both the paper and code from Gauthier et al.
        vNAD = exp(iVT * X_THD * (1.0 - D_THD) * ΔμH)
        fNAD = nad_m / KM_NAD_THD
        fvNAD = fNAD * vNAD
        fNADH = nadh_m / KM_NADH_THD
        fNADP = nadp_m / KM_NADP_THD
        vNADP = exp(iVT * X_THD * D_THD * ΔμH)
        fvNADP = vNADP * fNADP
        fNADPH = nadph_m / KM_NADPH_THD
        denom = 1 + fNAD + fNADH + fNADP + fNADPH + (fNADPH + fvNAD) * (fNADPH + fvNADP)
        num = KF_THD * fNADH * fvNADP - KB_THD * fvNAD * fNADPH
        vthd = ET_THD * num / denom
    end

    return ODESystem([vIDH2 ~ vidh2, vTHD ~ vthd], t; name)
end
