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

		### IDH (Isocitrte dehydrogenase)
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
		KF_SL = 0.0284/ms
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
		KCAT_MDH = 0.1259/ms
		ET_MDH = 0.154mM
		KM_MAL_MDH = 1.493mM
		KI_OAA_MDH = 3.1E-3mM
		KM_NAD_MDH = 0.2244mM

		### AAT (alanine aminotransferase)
		KF_AAT = 0.0217/(ms*mM)
		KEQ_AAT = 6.6
		KASP_AAT = 1.5E-6/ms
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
		oaa(t) = 2.683e-8mM
		cit(t) # Conserved
		isoc(t) = 0.486mM		# isocitrate
		akg(t) = 0.00019mM		# alpha-ketoglutarate
		scoa(t) = 0.286mM		# succinyl-CoA
		suc(t) = 0.00011mM      # succinate
        fum(t) = 0.00528mM      # fumarate
		mal(t) = 0.002782mM		# malate
	end

	vcs = KCAT_CS * ET_CS * hil(ACCOA, KM_ACCOA_CS) * hil(oaa, KM_OAA_CS)

	vidh = let
		fa = hilr(adp_m, KM_ADP_IDH) * hilr(ca_m, KM_CA_IDH)
		f_isoc = NaNMath.pow(KM_ISOC_IDH / isoc, NI_ISOC_IDH)
		f_nad = (1 + nadh_m/KI_NADH_IDH) * (KM_NAD_IDH / nad_m)
		f_h = h_m / KH1_IDH + KH2_IDH / h_m
		vidh = KCAT_IDH * ET_IDH / (f_h + (1 + fa * f_isoc) * (1 + f_nad))
	end

	vkgdh = let
		f_a = hilr(mg_m, KM_MG_KGDH) * hilr(ca_m, KM_CA_KGDH)
		f_akg = NaNMath.pow(KM_AKG_KGDH / akg, NI_AKG_KGDH)
		f_nad = KM_NAD_KGDH / nad_m
		f_h = 1 + h_m / KH1_KGDH + KH2_KGDH / h_m
		vkgdh = ET_KGDH * KCAT_KGDH / (f_h + f_a * (f_akg + f_nad))
	end

	vmdh = let
		f_ha = K_OFFSET_MDH + hil(KH1_MDH * hil(KH2_MDH, h_m), h_m)
		f_hi = hil(h_m * hil(h_m, KH4_MDH), KH3_MDH)^2
		f_oaa = hilr(oaa, KI_OAA_MDH)
		f_mal = hil(mal * f_oaa, KM_MAL_MDH)
		f_nad = hil(nad_m, KM_NAD_MDH)
		vmdh = KCAT_MDH * ET_MDH * f_ha * f_hi * f_nad * f_mal
	end

	vsl = let
		if use_mg
			poly_atp = 1 + h_m/KA_ATP + mg_m/KMG_ATP
			poly_adp = 1 + h_m/KA_ADP + mg_m/KMG_ADP
			atp4 = atp_m / poly_atp
			hatp = atp4 * (h_m/KA_ATP)
			pi_poly = 1 + h_m/KA_PI
			suc_poly = 1 + h_m/KA_SUC
			keapp = KEQ * (poly_atp * suc_poly) / (poly_adp * pi_poly)
			KF_SL * (scoa * adp_m - suc * (atp4 + hatp) * COA / keapp)
		else
			KF_SL * (scoa * adp_m - suc * atp_m * COA / KEQ_SL)
		end
	end

	eqs = [
		TCA_T ~ cit + isoc + oaa + akg + scoa + suc + fum + mal,
		vCS ~ vcs,
		vACO ~ KF_ACO * (cit - isoc/KEQ_ACO),
		vIDH ~ vidh,
		vKGDH ~ vkgdh,
		vSL ~ vsl,
		vFH ~ KF_FH * (fum - mal / KEQ_FH),
		vMDH ~ vmdh,
		vAAT ~ KF_AAT * oaa * GLU * hil(KASP_AAT * KEQ_AAT, akg * KF_AAT),
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
