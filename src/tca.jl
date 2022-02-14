# TCA cycle model, default parameters values from Gauthier et al. (2013)
using Parameters
include("common.jl")
include("concentrations.jl")

# CS (Citrate synthase)
@with_kw struct CSParams
	# KCAT = 1.5891E-4 # Gauthier (2013), cellular model
	KCAT = 2.3523E-4 # Gauthier (2013), mitochondrial model
	# KCAT = 0.05 # Li, 2014
    ET = 0.4
    VMAX = KCAT * ET
    KM_ACCOA = 0.0126
    KM_OAA = 6.4E-4
	ACCOA = 0.1  # De Oliveira (2016)
	# ACCOA = 1.0  # Li, 2015
    F_ACCOA = _mm(ACCOA, KM_ACCOA)
end

# Citrate synthase rate (scalar)
vcs(oaa, VMAX, F_ACCOA, KM_OAA) = VMAX * F_ACCOA * _mm(oaa, KM_OAA)
function vcs(oaa, pCS::CSParams)
	@unpack  VMAX, F_ACCOA, KM_OAA = pCS
	vCS = vcs(oaa, VMAX, F_ACCOA, KM_OAA)
end

# ACO (aconitase)
@with_kw struct ACOParams
    KF = 7.8959E-5
    KEQ = 2.22
	TCA_T = 1.3  # Total TCA metabolite pool (mM)
end

# Aconitase rate (scalar and inplace vector version)
vaco(cit, isoc, KF, KEQ) = KF * (cit - isoc / KEQ)
function vaco(isoc, akg, scoa, suc, fum, mal, oaa, pACO::ACOParams)
	@unpack KF, KEQ, TCA_T = pACO
	cit = TCA_T - isoc - akg - scoa - suc - fum - mal - oaa
	vACO = vaco(cit, isoc, KF, KEQ)
end

# IDH (Isocitrte dehydrogenase)
@with_kw struct IDHParams
    KI_NADH = 0.19
    KCAT = 0.5350
    ET = 0.109
    VMAX = KCAT * ET
    KH1 = 1E-5
    KH2 = 9E-4
    KM_ISOC = 1.52
    NI_ISOC = 2
    KM_NAD = 0.923
    KM_ADP = 0.62
    KM_CA = 5E-4
end

# Isocitrte dehydrogenase (IDH) rate
function vidh(isoc, nadh, nad, adp_m, ca_m, h_m, KM_ADP, KM_CA,
	          KM_NAD, KI_NADH, KM_ISOC, NI_ISOC, VMAX, KH1, KH2)
	fa = _mm(KM_ADP, adp_m) * _mm(KM_CA, ca_m)
	f_isoc = (KM_ISOC / isoc)^NI_ISOC
	f_nad = _mm_reci(KI_NADH, nadh) * (KM_NAD / nad)
	f_h = h_m / KH1 + KH2 / h_m
	vIDH = VMAX / (f_h + (1 + fa * f_isoc) *  (1 + f_nad))
end
# With correction from dissociation
function vidh(isoc, nadh, nad, adp_m, ca_m, h_m, mg_m, KM_ADP, KM_CA,
	          KM_NAD, KI_NADH, KM_ISOC, NI_ISOC, VMAX, KH1, KH2)
	(adp3, _, _, _) = breakdown_adp(adp_m, h_m, mg_m)
	vidh(isoc, nadh, nad, adp_m, ca_m, h_m, KM_ADP, KM_CA,
		 KM_NAD, KI_NADH, KM_ISOC, NI_ISOC, VMAX, KH1, KH2)
end

function vidh(isoc, nadh, nad, adp_m, ca_m, h_m, mg_m, pIDH::IDHParams)
	@unpack KM_ADP, KM_CA, KM_NAD, KI_NADH, KM_ISOC, NI_ISOC, VMAX, KH1, KH2 = pIDH
	vidh(isoc, nadh, nad, adp_m, ca_m, h_m, mg_m, KM_ADP, KM_CA,
		 KM_NAD, KI_NADH, KM_ISOC, NI_ISOC, VMAX, KH1, KH2)
end

# KGDH (alpha-ketoglutarate dehydrogenase)
@with_kw struct KGDHParams
    ET = 0.5
    KCAT = 0.0179
    VMAX = KCAT * ET
    KM_AKG = 30
    KM_NAD = 38.7
    KH1 = 4E-5
    KH2 = 7E-5
    KM_MG = 0.0308
    KM_CA = 1.5E-4
    NI_AKG = 1.2
end

# KGDH (alpha-ketoglutarate dehydrogenase) rates
function vkgdh(akg, nadh, nad, ca_m, mg_m, KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, VMAX, f_h=1)
	f_a = _mm(KM_MG, mg_m) * _mm(KM_CA, ca_m)
    f_akg = (KM_AKG / akg) ^ NI_AKG
    f_nad = KM_NAD / nad
    vKGDH = VMAX / (f_h + f_a * (f_akg + f_nad))
end

# KGDH (alpha-ketoglutarate dehydrogenase) rates, with proton activation
function vkgdh(akg, nadh, nad, h_m, ca_m, mg_m, KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, KH1, KH2, VMAX)
	f_h = 1 + h_m / KH1 + KH2 / h_m
	vkgdh(akg, nadh, nad, ca_m, mg_m, KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, VMAX, f_h)
end

function vkgdh(akg, nadh, nad, h_m, ca_m, mg_m, pKGDH::KGDHParams)
	@unpack KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, KH1, KH2, VMAX = pKGDH
	vkgdh(akg, nadh, nad, h_m, ca_m, mg_m, KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, KH1, KH2, VMAX)
end

function vkgdh(akg, nadh, nad, ca_m, mg_m, pKGDH::KGDHParams)
	@unpack KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, VMAX = pKGDH
	vkgdh(akg, nadh, nad, ca_m, mg_m, KM_MG, KM_CA, KM_AKG, NI_AKG, KM_NAD, VMAX)
end


# SL (Succinyl-coA lyase)
@with_kw struct SLParams
	COA = 0.02
    KF = 0.0284
    KEQ = 3.115
end

# SL (succinyl-CoA lyase) rate
vsl(scoa, suc, adp_m, atp_m, COA, KF, KEQ) = KF * (scoa * adp_m - suc * atp_m * COA / KEQ)
vsl(scoa, suc, adp_m, atp_m, pSL::SLParams) = vsl(scoa, suc, adp_m, atp_m, pSL.COA, pSL.KF, pSL.KEQ)

# With corrections from dissociation
function vsl(scoa, suc, adp_m, atp_m, h_m, pi_m, mg_m, COA, KF, KEQ)
	(atp4, hatp, _, poly_atp) = breakdown_atp(atp_m, h_m, mg_m)
    (_, _, _, poly_adp) = breakdown_adp(adp_m, h_m, mg_m)
	ke_app = KEQ * (poly_atp * _suc_poly(h_m)) / (poly_adp * _pi_poly(h_m))
	vSL = vsl(scoa, suc, adp_m * pi_m, atp4 + hatp, COA, KF, ke_app)
end

function vsl(scoa, suc, adp_m, atp_m, h_m, pi_m, mg_m, pSL::SLParams)
	@unpack COA, KF, KEQ = pSL
	vSL = vsl(scoa, suc, adp_m, atp_m, h_m, pi_m, mg_m, COA, KF, KEQ)
 end

# FH (Fumarate hydrase) parameters
@with_kw struct FHParams
    KF = 8.3E-3
    KEQ = 1.0
end

# FH (fumarate hydratase) rates
vfh(fum, mal, KF, KEQ) = KF * (fum - mal / KEQ)
vfh(fum, mal, pFH::FHParams) = vfh(fum, mal, pFH.KF, pFH.KEQ)

# MDH (Malate dehydrogenase)
@with_kw struct MDHParams
    KH1 = 1.131E-5
    KH2 = 26.7
    KH3 = 6.68E-9
    KH4 = 5.62E-6
    K_OFFSET = 3.99E-2
    KCAT = 0.1259
    ET = 0.154
    VMAX = KCAT * ET
    KM_MAL = 1.493
    KI_OAA = 3.1E-3
    KM_NAD = 0.2244
end

# MDH (malate dehydrogenase) rates
function vmdh(mal, oaa, nadh, nad, h_m, K_OFFSET, KH1, KH2, KH3, KH4, KM_MAL, KI_OAA, KM_NAD, VMAX)
	f_ha = K_OFFSET + _mm(KH1 * _mm(KH2, h_m), h_m)
	f_hi = _mm(h_m * _mm(h_m, KH4), KH3)^2
	f_oaa = _mm(KI_OAA, oaa)
    f_mal = _mm(mal * f_oaa, KM_MAL)
    f_nad = _mm(nad, KM_NAD)
    vMDH = VMAX * f_ha * f_hi * f_nad * f_mal
end

function vmdh(mal, oaa, nadh, nad, h_m, pMDH::MDHParams)
	@unpack K_OFFSET, KH1, KH2, KH3, KH4, KM_MAL, KI_OAA, KM_NAD, VMAX = pMDH
	vmdh(mal, oaa, nadh, nad, h_m, K_OFFSET, KH1, KH2, KH3, KH4, KM_MAL, KI_OAA, KM_NAD, VMAX)
end

# AAT (Transaminase)
@with_kw struct AATParams
	KF = 0.0217
    KEQ = 6.6
    KASP = 1.5E-6
    GLU = 30.0
end

# AAT (alanine aminotransferase) rates. scalar version
vaat(akg, oaa, KF, KEQ, KASP, GLU) = KF * oaa * GLU * _mm(KASP * KEQ, akg * KF)
function vaat(akg, oaa, pAAT::AATParams)
	@unpack KF, KEQ, KASP, GLU = pAAT
	vAAT = vaat(akg, oaa, KF, KEQ, KASP, GLU)
end

# TCA cyle rates, with corrections from dissociation, scalar version
function tca_system(isoc, akg, scoa, suc, fum, mal, oaa, nadh, adp_m, pi_m, h_m, ca_m, vSDH,
				   conc::Concentrations, pSL::SLParams, pCS::CSParams, pACO::ACOParams, pIDH::IDHParams,
				   pKGDH::KGDHParams, pFH::FHParams, pMDH::MDHParams, pAAT::AATParams)
	@unpack AXP_M_T, NAD_T, mg_m = conc
	nad = NAD_T - nadh
	atp_m = AXP_M_T - adp_m
	vCS  = vcs(oaa, pCS)
	vACO = vaco(isoc, akg, scoa, suc, fum, mal, oaa, pACO)
	vIDH = vidh(isoc, nadh, nad, adp_m, ca_m, h_m, mg_m, pIDH)
	vKGDH = vkgdh(akg, nadh, nad, h_m, ca_m, mg_m, pKGDH)
	vSL = vsl(scoa, suc, adp_m, atp_m, h_m, pi_m, mg_m, pSL)
	vFH = vfh(fum, mal, pFH)
	vMDH = vmdh(mal, oaa, nadh, nad, h_m, pMDH)
	vAAT = vaat(akg, oaa, pAAT)

	d_nadh_tca = vIDH + vKGDH + vMDH
	d_isoc = vACO - vIDH
	d_akg = vIDH - vKGDH + vAAT
	d_scoa = vKGDH - vSL
	d_suc = vSL - vSDH
	d_fum = vSDH - vFH
	d_mal = vFH - vMDH
	d_oaa = vMDH - vCS - vAAT
	return d_isoc, d_akg, d_scoa, d_suc, d_fum, d_mal, d_oaa, d_nadh_tca, vSL
end
