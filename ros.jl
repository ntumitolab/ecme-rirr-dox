#=
ROS scavenging (detoxification) system
=#
using Parameters
include("common.jl")
include("concentrations.jl")

# SOD(superoxide dismutase) parameters
@with_kw struct SODParams
    K1 = 1200.0  # 2nd order rate constant of SOD (1/mMms)
    K3 = 24.0  # 2nd order rate constant of SOD (1/mMms)
    K5 = 2.4E-4  # 1st order rate constant of SOD (1/ms)
    KI_H2O2 = 0.5  # Inhibition constant of H2O2 (mM)
    ET = 0.0003  # SOD concentration (mM)
    VMAX = ET * K5
end

# Rate of superoxide dismutase (mM/ms),
function vsod(sox, h2o2, K1, K3, K5, KI_H2O2, VMAX)
	f_h2o2 = K3 * _mm_reci(KI_H2O2, h2o2)
	f_sox  = K1 * sox
	denom = K5 * (2 * K1 + f_h2o2) + f_h2o2 * f_sox
	vSOD = 2 * VMAX * f_sox * (K1 + f_h2o2) / denom
end

function vsod(sox, h2o2, pSOD::SODParams)
	@unpack K1, K3, K5, KI_H2O2, VMAX = pSOD
	vsod(sox, h2o2, K1, K3, K5, KI_H2O2, VMAX)
end

# GPX (Glutathione peroxidase) parameters
@with_kw struct GPXParams
    𝚽1 = 5E-5  # Rate constant of GPX
    𝚽2 = 27.0  # Rate constant of GPX
    ET = 9.77E-5  # GPX concentration (mM)
end

@with_kw struct TPXParams
	𝚽1 = 3.83  # Rate constant of TPX
    𝚽2 = 1.85  # Rate constant of TPX
    ET = 3E-3  # TPX concentration (mM)
end

#=
Reaction rate of glutathione peroxidase (GPX) and thioredoxin peroxidase (TPX)
H2O2 + 2GSH -> 2H2O + GSSH
H2O2 + TRXSH2 -> 2H2O + TRXSS
=#
vgpx(h2o2, gsh, ET, 𝚽1, 𝚽2) = ET * h2o2 * gsh / (𝚽1 * gsh + 𝚽2 * h2o2)
vgpx(h2o2, gsh, pGPX::GPXParams) = vgpx(h2o2, gsh, pGPX.ET, pGPX.𝚽1, pGPX.𝚽2)
vtpx(h2o2, trxsh2, ET, 𝚽1, 𝚽2) = vgpx(h2o2, trxsh2, ET, 𝚽1, 𝚽2)
vtpx(h2o2, trxsh2, pTPX::TPXParams) = vtpx(h2o2, trxsh2, pTPX.ET, pTPX.𝚽1, pTPX.𝚽2)

# GR (glutathion reductace) parameters
@with_kw struct GRParams
    K1 = 2.5E-3  # Catalytic constant of GR
    ET = 2.3E-3  # Mitochondiral concentration of GR (mM)
    VMAX = K1 * ET
    KM_NADPH = 6.54E-2  # Michaelis constant for NADPH of GR (mM)
    KM_GSSG = 20.6E-3  # Michaelis constant for oxidized GSH of GR (mM)
    GSH_T = 1.0  # Glutathione pool (mM)
end

# Rate of glutathione reductase (GR) (mM/ms)
# GSSG + NADPH + H+ -> 2GSH + NADP+
vgr(gssg, nadph, VMAX, KM_GSSG, KM_NADPH) = VMAX * _mm(nadph, KM_NADPH) * _mm(gssg, KM_GSSG)
vgr(gsh, nadph, VMAX, KM_GSSG, KM_NADPH, GSH_T) = vgr(0.5 * (GSH_T - gsh), nadph, VMAX, KM_GSSG, KM_NADPH)
function vgr(gsh, nadph, pGR::GRParams)
	@unpack VMAX, KM_GSSG, KM_NADPH, GSH_T = pGR
	vGR = vgr(gsh, nadph, VMAX, KM_GSSG, KM_NADPH, GSH_T)
end

# Thioredoxin reductase (TR) parameters
@with_kw struct TRParams
    K1 = 22.75E-3  # Catalytic constant of GR
    ET = 3.5E-4  # Mitochondiral concentration of GR (mM)
    VMAX = K1 * ET
    KM_NADPH = 0.065  # Michaelis constant for NADPH of GR (mM)
    KM_TRXSS = 0.035  # Michaelis constant for oxidized thioredoxin of TR (mM)
    TRX_T = 0.025  # Thioredoxin pool(mM)
end

# Thioredoxin reductase (TR) (mM/ms)
# TRXSS + NADPH -> TRXSH2 + NADP
vtr(trxsh2, nadph, VMAX, KM_TRXSS, KM_NADPH, TRX_T) = vgr(TRX_T - trxsh2, nadph, VMAX, KM_TRXSS, KM_NADPH)
function vtr(trxsh2, nadph, pTR::TRParams)
	@unpack VMAX, KM_TRXSS, KM_NADPH, TRX_T = pTR
	vTR = vtr(trxsh2, nadph, VMAX, KM_TRXSS, KM_NADPH, TRX_T)
end

#=
Rate of glutathione reductase (GR) (mM/ms)
GSSH + NADPH -> 2GSH + NADP
as well as for thioredoxin reductase (TR) (mM/ms)
TRXSS + NADPH -> TRXSH2 + NADP
=#

# Rate of H2O2 diffusion
vdiff_h2o2(h2o2_m, h2o2_i, C_DIFF) = _vdiff(h2o2_m, h2o2_i, C_DIFF)

# Catalase (CAT)
@with_kw struct CATParams
    K1 = 17.0  # Rate constant of CAT (1/mMms)
    ET = 1E-6  # Total pool of CAT (mM)
    VMAX = ET * K1
    FR = 0.05  # H2O2 inhibition factor of CAT
end

# Rate of intracellular catalase (CAT), scalar version
v_cat(h2o2_i, VMAX, FR) = 2 * VMAX * h2o2_i * exp(-FR * h2o2_i)
v_cat(h2o2_i, pCAT::CATParams) = v_cat(h2o2_i, pCAT.VMAX, pCAT.FR)

# IMAC (Inner mitochondrial anion channel) from Cortassa et al. (2004)
@with_kw struct IMACParams
	A = 1E-3  # Basal IMAC conductance
    B = 1E4  # Activation factor by cytoplasmic superoxide
    KCC_SOX = 0.01  # Michaelis constant for cytoplasmic superoxide of IMAC
    GL = 3.5E-8  # Leak conductance of IMAC
    G_MAX = 3.9085E-6  # Maximal conductance of IMAC
    κ = 0.07  # Steepness factor
    DPSI_OFFSET = 4.0  # Potential at half saturation
    J = 0.1  # Fraction of IMAC conductance
end

#=
Rate of mitochondrial inner membrane anion channel (IMAC) (mM/ms),
as well as superoxide diffusion from mitochondira to cytosol (VtrROS)
=#
# Conductance of IMAC (scalar version)
function _g_imac(sox_i, dpsi, A, B, KCC_SOX, GL, G_MAX, κ, DPSI_OFFSET)
	basal = A + B * _mm(sox_i, KCC_SOX)
	voltFactor = GL + G_MAX / (1 + exp(κ * (DPSI_OFFSET + dpsi)))
	gIMAC = voltFactor * basal
end

# scalar version of IMAC and superoxide pass-through (vTrROS)
function vimac_vtrros(dpsi, sox_m, sox_i, A, B, KCC_SOX, GL, G_MAX, κ, DPSI_OFFSET, J)
    gimac = _g_imac(sox_i, dpsi, A, B, KCC_SOX, GL, G_MAX, κ, DPSI_OFFSET)
	# Workaround negative ROS problem
	sox_i = max(1e-120, sox_i)
	sox_m = max(1e-120, sox_m)

    vTrROS = J * gimac * (dpsi + _esox(sox_i, sox_m))
    vIMAC = gimac * dpsi
    return vIMAC, vTrROS
end

function vimac_vtrros(dpsi, sox_m, sox_i, pIMAC::IMACParams)
	@unpack A, B, KCC_SOX, GL, G_MAX, κ, DPSI_OFFSET, J = pIMAC
	vimac_vtrros(dpsi, sox_m, sox_i, A, B, KCC_SOX, GL, G_MAX, κ, DPSI_OFFSET, J)
end

# nadph_m-producing isocitrate dehydrogenase (IDH2) parameters
@with_kw struct IDH2Params
    KM_H = 0.5  # Dissociation constant for H+ of IDH2
    KM_ISOC = 0.045  # Dissociation constant for isocitrate of IDH2
    KM_NADP = 0.046  # Michealis constant for NADP of IDH2
    KI_NADP = 2E-6  # Inhibition constant for NADP of IDH2
    KM_NADPH = 1.2E-2  # Michealis constant for nadph_m of IDH2
    KM_AKG = 0.08  # Michealis constant for AKG of IDH2
    VF = 8.72E-2  # Max forward rate of IDH2
    VB = 5.45E-3  # Max backward rate of IDH2
end

# Rate of isocitrate dehydrogenase II (producing nadph_m). Scalar parameters.
function vidh2(nadph_m, nadp_m, isoc, akg, h_m, KM_H, KM_ISOC, KM_NADP, KI_NADP, KM_AKG, KM_NADPH, VF, VB)
	f_h = _mm(KM_H, h_m)
    f_isoc = isoc / KM_ISOC
    f_nadp = nadp_m / (KM_NADP * _mm(nadp_m, KI_NADP))
    f_akg = akg / KM_AKG
    f_nadph = nadph_m / KM_NADPH
	denom = (1 + f_isoc + f_akg) * (1 + f_nadp + f_nadph)
	num = (VF * f_isoc * f_nadp - VB * f_akg * f_nadph)
	vIDH2 = num * f_h / denom
end

function vidh2(nadph_m, nadp_m, isoc, akg, h_m, pIDH2::IDH2Params)
	@unpack KM_H, KM_ISOC, KM_NADP, KI_NADP, KM_AKG, KM_NADPH, VF, VB = pIDH2
	vidh2(nadph_m, nadp_m, isoc, akg, h_m, KM_H, KM_ISOC, KM_NADP, KI_NADP, KM_AKG, KM_NADPH, VF, VB)
end

# transhydrogenase (THD) parameters
@with_kw struct THDParams
    KM_NADH = 0.01  # Michealis constant for NADH of THD
    KM_NAD = 0.125  # Michealis constant for NAD of THD
    KM_NADP = 0.02  # Michealis constant for NADP of THD
    KM_NADPH = 0.02  # Michealis constant for nadph_m of THD
    ET = 1E-3  # THD concnetration (mM)
    KF = 1.17474  # Max forward rate of THD
    VF = ET * KF
    KB = 17.2756  # Max backward rate of THD
    VB = ET * KB
    D = 0.5  # Voltage assymetry factor of THD
    X = 0.1  # Voltage dependence factor of THD
end

# Rate of transhydrogenase (producing nadph_m). Scalar version
function vthd(nadh, nad, nadph_m, nadp_m, ΔμH, X, D, KM_NAD, KM_NADH, KM_NADP, KM_NADPH, VF, VB)
	# Corrected an error from both the paper and code from Gauthier et al.
	vNAD = _ra(X * (1.0 - D) * ΔμH)
	fNAD = nad / KM_NAD
	fvNAD = fNAD * vNAD
	fNADH = nadh / KM_NADH
	fNADP = nadp_m / KM_NADP
	vNADP = _ra(X * D * ΔμH)
	fvNADP = vNADP * fNADP
	fNADPH = nadph_m / KM_NADPH
	denom = one(fNAD) + fNAD + fNADH + fNADP + fNADPH + (fNADPH + fvNAD) * (fNADPH + fvNADP)
	num = VF * fNADH * fvNADP - VB * fvNAD * fNADPH
	vTHD = num / denom
end

function vthd(nadh, nad, nadph_m, nadp_m, ΔμH, pTHD::THDParams)
	@unpack X, D, KM_NAD, KM_NADH, KM_NADP, KM_NADPH, VF, VB = pTHD
	vthd(nadh, nad, nadph_m, nadp_m, ΔμH, X, D, KM_NAD, KM_NADH, KM_NADP, KM_NADPH, VF, VB)
end

function ros_system(nadh, nadph_m, sox_m, sox_i, h2o2_m, h2o2_i, gsh_m, gsh_i,
	trxsh2_m, trxsh2_i, dpsi, ΔμH, isoc, akg, h_m, vROS, conc::Concentrations,
	pSODi::SODParams, pSODm::SODParams, pGPXm::GPXParams, pGPXi::GPXParams, pTPXm::TPXParams,
	pTPXi::TPXParams, pGRm::GRParams, pGRi::GRParams, pTRm::TRParams, pTRi::TRParams, pCAT::CATParams,
	pIMAC::IMACParams, pIDH2::IDH2Params, pTHD::THDParams, C_DIFF_H2O2, V_MITO_V_MYO)
	@unpack NAD_T, NADP_M_T, nadph_i = conc
	nad = NAD_T - nadh
	nadp_m = NADP_M_T - nadph_m
	vIDH2 = vidh2(nadph_m, nadp_m, isoc, akg, h_m, pIDH2)
	vTHD = vthd(nadh, nad, nadph_m, nadp_m, ΔμH, pTHD)
	vIMAC, vTrROS = vimac_vtrros(dpsi, sox_m, sox_i, pIMAC)
	vCAT = v_cat(h2o2_i, pCAT)
	vDiffH2O2 = vdiff_h2o2(h2o2_m, h2o2_i, C_DIFF_H2O2)
	vTR_i = vtr(trxsh2_i, nadph_i, pTRi)
	vTR_m = vtr(trxsh2_m, nadph_m, pTRm)
	vGR_i = vgr(gsh_i, nadph_i, pGRi)
	vGR_m = vgr(gsh_m, nadph_m, pGRm)
	vTPX_i = vtpx(h2o2_i, trxsh2_i, pTPXi)
	vTPX_m = vtpx(h2o2_m, trxsh2_m, pTPXm)
	vGPX_i = vgpx(h2o2_i, gsh_i, pGPXi)
	vGPX_m = vgpx(h2o2_m, gsh_m, pGPXm)
	vSOD_i = vsod(sox_i, h2o2_i, pSODi)
	vSOD_m = vsod(sox_m, h2o2_m, pSODm)

	d_nadph_m = vIDH2 + vTHD - vGR_m - vTR_m
	d_sox_m = vROS - vTrROS - vSOD_m
	d_sox_i = V_MITO_V_MYO * vTrROS - vSOD_i
	d_h2o2_m = 0.5 * vSOD_m - vDiffH2O2 - vTPX_m - vGPX_m
	d_h2o2_i = 0.5 * vSOD_i + V_MITO_V_MYO * vDiffH2O2 - vTPX_i - vGPX_i - vCAT
	d_gsh_m  = vGR_m - vGPX_m
	d_gsh_i  = vGR_i - vGPX_i
	d_trxsh2_m = vTR_m - vTPX_m
	d_trxsh2_i = vTR_i - vTPX_i
	return ( d_nadph_m, d_sox_m, d_sox_i, d_h2o2_m, d_h2o2_i, d_gsh_m, d_gsh_i,
		d_trxsh2_m, d_trxsh2_i, vIDH2, vTHD, vIMAC)
end
