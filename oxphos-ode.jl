#=
Detailed Oxidative phosphorylation model by Gauthier et al. (2013)
With ROS generation
Default parameter values from Kembro et al. and Gauthier et al.
Some are adjusted by An-chi Wei to prevent complex number
=#
using Parameters
include("common.jl")
include("concentrations.jl")

const E_FMN = -0.375  # FMN redox potential (mV)
const E_SOX = -0.15  # Superoxide redox potential (mV)
const E_ROS = _ra(E_SOX - E_FMN)  # Voltage factor between Superoxide and FMN

# Doxorubicin effects
@with_kw struct DOXParams
    DOX = 0.0  # Doxorubicin (DOX) concentration (mM)
    NI = 3  # Hills coefficient of DOX inhibition on OXPHOS complex
    KI_COX_C1 = 400e-3  # DOX inhibition concentration (IC50) on complex I
    KI_COX_C2 = 2000e-3 # DOX inhibition concentration (IC50) on complex II
    KI_COX_C3 = 185e-3  # DOX inhibition concentration (IC50) on complex III
    KI_COX_C4 = 165e-3  # DOX inhibition concentration (IC50) on complex IV
    K_RC = 1E3 / 15  # DOX redox cycling constant
    C1_INHIB = _hills(KI_COX_C1, DOX, NI)  # complex I inhibitory scale
    C2_INHIB = _hills(KI_COX_C2, DOX, NI) # complex II inhibitory scale
    C3_INHIB = _hills(KI_COX_C3, DOX, NI)  # complex III inhibitory scale
    C4_INHIB = _hills(KI_COX_C4, DOX, NI)  # complex IV inhibitory scale
    E_LEAK_C1 = 1 + K_RC * DOX  # Electron leak scaling factor from complex I
end

# OXPHOS (complex I)
@with_kw struct C1Params
    SCALE = 1E3  # Convert 1/s to 1/ms
	ρC1 = 24.9302  # Concentration of complex I (mM), from Gauthier et al. (2013)
    dpsi_B = 50.0  # Phase boundary potential (mV)
    K12 = 6.3396E11 / SCALE
    K21 = 5 / SCALE
    K56 = 100 / SCALE
    K65 = 2.5119E13 / SCALE
    K61 = 3.329E7 / SCALE
    K16 = 432.7642 / SCALE
    K23 = 3886.7 / SCALE
    K32 = 9.1295E6 / SCALE
    K34 = 639.1364 / SCALE
    K43 = 3.2882 / SCALE
    K47 = 1.5962E7 / SCALE
    K74 = 65.2227 / SCALE
    K75 = 2.4615E4 / SCALE
    K57 = 1.1667E3 / SCALE
    K42 = 6.0318 / SCALE
end

function c1_rates(nadh, nad, Q_n, QH2_n, dpsi, h_m, sox_m,
                  dpsi_B, K12, K21, K65, K56, K61, K16, K23, K32,
                  K34, K43, K47, K74, K57, K75, K42, ρC1,
                  MT_PROT=1, E_LEAK_C1=1, C1_INHIB=1, h_i=1E-4, O2=6E-3)
    C1_CONC = ρC1 * MT_PROT
	# A hack to wrokaround negative concentration
	nadh = max(0.0, nadh)
	nad = max(0.0, nad)
	Q_n = max(0.0, Q_n)
	QH2_n = max(0.0, QH2_n)
    a = _ra(dpsi - dpsi_B)
    a12 = K12 * h_m^2
    a21 = K21
    a65 = K65 * h_i^2
    a56 = K56
    a61 = K61 / a
    a16 = K16 * a
    a23 = K23 * sqrt(nadh)
    a32 = K32
    a34 = K34
    a43 = K43 * sqrt(nad)
    a47 = C1_INHIB * K47 * sqrt(Q_n * h_m)
    a74 = K74
    a57 = C1_INHIB * K57 * sqrt(QH2_n)
    a75 = K75
    a42 = K42 * E_LEAK_C1 * O2
    a24 = K42 * O2 * E_ROS * sox_m
	a52 = a25 = 0

    # Fraction of each state in Complex I from the matlab code, derived from KA pattern
	e1 = a21*a32*a42*a56*a61*a74 + a21*a32*a42*a56*a61*a75 + a21*a32*a42*a57*a61*a74 + a21*a32*a42*a57*a65*a74 + a21*a32*a43*a56*a61*a74 + a21*a32*a43*a56*a61*a75 + a21*a32*a43*a57*a61*a74 + a21*a32*a43*a57*a65*a74 + a21*a32*a47*a56*a61*a75 + a21*a34*a42*a56*a61*a74 + a21*a34*a42*a56*a61*a75 + a21*a34*a42*a57*a61*a74 + a21*a34*a42*a57*a65*a74 + a21*a34*a47*a56*a61*a75 + a23*a34*a47*a56*a61*a75 + a24*a32*a47*a56*a61*a75 + a24*a34*a47*a56*a61*a75
	e2 = a12*a32*a42*a56*a61*a74 + a12*a32*a42*a56*a61*a75 + a12*a32*a42*a57*a61*a74 + a12*a32*a42*a57*a65*a74 + a12*a32*a43*a56*a61*a74 + a12*a32*a43*a56*a61*a75 + a12*a32*a43*a57*a61*a74 + a12*a32*a43*a57*a65*a74 + a12*a32*a47*a56*a61*a75 + a12*a34*a42*a56*a61*a74 + a12*a34*a42*a56*a61*a75 + a12*a34*a42*a57*a61*a74 + a12*a34*a42*a57*a65*a74 + a12*a34*a47*a56*a61*a75 + a16*a32*a42*a57*a65*a74 + a16*a32*a43*a57*a65*a74 + a16*a34*a42*a57*a65*a74
	e3 = a12*a23*a42*a56*a61*a74 + a12*a23*a42*a56*a61*a75 + a12*a23*a42*a57*a61*a74 + a12*a23*a42*a57*a65*a74 + a12*a23*a43*a56*a61*a74 + a12*a23*a43*a56*a61*a75 + a12*a23*a43*a57*a61*a74 + a12*a23*a43*a57*a65*a74 + a12*a23*a47*a56*a61*a75 + a12*a24*a43*a56*a61*a74 + a12*a24*a43*a56*a61*a75 + a12*a24*a43*a57*a61*a74 + a12*a24*a43*a57*a65*a74 + a16*a21*a43*a57*a65*a74 + a16*a23*a42*a57*a65*a74 + a16*a23*a43*a57*a65*a74 + a16*a24*a43*a57*a65*a74
	e4 = a12*a23*a34*a56*a61*a74 + a12*a23*a34*a56*a61*a75 + a12*a23*a34*a57*a61*a74 + a12*a23*a34*a57*a65*a74 + a12*a24*a32*a56*a61*a74 + a12*a24*a32*a56*a61*a75 + a12*a24*a32*a57*a61*a74 + a12*a24*a32*a57*a65*a74 + a12*a24*a34*a56*a61*a74 + a12*a24*a34*a56*a61*a75 + a12*a24*a34*a57*a61*a74 + a12*a24*a34*a57*a65*a74 + a16*a21*a32*a57*a65*a74 + a16*a21*a34*a57*a65*a74 + a16*a23*a34*a57*a65*a74 + a16*a24*a32*a57*a65*a74 + a16*a24*a34*a57*a65*a74
	e5 = a12*a23*a34*a47*a61*a75 + a12*a23*a34*a47*a65*a75 + a12*a24*a32*a47*a61*a75 + a12*a24*a32*a47*a65*a75 + a12*a24*a34*a47*a61*a75 + a12*a24*a34*a47*a65*a75 + a16*a21*a32*a42*a65*a74 + a16*a21*a32*a42*a65*a75 + a16*a21*a32*a43*a65*a74 + a16*a21*a32*a43*a65*a75 + a16*a21*a32*a47*a65*a75 + a16*a21*a34*a42*a65*a74 + a16*a21*a34*a42*a65*a75 + a16*a21*a34*a47*a65*a75 + a16*a23*a34*a47*a65*a75 + a16*a24*a32*a47*a65*a75 + a16*a24*a34*a47*a65*a75
	e6 = a12*a23*a34*a47*a56*a75 + a12*a24*a32*a47*a56*a75 + a12*a24*a34*a47*a56*a75 + a16*a21*a32*a42*a56*a74 + a16*a21*a32*a42*a56*a75 + a16*a21*a32*a42*a57*a74 + a16*a21*a32*a43*a56*a74 + a16*a21*a32*a43*a56*a75 + a16*a21*a32*a43*a57*a74 + a16*a21*a32*a47*a56*a75 + a16*a21*a34*a42*a56*a74 + a16*a21*a34*a42*a56*a75 + a16*a21*a34*a42*a57*a74 + a16*a21*a34*a47*a56*a75 + a16*a23*a34*a47*a56*a75 + a16*a24*a32*a47*a56*a75 + a16*a24*a34*a47*a56*a75
	e7 = a12*a23*a34*a47*a56*a61 + a12*a23*a34*a47*a57*a61 + a12*a23*a34*a47*a57*a65 + a12*a24*a32*a47*a56*a61 + a12*a24*a32*a47*a57*a61 + a12*a24*a32*a47*a57*a65 + a12*a24*a34*a47*a56*a61 + a12*a24*a34*a47*a57*a61 + a12*a24*a34*a47*a57*a65 + a16*a21*a32*a42*a57*a65 + a16*a21*a32*a43*a57*a65 + a16*a21*a32*a47*a57*a65 + a16*a21*a34*a42*a57*a65 + a16*a21*a34*a47*a57*a65 + a16*a23*a34*a47*a57*a65 + a16*a24*a32*a47*a57*a65 + a16*a24*a34*a47*a57*a65
    denom = e1 + e2 + e3 + e4 + e5 + e6 + e7
    f2 = e2 / denom
    f4 = e4 / denom
    f7 = e7 / denom
	# Electrons flow through C1
    jResC1 = C1_CONC * (f4 * a47 - f7 * a74)
    # ROS produced in C1
    vROSC1 = C1_CONC * (f4 * a42 - f2 * a24)
    # Q reduction rate from C1
    vc1 = 0.5 * jResC1
    return vc1, vROSC1
end

function c1_rates(nadh, Q_n, QH2_n, dpsi, h_m, sox_m, conc::Concentrations, pC1::C1Params, pDOX::DOXParams, MT_PROT=1)
	@unpack  (dpsi_B, K12, K21, K65, K56, K61, K16, K23, K32,
	K34, K43, K47, K74, K57, K75, K42, ρC1) = pC1
	@unpack NAD_T, h_i, O2 = conc
	@unpack E_LEAK_C1, C1_INHIB = pDOX
	nad = _nad(nadh, conc)
	return c1_rates(nadh, nad, Q_n, QH2_n, dpsi, h_m, sox_m,
	         dpsi_B, K12, K21, K65, K56, K61, K16, K23, K32,
	         K34, K43, K47, K74, K57, K75, K42, ρC1,
	         MT_PROT, E_LEAK_C1, C1_INHIB, h_i, O2)
end

# OXPHOS (complex II)
@with_kw struct C2Params
    VMAX = 250 / 60E3  # Maximal rate of SDH = complex II (mM/ms)
    KM_Q = 0.6  # MM constant for CoQ
    KI_OAA = 0.15  # MM constant for OAA
	FAC_SUC = 0.085
	KI_FUM = 1.3
	KM_SUC = 0.03
	ET = 0.5
	KCAT = 3E-3
end

#=
Rate of complex II, aka succinate dehydrogenase (SDH)
SUC + Q -> (via FAD/FADH2) -> FUM + QH2
with a magical SUC_SCALE term from the matlab code
adapted from Demin et al. (2000)
=#
function vsdh(Q_n, QH2_n, suc, fum, oaa, KI_OAA, VMAX, KM_Q, FAC_SUC, C2_INHIB=1)
    f_q = _mm(_mm(Q_n, QH2_n), KM_Q)
    f_oaa = _mm(KI_OAA, oaa)
    f_suc = FAC_SUC * sqrt(suc / fum)
    vSDH = VMAX * C2_INHIB * f_suc * f_oaa * f_q
end

function vsdh(Q_n, QH2_n, suc, fum, oaa, pC2::C2Params, C2_INHIB=1)
	@unpack KI_OAA, VMAX, KM_Q, FAC_SUC = pC2
	vSDH = vsdh(Q_n, QH2_n, suc, fum, oaa, KI_OAA, VMAX, KM_Q, FAC_SUC, C2_INHIB)
end

# SDH rate combines Cortassa et al. (2003) and Demin et al. (2000)
function vsdh_mm(Q_n, QH2_n, suc, fum, oaa, pC2::C2Params, C2_INHIB=1)
	@unpack KI_OAA, VMAX, KM_Q, KI_FUM, KM_SUC = pC2
	f_q = _mm(_mm(Q_n, QH2_n), KM_Q)
	f_oaa = _mm(KI_OAA, oaa)
	f_fum = _mm(KI_FUM, fum)
	f_suc = _mm(suc * f_oaa * f_fum, KM_SUC)
	vSDH = C2_INHIB * VMAX * f_suc * f_q
end

# OXPHOS (complex IV)
@with_kw struct C4Params
	# Concentration of complex IV (mM)
    ρC4 = 0.325
    SCALE = 60E3  # 1/min to 1/ms
    δ₅ = 0.5
    K34 = 1.7667E28 / SCALE
    K43 = 1.7402 / SCALE
    K35 = 45000 / SCALE
    K36 = 2.8955E25 / SCALE
    K63 = 2.8955E10 / SCALE
    K37 = 1.7542E12 / SCALE
    K73 = 1.7542E4 / SCALE
end

# C4 rates (scalar version)
function c4_rates(cytc_ox, dpsi, h_m, conc::Concentrations, pC4::C4Params, C4_INHIB=1, MT_PROT=1)
	@unpack K34, K73, K43, K35, δ₅, K36, K37, K63, ρC4 = pC4
	@unpack h_i, O2 = conc
    C4_CONC = ρC4 * MT_PROT

    aδ = _ra(-δ₅ * dpsi)
    a1mδ = _ra((1 - δ₅) * dpsi)
	f_hm = h_m * aδ
	f_hi = h_i * a1mδ
	f_cr = cytc_rd = C4_CONC - cytc_ox
	f_co = cytc_ox * a1mδ
    # Transition rates
    a12 = K34 * f_cr^3 * f_hm^4  # a12 = K34 * exp(-δ₅ * 4 * vfrt) * cytc_rd^3 * h_m^4
    a14 = K73 * f_hi  # K73 * exp((1 - δ₅) * F_RT * dpsi) * h_i
    a21 = K43 * f_co^3 * f_hi  # K43 * exp((1 - δ₅) * 4 * vfrt) * cytc_ox^3 * h_i
    a23 = K35 * O2 * C4_INHIB
    a34 = K36 * f_cr * f_hm^3  # K36 * exp(-δ₅ * 3 * vfrt) * cytc_rd * h_m^3
    a41 = K37 * f_hm  # K37 * exp(-δ₅ * vfrt) * h_m
    a43 = K63 * f_co * f_hi^2  # K63 * exp((1 - δ₅) * 3 * vfrt) * cytc_ox * h_i^2

    # Weight of each state (from KA pattern)
    e1 = a21 * a41 * a34 + a41 * a34 * a23
    e2 = a12 * a41 * a34
    e3 = a23 * a12 * a41 + a43 * a14 * a21 + a23 * a43 * a12 + a23 * a43 * a14
    e4 = a14 * a34 * a21 + a34 * a23 * a12 + a34 * a23 * a14
    den = e1 + e2 + e3 + e4

    # Fraction of each state
    y = e1 / den
    yr = e2 / den
    yo = e3 / den
    yoh = e4 / den
    # Reaction rates
    v34 = C4_CONC * (y * a12 - yr * a21)
    v35 = C4_CONC * yr * a23
    v36 = C4_CONC * (a34 * yo - a43 * yoh)
    v37 = C4_CONC * (a41 * yoh - a14 * y)

    # cytc -> complex IV electron flux
    vE = 3 * v34 + v35
    # Proton flux by C4
	# proton pumped by C1 + C4
	vH = 4 * v34 + 3 * v36 + v37
	# proton pumped by C4 alone
    vHresC4 = v34 + 2 * v36 + v37
    # O2 consumption
    vO2 = v35
    return vE, vHresC4, vO2
end

# Complex V (ATP synthase) parameters
@with_kw struct C5Params
    ρF1 = 5.0  # Concentration of complex V (mM)
    P1 = 1.346E-8
    P2 = 7.739E-7
    P3 = 6.65E-15
    PA = 1.656E-8
    PB = 3.373E-10
    PC1 = 9.651E-17
    PC2 = 4.585E-17
    KEQ = 1.71E6  # Equilibrium constant of ATP synthase
	E150 = _ra(3 * 50)
end

# Without Binding corrections
function _vaf1(pi_m, adp_m, atp_m, KEQ)
	vaf1 = KEQ * atp_m / (pi_m * adp_m)
end

# Relative electrochemical activit of ATP <-> ADP + Pi, corrected by dissociation
function _vaf1(h_m, pi_m, adp_m, atp_m, mg_m, KEQ)
    # Dissociation
    (adp3_m, hadp_m, _, poly_adp_m) = breakdown_adp(adp_m, h_m, mg_m)
    (_, _, mgatp_m, poly_atp_m) = breakdown_atp(atp_m, h_m, mg_m)
	ke_app = KEQ * h_m * poly_atp_m * _h2o_poly(h_m) / (poly_adp_m * _pi_poly(h_m))
	vaf1 = _vaf1(pi_m, adp3_m + hadp_m, mgatp_m, ke_app)
end

function _c5_rates(ΔμH, vaf1, pC5::C5Params, MT_PROT=1, C5_INHIB=1)
	@unpack ρF1, P1, P2, P3, PA, PB, PC1, PC2, KEQ, E150 = pC5
	vh = _ra(3 * ΔμH)
    common = -ρF1 * MT_PROT * C5_INHIB / ((1 + P1 * vaf1) * E150 + (P2 + P3 * vaf1) * vh)
    vATPase = common * ((100 * PA + PC1 * E150) * vaf1 - (PA + PC2 * vaf1) * vh)
    vHu = 3 * common * (100 * PA * (1 + vaf1) - vh * (PA + PB))
	return (vATPase, vHu)
end

# ATP synthase (cpmplex V) and related rates from Kembro et. al (2013)
# vATPase: ATP production rate
# vHu : protons flux through ATP synthase
function c5_rates(dpsi, ΔμH, h_m, pi_m, adp_m, conc::Concentrations, pC5::C5Params, MT_PROT=1, C5_INHIB=1)
	@unpack KEQ = pC5
	@unpack mg_m = conc
    # Potential difference of ATP <-> ADP + Pi
	vaf1 = _vaf1(h_m, pi_m, adp_m, _atp_m(adp_m, conc), mg_m, KEQ)
	return _c5_rates(ΔμH, vaf1, pC5, MT_PROT, C5_INHIB)
end

# From Li et al. 2015
function c5_rates(ΔμH, pi_m, adp_m, conc::Concentrations, pC5::C5Params, MT_PROT=1, C5_INHIB=1)
	@unpack KEQ = pC5
	vaf1 = _vaf1(pi_m, adp_m, _atp_m(adp_m, conc), KEQ)
	return _c5_rates(ΔμH, vaf1, pC5, MT_PROT, C5_INHIB)
end

@with_kw struct ANTParams
	VMAX = 0.435  # Max rate of ANT (mM/ms), (from AC Wei)
    H = 0.5  # Voltage steepness
end

# Rate of ATP/ADP translocator (ANT)
function _vant(atp4_m, adp3_m, atp4_i, adp3_i, dpsi, H, VMAX)
	f_i = atp4_i / adp3_i
	f_m = adp3_m / atp4_m
	# Magnus_Keizer_1997
	# vANT = VMAX * (1 - f_i * f_m * _ra(-dpsi)) / ((1 + f_i * _ra(-H * dpsi)) * (1 + f_m))
	# Cortassa_2003
	vANT = VMAX * (1 - f_i * f_m) / ((1 + f_i * _ra(-H * dpsi)) * (1 + f_m))
	# Li_2015
	# vANT = VMAX * _ra(-dpsi) * (1 - f_i * f_m) / ((1 + f_i * _ra(-H * dpsi)) * (1 + f_m))
end

# Rate of ATP/ADP translocator (ANT), with auto correction from dissociation
function vant(atp_m, adp_m, atp_i, adp_i, mg_m, mg_i, h_m, h_i, dpsi, H, VMAX)
	(atp4_i, _, _, _) = breakdown_atp(atp_i, h_i, mg_i)
	(atp4_m, _, _, _) = breakdown_atp(atp_m, h_m, mg_m)
	(adp3_i, _, _, _) = breakdown_adp(adp_i, h_i, mg_i)
	(adp3_m, _, _, _) = breakdown_adp(adp_m, h_m, mg_m)
	vANT = _vant(atp4_m, adp3_m, atp4_i, adp3_i, dpsi, H, VMAX)
end

function vant(adp_m, atp_i, h_m, dpsi, conc::Concentrations, pANT::ANTParams)
	@unpack H, VMAX = pANT
	@unpack h_i, mg_m, mg_i = conc
	atp_m = _atp_m(adp_m, conc)
	adp_i = _adp_i(atp_i, conc)
	vANT = vant(atp_m, adp_m, atp_i, adp_i, mg_m, mg_i, h_m, h_i, dpsi, H, VMAX)
end

# Without correction
function vant(adp_m, atp_i, dpsi, conc::Concentrations, pANT::ANTParams)
	@unpack H, VMAX = pANT
	atp_m = _atp_m(adp_m, conc)
	adp_i = _adp_i(atp_i, conc)
	atp4_i = 0.25 * atp_i
	atp4_m = 0.025 * atp_m
	adp3_i = 0.45 * adp_i
	adp3_m = 0.17 * adp_m
	return _vant(atp4_m, adp3_m, atp4_i, adp3_i, dpsi, H, VMAX)
end

# Proton leak across inner mitochondrial membrane
vhleak(ΔμH, G_H_MITO) = G_H_MITO * ΔμH

# OXPHOS (complex III)
@with_kw struct C3Params
	SCALE = 60e3  # 1/min to 1/ms
	K03 = 9.9998E4 / SCALE
    KEQ3 = 0.6877
    K04 = 3.6402E3 / SCALE
    KEQ4_OX = 129.9853
    KEQ4_RD = 13.7484
    δ₁ = 0.5
    δ₂ = 0.5
    δ₃ = 0.5
    β = 0.5006
	α = 0.5 * (1 - β)
    γ = α
    KD_Q = 1.32E6 / SCALE
    K06 = 10000 / SCALE
    KEQ6 = 9.4546
    K07_OX = 800 / SCALE
    K07_RD = 100 / SCALE
    KEQ7_OX = 3.0748
    KEQ7_RD = 29.0714
    K08_OX = 5000 / SCALE
    K08_RD = 500 / SCALE
    KEQ8_OX = 129.9853
    KEQ8_RD = 9.4546
    K09 = 4.9949E4 / SCALE
    KEQ9 = 0.2697
    K010 = 1700 / SCALE
    KEQ10 = 1.4541
    K33 = 148148 / SCALE
    KEQ33 = 2.1145
    ρC3 = 0.325
    Q_T = 4.0  # Total CoQ pool
end

_v3(A1, Q1, A2, Q2, K0, KEQ) = K0 * (KEQ * A1 * Q1 - A2 * Q2)
_v3(fes_ox, QH2_p, fes_rd, Qdot_p, K03, KEQ3, FAC_PH, MYXOTHIAZOLE) = _v3(fes_ox * FAC_PH, QH2_p, fes_rd * MYXOTHIAZOLE, Qdot_p, K03, KEQ3)
_el(α, δ, dpsi) = _ra(-α * δ * dpsi)
_er(α, δ, dpsi) = _ra(α * (1 - δ) * dpsi)
_v4(bleft, Qdot_p, el, bright, Q_p, er, K04, KEQ4) = _v3(bleft, Qdot_p * el, bright, Q_p * er, K04, KEQ4)
_v6(b2, b3, dpsi, K06, KEQ6, β, δ₂) = _v3(b2, _el(β, δ₂, dpsi), b3, _er(β, δ₂, dpsi), K06, KEQ6)
_v7(bleft, Q_n, el, bright, Qdot_n, er, K07, KEQ7, ANTIMYCIN) = _v4(bleft, Q_n, el, bright, Qdot_n * ANTIMYCIN, er, K07, KEQ7)
_v8(bleft, Qdot_n, el, bright, QH2_n, er, K08, KEQ8, FAC_PH, ANTIMYCIN) = _v7(bleft, Qdot_n * FAC_PH^2, el, bright, QH2_n, er, K08, KEQ8, ANTIMYCIN)
_v9(fes_rd, cytc1_ox, fes_ox, cytc1_rd, K09, KEQ9) = _v3(fes_rd, cytc1_ox, fes_ox, cytc1_rd, K09, KEQ9)
_v10(O2, Qdot_p, sox_m, Q_p, K010, KEQ10) = _v3(O2, Qdot_p, sox_m, Q_p, K010, KEQ10)

# OXPHOS ODE system
# TODO:
# 1. conservation of Q pool (4.0 mM)
# 1.5 Fast equilibrium of (Q_n, Q_p) and (QH2_n, QH2_p) (v2 and v5 are much faster)
# 2. conservation of BH/BL pool (same as C3_CONC)
function oxphos_system(nadh, Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, b1, b2, b3, b4, fes_ox,
	cytc1_ox, cytc_ox, dpsi, ΔμH, h_m, sox_m, suc, fum, oaa, pi_m, adp_m, atp_i, conc::Concentrations,
	pC1::C1Params, pC2::C2Params, pC3::C3Params, pC4::C4Params, pC5::C5Params,
	pDOX::DOXParams, pANT::ANTParams, MT_PROT=1)
	@unpack O2, FAC_PH = conc
	@unpack C2_INHIB, C3_INHIB, C4_INHIB = pDOX
	vc1, vROSC1 = c1_rates(nadh, Q_n, QH2_n, dpsi, h_m, sox_m, conc, pC1, pDOX, MT_PROT)
	vc2 = vSDH = vsdh_mm(Q_n, QH2_n, suc, fum, oaa, pC2, C2_INHIB)
	vE, vHresC4, vO2 = c4_rates(cytc_ox, dpsi, h_m, conc, pC4, C4_INHIB, MT_PROT)
	vATPase, vHu = c5_rates(dpsi, ΔμH, h_m, pi_m, adp_m, conc, pC5, MT_PROT)
	vANT = vant(adp_m, atp_i, h_m, dpsi, conc, pANT)
	# Complex III ODE system
	@unpack (KD_Q, K03, KEQ3, α, δ₁, K04, KEQ4_OX, KEQ4_RD, β, δ₂, KEQ6, K06, γ, δ₃,
        KEQ7_OX, K07_OX, KEQ7_RD, K07_RD, KEQ8_OX, K08_OX, KEQ8_RD, K08_RD, K09,
        KEQ9, K010, KEQ10, K33, KEQ33, ρC3) = pC3

    C3_CONC = MT_PROT * ρC3
	C4_CONC = MT_PROT * pC4.ρC4
    fes_rd = C3_CONC - fes_ox
    cytc1_rd = C3_CONC - cytc1_ox
	cytc_rd = C4_CONC - cytc_ox
	MYXOTHIAZOLE = ANTIMYCIN = 1
	# v1 = Q reduction
	v1 = vc1 + vc2
	# v2 = QH2 diffusion (n-side -> p-side)
    v2 = KD_Q * (QH2_n - QH2_p)
    # v3 = QH2 to FeS
	v3 = _v3(fes_ox, QH2_p, fes_rd, Qdot_p, K03, KEQ3, FAC_PH, MYXOTHIAZOLE)
	# v4 = Qdot_p to bH
	el = _el(α, δ₁, dpsi)
	er = _er(α, δ₁, dpsi)
	v4_ox = _v4(b1, Qdot_p, el, b2, Q_p, er, K04, KEQ4_OX)
	v4_rd = _v4(b3, Qdot_p, el, b4, Q_p, er, K04, KEQ4_RD)
    # v5 = Q diffusion (p-side -> n-side)
    v5 = KD_Q * (Q_p - Q_n)
    # v6 = bH to bL
	v6 = _v6(b2, b3, dpsi, K06, KEQ6, β, δ₂)
    # v7 = bL to Qn; v8: bL to Qdot_n
	el = _el(γ, δ₃, dpsi)
	er = _er(γ, δ₃, dpsi)
	v7_ox = _v7(b3, Q_n, el, b1, Qdot_n, er, K07_OX, KEQ7_OX, ANTIMYCIN) * C3_INHIB
	v7_rd = _v7(b4, Q_n, el, b2, Qdot_n, er, K07_RD, KEQ7_RD, ANTIMYCIN) * C3_INHIB
	# v8 = bL to Qdot_n and proton transfer
	v8_ox = _v8(b3, Qdot_n, el, b1, QH2_n, er, K08_OX, KEQ8_OX, FAC_PH, ANTIMYCIN) * C3_INHIB
	v8_rd = _v8(b4, Qdot_n, el, b2, QH2_n, er, K08_RD, KEQ8_RD, FAC_PH, ANTIMYCIN) * C3_INHIB
    # v9 = fes -> cytc1
    v9 = _v9(fes_rd, cytc1_ox, fes_ox, cytc1_rd, K09, KEQ9)
    # v10 = Qdot + O2 -> O2- + Q  (ROS produced by complex III)
    vROSC3 = v10 = _v10(O2, Qdot_p, sox_m, Q_p, K010, KEQ10)

	# Redox reaction between cytc1 and cytc
	v33 = _v3(cytc1_rd, cytc_ox, cytc1_ox, cytc_rd, K33, KEQ33)

    # ODE system
	d_nadh_ox = -vc1
    d_Q_n = v5 - v7_ox - v7_rd - v1
    d_Qdot_n = v7_ox + v7_rd - v8_ox - v8_rd
    d_QH2_n = v8_ox + v8_rd - v2 + v1
    d_QH2_p = v2 - v3
    d_Qdot_p = v3 - v10 - v4_ox - v4_rd
    d_Q_p = v10 + v4_ox + v4_rd - v5
    d_b1 = v7_ox + v8_ox - v4_ox
    d_b2 = v4_ox + v7_rd + v8_rd - v6
    d_b3 = v6 - v4_rd - v7_ox - v8_ox
    d_b4 = v4_rd - v7_rd - v8_rd
    d_fes_ox = v9 - v3
    d_cytc1_ox = v33 - v9
	d_cytc_ox = vE - v33
	vHres = 4 * vc1 + 2 * v3 + vHresC4
	return (d_nadh_ox, d_Q_n, d_Qdot_n, d_QH2_n, d_QH2_p, d_Qdot_p, d_Q_p, d_b1, d_b2,
            d_b3, d_b4, d_fes_ox, d_cytc1_ox, d_cytc_ox, vHres, vROSC1, vSDH, vROSC3, vO2,
			vATPase, vHu, vANT)
end
