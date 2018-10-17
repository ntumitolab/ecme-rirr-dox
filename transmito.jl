#Exchangers and transporters on the inner mitochondirla membrane
using Parameters
include("common.jl")

# NHE (mitochondrial Na-H exchanger) parameters
@with_kw struct NHEParams
	# Intracelllular proton
	H_I = 1E-4
    # NHE forward rate constant
    K1_P = 0.0252
    # NHE backward rate constant
    K1_M = 0.0429
    # NHE forward rate constant
    K4_P = 0.16
    # NHE backward rate constant
    K4_M = 0.0939
    # Na dissociation constant
    KM_NA = 24.0
    # H+ dissociation constant
    KM_H = 1.584E-4
	# Proton inhibitory constant
	H_KI = 3.02E-6
    # Hill coefficient for proton binding
    NI = 3
    # NHE concentration
    ET = 0.00785
    # Dependence factor of mitochondrial proton
    DENOM = 1 + (H_KI / H_I)^ NI
	VMAX = ET / DENOM
end

# Mitochondrial Na-H exchanger (NHE). scalar version
function vnhe(na_m, na_i, h_m, h_i, K1_P, K1_M, K4_P, K4_M, KM_H, KM_NA, VMAX)
	β₁⁺ = K1_P * _mm(na_m * _mm(KM_H, h_m), KM_NA)
	β₁⁻ = K1_M * _mm(na_i * _mm(KM_H, h_i), KM_NA)
	β₂⁺ = K4_P * _mm(h_i * _mm(KM_NA, na_i), KM_H)
	β₂⁻ = K4_M * _mm(h_m * _mm(KM_NA, na_m), KM_H)
    vNHE = VMAX * (β₁⁺ * β₂⁺ - β₁⁻ * β₂⁻) / (β₁⁺ + β₂⁺ + β₁⁻ + β₂⁻)
end

function vnhe(na_m, na_i, h_m, h_i, pNHE::NHEParams)
	@unpack K1_P, K1_M, K4_P, K4_M, KM_H, KM_NA, VMAX = pNHE
	vNHE = vnhe(na_m, na_i, h_m, h_i, K1_P, K1_M, K4_P, K4_M, KM_H, KM_NA, VMAX)
end

# PIC (phosphate carrier)
@with_kw struct PICParams
    # PIC concentration (mM)
    ET = 4.9
    # Cytosolic Pi binding constant
    KM_PI_IN = 11.06
    # Mitochondrial Pi binding constant
    KM_PI_MT = 11.06
    # Cytosolic OH binding constant
    KM_OH_IN = 4.08E-5
    # Mitochondrial OH binding constant
    KM_OH_MT = 4.08E-5
    # Max forward rate of PIC (mM/ms)
    VF = 90 / 60E3
    VF_MAX = VF * ET
    # Max backward rate of PIC (mM/ms)
    VB = 90 / 60E3
    VB_MAX = VB * ET
end

# Phosphate carrier (PiC)
# Manual correction from dissociation
function _vpic(H₂PO₄⁻_m, H₂PO₄⁻_i, h_m, h_i, KM_PI_IN, KM_PI_MT, KM_OH_IN, KM_OH_MT, VF_MAX, VB_MAX)
	f_pi = H₂PO₄⁻_i / KM_PI_IN
    f_pm = H₂PO₄⁻_m / KM_PI_MT
    f_ohi = _oh(h_i) / KM_OH_IN
    f_ohm = _oh(h_m) / KM_OH_MT
	ϕf = f_pi * f_ohm
	ϕb = f_pm * f_ohi
    den = 1 + f_pi + f_pm + f_ohi + f_ohm + ϕf + ϕb
    vPIC = (VF_MAX * ϕf - VB_MAX * ϕb) / den
end

# Auto correction from dissociation
function vpic(pi_m, pi_i, h_m, h_i, KM_PI_IN, KM_PI_MT, KM_OH_IN, KM_OH_MT, VF_MAX, VB_MAX)
	H₂PO₄⁻_m = _H2PO4(pi_m, h_m)
	H₂PO₄⁻_i = _H2PO4(pi_i, h_i)
	vPIC = _vpic(H₂PO₄⁻_m, H₂PO₄⁻_i, h_m, h_i, KM_PI_IN, KM_PI_MT, KM_OH_IN, KM_OH_MT, VF_MAX, VB_MAX)
end

function vpic(pi_m, pi_i, h_m, h_i, pPIC::PICParams)
	@unpack KM_PI_IN, KM_PI_MT, KM_OH_IN, KM_OH_MT, VF_MAX, VB_MAX = pPIC
	vPIC = vpic(pi_m, pi_i, h_m, h_i, KM_PI_IN, KM_PI_MT, KM_OH_IN, KM_OH_MT, VF_MAX, VB_MAX)
end

# Mitochondrial NCE parameters
@with_kw struct MNCEParams
    VMAX = 4.665E-5
    B = 0.5
    KM_NA = 9.4
    KM_CA = 3.75E-4
    N = 3
end

# Mitochondrial Ca uniporter
@with_kw struct MCUParams
	VMAX = 1.2295E-4
    DPSI_OFFSET = 91
    KACT = 3.8E-4
    KTRANS = 0.019
    L = 110
    N = 2.8
end

# Mitochondrial Calcium uniporter (MCU) flux (mM/ms)
# scalar version
function vuni(ca_i, dpsi, DPSI_OFFSET, KTRANS, KACT, N, L, VMAX)
	v2frt = 2 * F_RT * (dpsi - DPSI_OFFSET)
    f_ca = ca_i / KTRANS
    f_ca_p1 = f_ca + one(f_ca)
	f_act =  L * (_mm(KACT, ca_i) ^ N)
	vUni = VMAX * f_ca * f_ca_p1^3 * v2frt / ((f_ca_p1^4 + f_act) * (-expm1(-v2frt)))
end

function vuni(ca_i, dpsi, pMCU::MCUParams)
	@unpack DPSI_OFFSET, KTRANS, KACT, N, L, VMAX = pMCU
	vUni = vuni(ca_i, dpsi, DPSI_OFFSET, KTRANS, KACT, N, L, VMAX)
end

# Mitochondrial NCE (mNCE) flux (mM/ms)
function vnaca(ca_i, ca_m, na_i, dpsi, KM_NA, KM_CA, N, B, VMAX, DPSI_OFFSET)
	f_na = _mm(na_i, KM_NA) ^ N
    f_ca = _mm(ca_m, KM_CA)
	ϕ_ca = ca_m / ca_i
	vNaCa = VMAX * _ra(B * (dpsi - DPSI_OFFSET)) * ϕ_ca * f_na * f_ca
end

function vnaca(ca_i, ca_m, na_i, dpsi, pNCE::MNCEParams, pMCU::MCUParams)
	@unpack DPSI_OFFSET = pMCU
	@unpack KM_NA, KM_CA, N, B, VMAX = pNCE
	vNaCa = vnaca(ca_i, ca_m, na_i, dpsi, KM_NA, KM_CA, N, B, VMAX, DPSI_OFFSET)
end
