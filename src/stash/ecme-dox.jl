# The simplified model of ECME-RIRR plus ROS generation by ETC ODE system
using Parameters, StaticArrays

include("force.jl")
include("ik2009.jl")
include("ina.jl")
include("inaca2009.jl")
include("ins_bg.jl")
include("katp.jl")
include("mtdna.jl")
include("oxphos-ode.jl")
include("pump.jl")
include("ros.jl")
include("tca.jl")
include("transmito.jl")
include("cabuf.jl")
include("concentrations.jl")

# Initial conditions (largely from Zhou et al. (2009))
u0_tup = (
    vm = -87.28,  # Sarcolemmal membrane potential (mV)
    dpsi = 175.318,  # Mitochondrial membrane potential (mV)
    m_na = 0.0327,  # Fast Na gating (activation)
    h_na = 0.9909,  # Fast Na gating (inactivation)
    j_na = 0.9941,  # Fast Na gating (slow inactivation)
    x_k = 1.1212e-4,  # K gating (activation)
    na_i = 8.2143,  # Cytoplasmic Na (mM)
    k_i = 150.8,  # Cytoplasmic K (mM)
    ca_i = 5.8653e-5,  # Cytoplasmic Ca (mM)
    ca_jsr = 0.1948,  # Junctional SR Ca (mM)
    ca_nsr = 0.1948,  # Network SR Ca (mM)
    ca_ss = 7.057e-5,  # Subspace calcium
    ltr_ca = 8.949E-3,  # Ca bound to low affinity troponin sites (mM)
    htr_ca = 0.1321,  # Ca bound to high affinity troponin sites (mM)
    # CICR states
    po1 = 8.309e-5,
    po2 = 7.12e-11,
    pc2 = 0.2471,
    c0 = 0.9991,
    c1 = 8.175e-5,
    c2 = 2.508e-9,
    c3 = 3.421e-14,
    c4 = 1.749e-20,
    o = 2.624e-20,
    cca0 = 1.1328e-3,
    cca1 = 1.7591e-8,
    cca2 = 6.3826e-13,
    cca3 = 1.815e-15,
    y_ca = 0.9479,
    # Tropomyosin cross bridge states
    p0 = 2.601e-5,
    p1 = 2.248e-5,
    p2 = 4.199e-5,
    p3 = 3.657e-5,
    n1 = 2.243e-5,
    adp_i = 8 - 7.8707,  # Cytosolic ATP(mM), linked to EC coupling
    adp_ic = 8 - 7.7107,  # Cytosolic ATP(mM), not linked to EC coupling
    crp_i = 5.0297,  # Cytosolic creatine phosphate(mM), linked to EC coupling
    crp_ic = 5.1291,  # Cytosolic creatine phosphate(mM), not linked to EC coupling
    # TCA cycle states
    ca_m = 2.190477831750499e-5,
    adp_m = 0.005788892127721868,
    nadh = 9.267368089681039,
    isoc = 0.4858673722382066,
    akg = 0.0001907868710833752,
    scoa = 0.28625004579156677,
    suc = 0.00011243840322880198,
    fum = 0.005275867239849815,
    mal = 0.002781903362877464,
    oaa = 2.683029248919834e-8,
    # ROS
    sox_m = 0.0004966725362784072,
    sox_i = 8.570060710246488e-9,
    h2o2_i = 2.061879987223193e-9,
    gsh_i = 0.9999899536989296,
    # ETC states
    Q_n = 1.3268400422082767,
    Qdot_n = 0.06757658200311795,
    QH2_n = 0.6002709581315435,
    QH2_p = 0.6001505595245672,
    Qdot_p = 0.07820141736951611,
    Q_p = 1.3269604407629754,
    b1 = 0.11362728970178769,
    b2 = 0.14764414671393097,
    b3 = 0.04917170786399363,
    b4 = 0.014656855720284924,
    fes_ox = 0.054765936912457895,
    cytc1_ox = 0.16423864061514593,
    cytc_ox = 0.1081945884323241)

nameLUT = Dict((k=>i) for (i, k) in enumerate(keys(u0_tup)))

# TODO: conservation of Q pool
Q_POOL = u0_tup[:Q_n] + u0_tup[:Qdot_n] + u0_tup[:QH2_n] + u0_tup[:QH2_p] + u0_tup[:Qdot_p] + u0_tup[:Q_p]
# TODO: b pool = C3 conc
B_POOL = u0_tup[:b1] + u0_tup[:b2] + u0_tup[:b3] + u0_tup[:b4]

conc = Concentrations(nadph_i=1.0,
                      pi_m=2.0,
                      mg_m=0.4,
                      mg_i=3.1,
                      ΣA_m=1.5,
                      ΣNAD_m=10.0,
                      ΣA_i=8.0,
                      O2=6E-3)

@with_kw struct Params
    conc = conc
    pSODi = SODParams(K1=1.2E3, K3=24.0, K5=2.5E-4, ET=1.43E-3, KI_H2O2=0.5)
    pGPXi = GPXParams(ET=0.01, 𝚽1=5e-3, 𝚽2=0.75)
    pGRi = GRParams(K1=5E-3, ET=0.01, KM_GSSG=0.06, KM_NADPH=0.015, GSH_T=1.0)
    pCAT = CATParams(K1=17.0, ET=0.01, FR=0.05)
    pC1 = C1Params()
    pC2 = C2Params()
    pC3 = C3Params()
    pC4 = C4Params()
    pC5 = C5Params(ρF1=0.05)
    pANT = ANTParams(VMAX=5E-3)
    pMNCE = MNCEParams(VMAX=7.3E-4)
    pMCU = MCUParams(VMAX=0.028)
    pIMAC = IMACParams()
    pCS = CSParams(ACCOA=1.0, KCAT=0.5, ET=0.4, KM_ACCOA = 0.0126, KM_OAA = 6.4E-4)
    pACO = ACOParams(TCA_T=1.0, KF=1.25E-2, KEQ=2.22)
    pIDH = IDHParams(KM_ADP=0.62, KM_CA=5E-4, KI_NADH=0.19, KCAT=0.05, ET=0.109, KH1=8.1e-5, KH2=5.98E-5, KM_ISOC=1.52, NI_ISOC=2, KM_NAD=0.923)
    pKGDH = KGDHParams(KM_MG=0.0308, KM_CA=1.27E-3, ET=0.5, KCAT=0.075, KM_AKG=1.94, KM_NAD=38.7, NI_AKG=1.2)
    pSL = SLParams(KF=5E-3, KEQ=3.115, COA=0.02)
    pFH = FHParams(KF=3.32E-3, KEQ=1.0)
    pMDH = MDHParams(KH1=1.131E-5, KH2=26.7, KH3=6.68E-9, KH4=5.62E-6, K_OFFSET=3.99E-2, KCAT=0.11, ET=0.15, KM_MAL=1.49, KI_OAA=3.1E-3, KM_NAD=0.22)
    pAAT = AATParams(KF=6.44E-4, KEQ=6.6, KASP=1.5E-6, GLU=10.0)
    pDOX = DOXParams()
    pMTDNA = MTDNAParams()
    # G_H_MITO = 2E-6  # Mitochondrial membrane conductance (Gauthier 2013)
    G_H_MITO = 1E-8  # (Li 2015)
    δCA = 3E-4  # Ca buffering factor
    R_TR= inv(9.09) # Rate of Ca diffusion from NSR to JSR (1/ms)
    R_XFER = inv(0.5747)  # Rate of Ca diffusion from subspace to cytosol (1/ms)
    G_NA = 12.8  # Sodium conductance
    V_NSR = 1.4E-6  # Network SR Volume (μL)
    V_JSR = 0.16E-6  # Junctional SR Volume (μL)
    V_SS = 0.495E-9  # Total subspace Volume (μL)
    V_MYO = 25.84E-6  # Cytoplasm Volume (μL)
    V_MITO = 15.89E-6  # Mitochondrial Volume (μL)
    A_CAP = 1.534E-4  # Cell capacitance area (cm²)
    A_CAP_V_MYO_F = A_CAP / (V_MYO * F)
    V_MITO_V_MYO = V_MITO / V_MYO
    pIK = IKParams()
    pCK = CKParams()
    pXB = XBParams()
    pLCC = LCCParams()
    pRYR = RYRParams()
    pSERCA = SERCAParams(KMF_CA=2.4e-4, KMR_CA=1.64269, NFB=1.4, NRB=1)
    pPMCA = PMCAParams()
    pNS = NSParams()
    pNKA = NKAParams()
    pNCX = NCXParams()
    pKATP = KATPParams()
    pCSQN = CABUFParams(KM_CA=0.8, ET=5.0)
    pCMDN = CABUFParams(KM_CA=2.38E-3, ET=5E-2)
    pBG = BGParams()
    BCL = 0
end

function _istim(t, bcl=0, peroid=0.5, strength=-80.0)
    if (bcl > 0) && (rem(t, bcl) <= peroid)
        return strength
    else
        return zero(strength)
    end
end

function ecme_dox(u, p, t)
    (vm, dpsi, m_na, h_na, j_na, x_k, na_i, k_i, ca_i, ca_jsr, ca_nsr, ca_ss, ltr_ca, htr_ca,
     po1, po2, pc2, c0, c1, c2, c3, c4, o, cca0, cca1, cca2, cca3, y_ca, p0, p1, p2, p3, n1,
     adp_i, adp_ic, crp_i, crp_ic, ca_m, adp_m, nadh, isoc, akg, scoa, suc, fum, mal, oaa,
     sox_m, sox_i, h2o2_i, gsh_i, Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, b1, b2, b3, b4,
     fes_ox, cytc1_ox, cytc_ox) = u

    @unpack conc = p
	@unpack ca_o, na_o, k_o, ΣA_i, h_i, h_m, mg_i, mg_m, ΔpH, FAC_PH, O2, nadph_i, pi_m, ΣA_i = conc
    vfrt = vm * F_RT
    evfrtm1 = expm1(vfrt)
	ΔμH = _ΔμH(dpsi, ΔpH)  # proton motive force
    atp_i = ΣA_i - adp_i
    atp_ic = ΣA_i - adp_ic

    # CICR
    @unpack pLCC, pRYR, R_TR, R_XFER = p
	(d_c0, d_c1, d_c2, d_c3, d_c4, d_o, d_cca0, d_cca1, d_cca2,
     d_cca3) = lcc_system(c0, c1, c2, c3, c4, o, cca0, cca1, cca2, cca3, ca_ss, vm, pLCC)
    d_y_ca = dyca(y_ca, vm)
    (d_po1, d_po2, d_pc2) = ryr_system(po1, po2, pc2, ca_ss, pRYR)
    (iCaL, iCaK) = ical(y_ca, o, vfrt, evfrtm1, k_i, k_o, ca_o, pLCC)
    jRel = jrel(po1, po2, ca_jsr, ca_ss, pRYR.R_RYR)
    jTr = R_TR * (ca_nsr - ca_jsr)
    jXfer = R_XFER * (ca_ss - ca_i)
    # Myofibril
    @unpack pCK, pXB = p
    (d_p0, d_p1, d_p2, d_p3, d_n1, d_ltr_ca, d_htr_ca) = xb_system(p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i, pXB)
	vAM = vam(p0, p1, p2, atp_i, adp_i, pXB)
    jTrpn = d_ltr_ca + d_htr_ca

    # CK shuttle
    (d_atp_in_ck, d_atp_ic, d_crp_i, d_crp_ic) = ck_system(atp_i, atp_ic, crp_i, crp_ic, pCK, ΣA_i)

    # Other ion channels
    @unpack G_CAB, G_NAB = p.pBG
    @unpack G_KP, P_NA_K, G_K, G_K1 = p.pIK
    @unpack pNKA, pPMCA, pSERCA, pKATP, G_NA, pNS, pNCX = p
    d_m_na = d_mna(m_na, vm)
    d_h_na = d_hna(h_na, vm)
    d_j_na = d_jna(j_na, vm)
    eNa = _ena(na_o, na_i)
    iNa = ina(vm, eNa, m_na, h_na, j_na, G_NA)
    d_x_k = dxk(x_k, vm)
    iK = ik(vm, _enak(k_o, na_o, k_i, na_i, P_NA_K), x_k, G_K)
	eK = _ek(k_o, k_i)
    iKp = ikp(vm, eK, G_KP)
    iK1 = ik1(vm, eK, G_K1)
	iKatp = ikatp(vm, eK, na_i, atp_i, adp_i, mg_i, pKATP)
    iNsNa = insna(na_i, ca_i, vfrt, evfrtm1, na_o, pNS)
    iNaCa = inaca(na_i, ca_i, vfrt, evfrtm1, ca_o, pNCX)
    iCaB = icab(vm, _eca(ca_o, ca_i), G_CAB)
    iNaB = inab(vm, eNa, G_NAB)
    iNaK = inak(na_i, atp_i, adp_i, vm, evfrtm1, pNKA)
    iPCa = ipca(ca_i, atp_i, adp_i, pPMCA)
    jUp = jup(ca_i, ca_nsr, atp_i, adp_i, pSERCA)

    # Ion fluxes across mitochondria
    @unpack pMCU, G_H_MITO, pMNCE = p
    vUni = vuni(ca_i, dpsi, pMCU)
    vNaCa = vnaca(ca_i, ca_m, na_i, dpsi, pMNCE, pMCU)
    vHleak = vhleak(ΔμH, G_H_MITO)

    # ETC
    @unpack pC1, pC2, pC3, pC4, pC5, pDOX, pANT = p
    @unpack C2_INHIB, C3_INHIB, C4_INHIB = pDOX
    @unpack MT_PROT = p.pMTDNA
    vANT = vant(adp_m, atp_i, dpsi, conc, pANT)
    vc1, vROSC1 = c1_rates(nadh, Q_n, QH2_n, dpsi, h_m, sox_m, conc, pC1, pDOX, MT_PROT)
    vc2 = vSDH = vsdh_mm(Q_n, QH2_n, suc, fum, oaa, pC2, C2_INHIB)
    vE, vHresC4, vO2 = c4_rates(cytc_ox, dpsi, h_m, conc, pC4, C4_INHIB, MT_PROT)
    (vATPase, vHu) = c5_rates(ΔμH, pi_m, adp_m, conc, pC5, MT_PROT)
    # Complex III rates
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

    # ODE system of complex III
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

	vROS = vROSC1 + vROSC3
    vHres = 4 * vc1 + 2 * v3 + vHresC4

    # TCA cycle rates
    @unpack pSL, pCS, pACO, pIDH, pKGDH, pFH, pMDH, pAAT = p
    nad = _nad(nadh, conc)
	vCS  = vcs(oaa, pCS)
	vACO = vaco(isoc, akg, scoa, suc, fum, mal, oaa, pACO)
	vIDH = vidh(isoc, nadh, nad, adp_m, ca_m, h_m, mg_m, pIDH)
	vKGDH = vkgdh(akg, nadh, nad, ca_m, mg_m, pKGDH)
	vSL = vsl(scoa, suc, adp_m, _atp_m(adp_m, conc), pSL)
	vFH = vfh(fum, mal, pFH)
	vMDH = vmdh(mal, oaa, nadh, nad, h_m, pMDH)
	vAAT = vaat(akg, oaa, pAAT)

    # ODE system of TCA cycle
	d_nadh = -vc1 + vIDH + vKGDH + vMDH
	d_isoc = vACO - vIDH
	d_akg = vIDH + vAAT - vKGDH
	d_scoa = vKGDH - vSL
	d_suc = vSL - vSDH
	d_fum = vSDH - vFH
	d_mal = vFH - vMDH
	d_oaa = vMDH - vCS - vAAT
    d_adp_m = vANT - vSL - vATPase

    @unpack (pSODi, pGPXi, pGRi, pCAT, pIMAC, V_MITO_V_MYO) = p
    vSOD_i = vsod(sox_i, h2o2_i, pSODi)
    vGPX_i = vgpx(h2o2_i, gsh_i, pGPXi)
    vGR_i = vgr(gsh_i, nadph_i, pGRi)
    vCAT = v_cat(h2o2_i, pCAT)
    vIMAC, vTrROS = vimac_vtrros(dpsi, sox_m, sox_i, pIMAC)

    # ODE system of ROS scavenging
    d_sox_m = vROS - vTrROS
    d_sox_i = V_MITO_V_MYO * vTrROS - vSOD_i
	d_h2o2_i = 0.5 * vSOD_i - vGPX_i - vCAT
	d_gsh_i  = vGR_i - vGPX_i

	# Remaining ODE equations
	@unpack (δCA, A_CAP_V_MYO_F, V_MYO, V_SS, V_NSR, V_JSR, V_MITO_V_MYO, BCL) = p
	iStim = _istim(t, BCL)
	d_vm = -CM_INV * (iNa + iCaL + iK + iK1 + iKp + iNaCa + iNaK + iNsNa + iPCa + iCaB + iNaB + iStim + iKatp + iCaK)
	d_na_i = -A_CAP_V_MYO_F * (iNa + iNaB + iNsNa + 3 * (iNaCa + iNaK)) # + V_MITO_V_MYO * (vNHE - 3 * vNaCa)
    d_k_i = -A_CAP_V_MYO_F * (iK + iK1 + iKp + iKatp + iStim - 2 * iNaK + iCaK)
	d_atp_i = d_atp_in_ck + V_MITO_V_MYO * vANT - vAM - 0.5 * jUp - A_CAP_V_MYO_F * (iPCa + iNaK)

	@unpack pCSQN, pCMDN, A_CAP = p
	d_ca_i = _β_ca(ca_i, pCMDN) * (jXfer - jUp - jTrpn - 0.5 * A_CAP_V_MYO_F * (iPCa + iCaB - 2 * iNaCa) + V_MITO_V_MYO * (vNaCa - vUni))
	d_ca_nsr = _β_ca(ca_nsr, pCSQN) * (V_MYO * jUp - V_JSR * jTr) / V_NSR
    d_ca_jsr = _β_ca(ca_jsr, pCSQN) * (jTr - jRel)
    d_ca_ss = _β_ca(ca_ss, pCMDN) * (V_JSR * jRel - V_MYO * jXfer - 0.5 * A_CAP / F * iCaL) / V_SS
    d_ca_m = δCA * (vUni - vNaCa)

    d_dpsi = CM_MITO_INV * (vHres - vHu - vANT - vHleak - vNaCa - 2 * vUni - vIMAC)
    # @show t, na_i, m_na, h_na, j_na, d_m_na, d_h_na, d_j_na, iNa

	# Disable rates for debugging
    #=
	d_vm = 0
	d_m_na = d_h_na = d_j_na = d_na_i = 0
	d_x_k = d_k_i = 0
	d_p0 = d_p1 = d_p2 = d_p3 = d_n1 = 0
	d_ca_i = d_ca_jsr = d_ca_nsr = d_ca_ss = d_ltr_ca = d_htr_ca = 0
	d_atp_i = d_atp_ic = d_crp_i = d_crp_ic = 0
    d_po1 = d_po2 = d_pc2 = d_c0 = d_c1 = d_c2 = d_c3 = d_c4 = d_o = d_cca0 = d_cca1 = d_cca2 = d_cca3 = d_y_ca = 0
    =#

    d_adp_i = -d_atp_i
    d_adp_ic = -d_atp_ic
    du = @SVector[d_vm, d_dpsi, d_m_na, d_h_na, d_j_na, d_x_k, d_na_i, d_k_i, d_ca_i, d_ca_jsr, d_ca_nsr, d_ca_ss, d_ltr_ca, d_htr_ca,
     d_po1, d_po2, d_pc2, d_c0, d_c1, d_c2, d_c3, d_c4, d_o, d_cca0, d_cca1, d_cca2, d_cca3, d_y_ca, d_p0, d_p1, d_p2, d_p3, d_n1,
     d_adp_i, d_adp_ic, d_crp_i, d_crp_ic, d_ca_m, d_adp_m, d_nadh, d_isoc, d_akg, d_scoa, d_suc, d_fum, d_mal, d_oaa,
     d_sox_m, d_sox_i, d_h2o2_i, d_gsh_i, d_Q_n, d_Qdot_n, d_QH2_n, d_QH2_p, d_Qdot_p, d_Q_p, d_b1, d_b2, d_b3, d_b4,
     d_fes_ox, d_cytc1_ox, d_cytc_ox,
     iCaL, iCaK, jRel, jTr, jXfer, vAM, jTrpn, iNa, iK, iKp, iK1, iKatp, iNsNa, iNaCa, iCaB,
     iNaB, iNaK, iPCa, jUp, vUni, vNaCa, vHleak, vANT, vc1, vROSC1, vSDH, vE, vHresC4, vO2,
     vATPase, vHu, vROSC3, vHres, vCS, vACO, vIDH, vKGDH, vSL, vFH, vMDH, vAAT, vSOD_i, vGPX_i,
     vGR_i, vCAT, vIMAC, vTrROS]
end

rateSymbols = [:d_vm, :d_dpsi, :d_m_na, :d_h_na, :d_j_na, :d_x_k, :d_na_i, :d_k_i, :d_ca_i, :d_ca_jsr, :d_ca_nsr, :d_ca_ss, :d_ltr_ca, :d_htr_ca,
 :d_po1, :d_po2, :d_pc2, :d_c0, :d_c1, :d_c2, :d_c3, :d_c4, :d_o, :d_cca0, :d_cca1, :d_cca2, :d_cca3, :d_y_ca, :d_p0, :d_p1, :d_p2, :d_p3, :d_n1,
 :d_adp_i, :d_adp_ic, :d_crp_i, :d_crp_ic, :d_ca_m, :d_adp_m, :d_nadh, :d_isoc, :d_akg, :d_scoa, :d_suc, :d_fum, :d_mal, :d_oaa,
 :d_sox_m, :d_sox_i, :d_h2o2_i, :d_gsh_i, :d_Q_n, :d_Qdot_n, :d_QH2_n, :d_QH2_p, :d_Qdot_p, :d_Q_p, :d_b1, :d_b2, :d_b3, :d_b4,
 :d_fes_ox, :d_cytc1_ox, :d_cytc_ox,
 :iCaL, :iCaK, :jRel, :jTr, :jXfer, :vAM, :jTrpn, :iNa, :iK, :iKp, :iK1, :iKatp, :iNsNa, :iNaCa, :iCaB,
 :iNaB, :iNaK, :iPCa, :jUp, :vUni, :vNaCa, :vHleak, :vANT, :vc1, :vROSC1, :vSDH, :vE, :vHresC4, :vO2,
 :vATPase, :vHu, :vROSC3, :vHres, :vCS, :vACO, :vIDH, :vKGDH, :vSL, :vFH, :vMDH, :vAAT, :vSOD_i, :vGPX_i,
 :vGR_i, :vCAT, :vIMAC, :vTrROS]

rateMap = Dict((k => i) for (i, k) in enumerate(rateSymbols))

# Ignore the intermediate rates
rhs(u, p, t) = ecme_dox(u, p, t)[1:length(u)]

reject_step(u, p, t) = any( x -> (x < 0.0), @view u[3:end])
