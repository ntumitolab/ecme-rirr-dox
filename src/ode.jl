using Parameters
using LabelledArrays
using DifferentialEquations

import .Utils: mM, μM, Molar, ms, mV, μA, cm², iVT, p_one, iCM, A_CAP_V_MYO_F, A_CAP_V_SS_F, V_NSR, V_MYO, V_JSR, V_SS

@with_kw struct CMCParams{R}
    iStim::R = 0.0μA / cm² # External stimulus current
    k_o::R = 5.4mM         # Extracellular potassium
    na_o::R = 140.0mM      # Extracellular sodium
    ca_o::R = 2.0mM        # Extracellular calcium
    mg_i::R = 3.1mM        # Cytosolic magnesium (Li-2015)
    pi_i::R = 2.0mM        # Cytosolic inorganic phosphate (Li-2015)
    PH_I::R = 7.0          # Cytosolic pH
    h_i::R = exp10(-PH_I) * Molar # Cytosolic proton concentration
    R_TR::R = inv(9.09ms)   # Rate of Ca diffusion from NSR to JSR (1/ms)
    R_XFER::R = inv(0.5747ms)  # Rate of Ca diffusion from subspace to cytosol
    G_NAB::R = 3.22E-3mS / cm²
    G_CAB::R = 5.45E-4mS / cm²
    pCK = CK()
    pNKA = NKA(k_o = k_o, na_o = na_o)
    pPMCA = PMCA()
    pSERCA = SERCA()
    pCSQN = CaBuffer(KM_CA = 0.8mM, ET = 5.0mM)
    pCMDN = CaBuffer(KM_CA = 2.38μM, ET = 50μM)
    pLCC = LCC()
    pRyR = RyR()
    pINa = INa()
    pIK = IK(k_o = k_o)
    pNCX = NCX()
    pNSNa = NSNa()
    pMyo = MyoFibril()
    atp_i::R = 7.9
    adp_i::R = 8.0 - atp_i
    crp_i::R = 5.0
    cr_i::R = 25.0 - crp_i
end

function build_u0(t = 0;
    vm = -87.28mV,          # Sarcolemmal membrane potential
    m_na = 0.0327,          # Fast Na gating (activation)
    h_na = 0.9909,          # Fast Na gating (inactivation)
    j_na = 0.9941,          # Fast Na gating (slow inactivation)
    x_k = 1.1212e-4,        # K channel gating (activation)
    na_i = 8.2143mM,        # Cytoplasmic Na
    k_i = 150.8mM,          # Cytoplasmic K 
    ca_i = 5.8653e-5mM,     # Cytoplasmic Ca
    ca_jsr = 0.1948mM,      # Junctional SR Ca
    ca_nsr = 0.1948mM,      # Network SR Ca (mM)
    ca_ss = 7.057e-5mM,     # Subspace Ca
    ltr_ca = 8.949E-3mM,    # Ca bound to low affinity troponin sites (mM)
    htr_ca = 0.1321mM,      # Ca bound to high affinity troponin sites (mM)
    # CICR states
    po1 = 8.309e-5,
    po2 = 7.12e-11,
    pc2 = 0.2471,
    c0 = 0.9991,
    c1 = 8.175e-5,
    c2 = 2.508e-9,
    c3 = 3.421e-14,
    c4 = 1.749e-20,
    o_lcc = 2.624e-20,
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
    n1 = 2.243e-5
    # adp_i = 8 - 7.8707,    # Cytosolic ATP(mM), linked to EC coupling
    # adp_ic = 8 - 7.7107,   # Cytosolic ATP(mM), not linked to EC coupling
    # crp_i = 5.0297,        # Cytosolic creatine phosphate(mM), linked to EC coupling
    # crp_ic = 5.1291,       # Cytosolic creatine phosphate(mM), not linked to EC coupling
)
    return LVector(
        vm = vm,
        m_na = m_na,
        h_na = h_na,
        j_na = j_na,
        x_k = x_k,
        na_i = na_i,
        k_i = k_i,
        ca_i = ca_i,
        ca_jsr = ca_jsr,
        ca_nsr = ca_nsr,
        ca_ss = ca_ss,
        ltr_ca = ltr_ca,
        htr_ca = htr_ca,
        po1 = po1,
        po2 = po2,
        pc2 = pc2,
        c0 = c0,
        c1 = c1,
        c2 = c2,
        c3 = c3,
        c4 = c4,
        o_lcc = o_lcc,
        cca0 = cca0,
        cca1 = cca1,
        cca2 = cca2,
        cca3 = cca3,
        y_ca = y_ca,
        p0 = p0,
        p1 = p1,
        p2 = p2,
        p3 = p3,
        n1 = n1,
        # adp_i = adp_i,
        # adp_ic = adp_ic,
        # crp_i = crp_i,
        # crp_ic = crp_ic
    )
end

function build_pacing_callbacks(tspan; period = 1000ms, duration = 0.5ms, strength = -80.0μA / cm², proposed_dt = 1ms)
    start_stimulus = (u, t, integrator) -> begin
        integrator.p = CMCParams(iStim = strength)
        set_proposed_dt!(integrator, proposed_dt)
    end
    end_stimulus = (u, t, integrator) -> begin
        integrator.p = CMCParams(iStim = zero(strength))
        set_proposed_dt!(integrator, proposed_dt)
    end
    start_callback = FunctionCallingCallback(start_stimulus, funcat = tspan[1]:period:tspan[2])
    end_callback = FunctionCallingCallback(end_stimulus, funcat = (tspan[1]+duration):period:tspan[2])
    return CallbackSet(start_callback, end_callback)
end

# Fixed ATP model
function model!(du, u, p::CMCParams, t)
    @unpack adp_i, atp_i, crp_i, cr_i = p
    @unpack vm = u
    vfrt = vm * iVT
    evfrtm1 = expm1(vfrt)
    evfrt = p_one(evfrtm1)

    # CICR
    @unpack pLCC, pRyR, R_TR, R_XFER, k_o, ca_o = p
    @unpack c0, c1, c2, c3, c4, o_lcc, cca0, cca1, cca2, cca3, y_ca, vm, k_i, ca_nsr, ca_jsr, ca_i, ca_ss, po1, po2, pc2 = u
    rates = lcc_system(c0, c1, c2, c3, c4, o_lcc, cca0, cca1, cca2, cca3, ca_ss, vm, pLCC)
    du.c0 = rates.d_c0
    du.c1 = rates.d_c1
    du.c2 = rates.d_c2
    du.c3 = rates.d_c3
    du.c4 = rates.d_c4
    du.o_lcc = rates.d_o
    du.cca0 = rates.d_cca0
    du.cca1 = rates.d_cca1
    du.cca2 = rates.d_cca2
    du.cca3 = rates.d_cca3
    du.y_ca = d_yca(y_ca, vm)
    rates = ryr_system(po1, po2, pc2, ca_ss, pRyR)
    du.po1 = rates.d_po1
    du.po2 = rates.d_po2
    du.pc2 = rates.d_pc2

    jRel = j_rel(po1, po2, ca_jsr, ca_ss, pRyR)
    jTr = R_TR * (ca_nsr - ca_jsr)
    jXfer = R_XFER * (ca_ss - ca_i)

    @unpack iCaL, iCaK = i_cal(y_ca, o_lcc, vm, k_i, k_o, ca_o, pLCC)

    # Calcium channels
    @unpack vm, na_i, ca_i = u
    @unpack pNCX, na_o, ca_o, G_CAB = p

    iNaCa = i_naca(vm, na_i, ca_i, na_o, ca_o, pNCX, vfrt)
    iCaB = i_cab(vm, ca_i, ca_o, G_CAB)

    # Myofibril cross-bridge

    @unpack pMyo, atp_i, adp_i = p
    @unpack p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i = u
    vAM = v_am(p0, p1, p2, atp_i, adp_i, pMyo)

    rates = xb_system(p0, p1, p2, p3, n1, ltr_ca, htr_ca, ca_i, pMyo)
    du.p0 = rates.d_p0
    du.p1 = rates.d_p1
    du.p2 = rates.d_p2
    du.p3 = rates.d_p3
    du.n1 = rates.d_n1
    du.ltr_ca = rates.d_ltr_ca
    du.htr_ca = rates.d_htr_ca
    jTrpn = rates.j_trpn

    # CK shuttle: ATP fixed by now, ignored for now

    # Sodium channels
    @unpack na_i, m_na, h_na, j_na, vm, ca_i = u
    @unpack pINa, pNSNa, G_NAB, na_o = p

    du.m_na = d_mna(m_na, vm)
    du.h_na = d_hna(h_na, vm)
    du.j_na = d_jna(j_na, vm)

    ΔvNa = vm - nernst(na_o, na_i)
    iNa = i_na(m_na, h_na, j_na, ΔvNa, pINa)
    iNsNa = i_nsna(vm, na_i, ca_i, na_o, pNSNa)
    iNaB = i_nab(ΔvNa, G_NAB)

    # Potassium channels
    @unpack k_i, na_i, x_k, vm = u
    @unpack adp_i, atp_i, k_o, na_o, pIK, mg_i = p

    du.x_k = d_xk(x_k, vm)
    ΔvK = vm - nernst(k_o, k_i)
    iK1 = i_k1(ΔvK, pIK)
    iK = i_k(x_k, vm, na_i, k_i, na_o, pIK)
    iKp = i_kp(vm, ΔvK, pIK)
    iKatp = i_katp(adp_i, atp_i, vm, na_i, mg_i, ΔvK, pIK)

    # Ion pumps
    @unpack ca_i, ca_nsr, na_i, vm = u
    @unpack adp_i, atp_i, pNKA, pPMCA, pSERCA = p
    iNaK = i_nak(na_i, atp_i, adp_i, vm, pNKA, evfrt)
    iPCa = i_pca(ca_i, atp_i, adp_i, pPMCA)
    jUp = j_up(ca_i, ca_nsr, atp_i, adp_i, pSERCA)

    # Remaining ODEs
    @unpack iStim = p
    du.vm = -iCM * (iNa + iCaL + iK + iK1 + iKp + iNaCa + iNaK + iNsNa + iPCa + iCaB + iNaB + iStim + iKatp + iCaK)
    du.na_i = -A_CAP_V_MYO_F * (iNa + iNaB + iNsNa + 3 * (iNaCa + iNaK))
    du.k_i = -A_CAP_V_MYO_F * (iK + iK1 + iKp + iKatp + iStim - 2 * iNaK + iCaK)
    # du.atp_i = 0  # ATP is fixed
    # du.atp_ic = 0  # ATP is fixed
    # du.crp_i = 0  # CrP is fixed
    # du.crp_ic = 0  # CrP is fixed

    @unpack pCSQN, pCMDN = p
    @unpack ca_i, ca_nsr, ca_jsr, ca_ss = u
    du.ca_i = β_ca(ca_i, pCMDN) * (jXfer - jUp - jTrpn - 0.5 * A_CAP_V_MYO_F * (iPCa + iCaB - 2 * iNaCa))
    du.ca_nsr = (V_MYO * jUp - V_JSR * jTr) / V_NSR
    du.ca_jsr = β_ca(ca_jsr, pCSQN) * (jTr - jRel)
    du.ca_ss = β_ca(ca_ss, pCMDN) * (V_JSR / V_SS * jRel - V_MYO / V_SS * jXfer - 0.5 * iCaL * A_CAP_V_SS_F)

    return (;
        jRel,
        jTr,
        jXfer,
        iCaL,
        iCaK,
        iNaCa,
        iCaB,
        vAM,
        jTrpn,
        iNa,
        iNsNa,
        iNaB,
        iK1,
        iK,
        iKp,
        iKatp,
        iNaK,
        iPCa,
        jUp
    )
end