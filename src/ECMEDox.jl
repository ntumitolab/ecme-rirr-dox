module ECMEDox

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using NaNMath

export build_model, build_u0

include("utils.jl")
include("ck.jl")
include("nakpump.jl")
include("cicr16.jl")
include("ina.jl")
include("ik2009.jl")
include("force.jl")
include("ica.jl")
include("tca.jl")
include("ros.jl")
include("oxphos.jl")
include("transmito.jl")

"Initial conditions after 1000 seconds of 1Hz pacing"
function build_u0(sys)
    return [
        sys.sox_i => 1.0222129545264641nM,
        sys.h2o2_i => 0.755888231084137nM,
        sys.gssg_i => 2.220885967007533μM,
        sys.Q_n => 1.7555277577733173mM,
        sys.Qdot_n => 0.250155168931016mM,
        sys.QH2_n => 0.10775990996722686mM,
        sys.QH2_p => 0.10771653709539612mM,
        sys.Qdot_p => 0.023269397495713263mM,
        sys.cytb_1 => 0.19471601197446434mM,
        sys.cytb_2 => 0.08033504811661059mM,
        sys.cytb_3 => 0.047309996276701564mM,
        sys.fes_ox => 0.20066008052559806mM,
        sys.cytc1_ox => 0.28338319859718647mM,
        sys.cytc_ox => 0.24903758981942076mM,
        sys.isoc => 0.05159054318098895mM,
        sys.akg => 0.051145197718677655mM,
        sys.scoa => 0.03508849487000582mM,
        sys.suc => 0.0019107469302081612mM,
        sys.fum => 0.1751906841603877mM,
        sys.mal => 0.15856757152954906mM,
        sys.oaa => 0.011576938766421891mM,
        sys.m_na => 0.0010899662356754164,
        sys.h_na => 0.9905304527574307,
        sys.j_na => 0.9892354695981638,
        sys.x_k => 0.0012736748778708119,
        sys.ltr_ca => 0.02049022737620247mM,
        sys.htr_ca => 0.1372356208938611mM,
        sys.x_p0 => 0.0036716635421201937,
        sys.x_p1 => 0.0033538440069086255,
        sys.x_p2 => 0.006345552748082013,
        sys.x_p3 => 0.005553300147495315,
        sys.x_n1 => 0.012559971511989917,
        sys.adp_ic => 0.34296593237569856mM,
        sys.crp_i => 4.3883068168196715mM,
        sys.crp_ic => 4.376369356028448mM,
        sys.po1_ryr => 0.0006784597633442083,
        sys.po2_ryr => 6.977898940703352e-9,
        sys.pc2_ryr => 0.5676773131946041,
        sys.c1_lcc => 9.238521820167919e-6,
        sys.c2_lcc => 3.2118717936223875e-11,
        sys.c3_lcc => 4.962840352619686e-17,
        sys.c4_lcc => 2.919740376980712e-23,
        sys.o_lcc => 5.055151129327962e-24,
        sys.cca0_lcc => 0.0034892823303371255,
        sys.cca1_lcc => 1.2939760235285155e-7,
        sys.cca2_lcc => 1.7994436509185547e-12,
        sys.cca3_lcc => 1.112151317794485e-17,
        sys.cca4_lcc => 1.0637011854321615e-22,
        sys.x_yca => 0.737396088676386,
        sys.na_i => 15.360882167675042mM,
        sys.k_i => 142.39854398555977mM,
        sys.ca_i => 0.00016134091543239482mM,
        sys.ca_nsr => 0.9944103502674689mM,
        sys.ca_jsr => 0.9760651184113716mM,
        sys.ca_ss => 0.000169823015720396mM,
        sys.ca_m => 0.0002427548693370271mM,
        sys.adp_i => 0.19176238251620303mM,
        sys.adp_m => 0.12562731414699957mM,
        sys.nadh_m => 0.7273229914481529mM,
        sys.dpsi => 146.53914595776085mV,
        sys.sox_m => 7.226067750253397e-9mM,
        sys.vm => -87.08987549736572mV,
    ]
end

function build_model(; name, use_mg=false, simplify=true, bcl=1second, istim=-80μA / cm², tstart=0second, tend=10second, duty=0.5ms)
    @parameters begin
        iStim(t) = 0μA / cm²    # Stimulation current
        DOX(t) = 0mM               # Doxorubicin concentration
        MT_PROT = 1             # OXPHOS protein content
        ΣA_m = 1.01mM           # Mitochondrial ATP + ADP pool (Gauthier-2013)
        ΣA_i = 8mM              # Cytosolic ATP + ADP pool (Li-2015)
        ΣNAD_m = 10mM           # Mitochondrial NAD + NADH pool # 1.0mM (Gauthier-2013)
        ΣNADP_m = 0.1mM         # Mitochondrial NADP + NADPH pool (Gauthier-2013)

        CM = 1μF / cm²          # Plasma membrane capacitance
        CM_MITO = 1.812μM / mV  # Inner mitochondrial membrane capacitance

        # Cell geometries
        A_CAP = 1.534E-4cm²
        V_MYO = 25.84pL
        V_MT = 15.89pL
        V_NSR = 1.4pL
        V_JSR = 0.16pL
        V_SS = 0.495E-3pL
        V_MITO_V_MYO = V_MT / V_MYO         # Ratio of mitochondrial matrix to cytosolic ion diffusion space
        # conversion factor from current (μA/cm²) to traslocation rate (mM/ms)
        A_CAP_V_MYO_F = A_CAP / (V_MYO * Faraday)
        A_CAP_V_SS_F = A_CAP / (V_SS * Faraday)

        # Calcium buffering
        KM_CA_CMDN = 2.38μM
        ET_CMDN = 50μM
        KM_CA_CSQN = 0.8mM
        ET_CSQN = 5.0mM
        δCA = 3E-4          # Mito. Ca buffering factor
        # concentrations
        nadph_i = 0.075mM   # Cytosolic NADPH (Gauthier-2013) # 1.0mM (Li-2015)
        pi_m = 8.6512mM     # Inorganic phosphate (Gauthier-2013 and Kembro-2013)
        k_o = 5.4mM         # Extracellular potassium
        na_o = 140.0mM      # Extracellular sodium
        ca_o = 2.0mM        # Extracellular calcium
        mg_m = 0.4mM        # Mitochondrial magnesium (li-2015)
        mg_i = 1.0mM        # Cytosolic magnesium (Gauthier-2013) # 3.1mM (Li-2015)
        pi_i = 3.0mM        # Cytosolic inorganic phosphate
        h_i = 10^(-7) * Molar   # Cytosolic proton concentration
        h_m = 10^(-7.6) * Molar # Mitochondrial proton concentration
    end

    @variables begin
        vm(t) = -87.28mV        # Sarcolemmal membrane potential
        dpsi(t) = 175.318mV     # Mitochondrial membrane potential
        atp_i(t) # Conserved
        adp_i(t) = 0.13mM
        atp_m(t) # Conserved
        adp_m(t) = 0.0058mM
        nad_m(t) # Conserved
        nadh_m(t) = 9.26mM
        na_i(t) = 8.2143mM      # Cytoplasmic Na
        k_i(t) = 150.8mM        # Cytoplasmic K
        ca_i(t) = 5.8653e-5mM   # Cytoplasmic Ca
        ca_jsr(t) = 0.1948mM    # Junctional SR Ca
        ca_nsr(t) = 0.1948mM    # Network SR Ca
        ca_ss(t) = 7.057e-5mM   # Subspace calcium
        ca_m(t) = 2.19e-5mM     # Mitochondrial calcium
        sox_m(t) = 0.5μM        # Mitochondrial superoxide
    end

    # Subsystems
    lccsys = get_lcc_sys(ca_ss, ca_o, k_i, k_o, vm)
    ryrsys = get_ryr_sys(ca_jsr, ca_ss)
    cksys = get_ck_sys(atp_i, adp_i)
    forcesys = get_force_sys(atp_i, adp_i, ca_i)
    jcasys = get_jca_sys(atp_i, adp_i, ca_i, ca_nsr, ca_jsr, ca_ss, ca_o, na_i, na_o, vm)
    iksys = get_ik_sys(na_i, na_o, k_i, k_o, mg_i, vm, atp_i, adp_i)
    inasys = get_ina_sys(na_i, na_o, ca_i, ca_o, vm)
    inaksys = get_inak_sys(atp_i, adp_i, vm, na_i, na_o, k_o)
    tcassys = get_tca_sys(atp_m, adp_m, nad_m, nadh_m, h_m, ca_m, mg_m; use_mg)
    @unpack suc, fum, oaa, vSL, vIDH, vKGDH, vMDH = tcassys
    etcsys = get_etc_sys(nad_m, nadh_m, dpsi, h_i, h_m, sox_m, suc, fum, oaa, DOX, MT_PROT)
    c5sys = get_c5_sys(dpsi, h_i, h_m, atp_i, adp_i, atp_m, adp_m, pi_m, MT_PROT; mg_i, mg_m, use_mg)
    @unpack vROS, vC1, vHres = etcsys
    rossys = get_ros_sys(dpsi, sox_m, nadph_i, V_MITO_V_MYO)
    mitocasys = get_mitoca_sys(na_i, ca_m, ca_i, dpsi)

    @unpack vTrROS, vIMAC = rossys
    @unpack INa, INsNa, ICaB, INaB = inasys
    @unpack IPMCA, Jup, INaCa, Jtr, Jxfer = jcasys
    @unpack IK, IK1, IKp, IKatp = iksys
    @unpack INaK = inaksys
    @unpack ICaK, ICaL = lccsys
    @unpack vCK_mito = cksys
    @unpack vANT, vC5, vHu, vHleak = c5sys
    @unpack Jrel = ryrsys
    @unpack vAm, Jtrpn = forcesys
    @unpack vUni, vNaCa = mitocasys

    eqs = [
        D(vm) * (-CM) ~ INa + ICaL + IK + IK1 + IKp + INaCa + INaK + INsNa + IPMCA + ICaB + INaB + iStim + IKatp + ICaK,
        D(na_i) ~ -A_CAP_V_MYO_F * (INa + INaB + INsNa + 3 * (INaCa + INaK)),
        D(k_i) ~ -A_CAP_V_MYO_F * (IK + IK1 + IKp + IKatp + iStim - 2 * INaK + ICaK),
        D(ca_i) ~ β_ca(ca_i, KM_CA_CMDN, ET_CMDN) * (Jxfer - Jup - Jtrpn - 0.5 * A_CAP_V_MYO_F * (IPMCA + ICaB - 2 * INaCa) + V_MITO_V_MYO * (vNaCa - vUni)),
        D(ca_nsr) ~ β_ca(ca_nsr, KM_CA_CSQN, ET_CSQN) * (V_MYO * Jup - V_JSR * Jtr) / V_NSR,
        D(ca_jsr) ~ β_ca(ca_jsr, KM_CA_CSQN, ET_CSQN) * (Jtr - Jrel),
        D(ca_ss) ~ β_ca(ca_ss, KM_CA_CMDN, ET_CMDN) * ((V_JSR * Jrel - V_MYO * Jxfer) / V_SS - 0.5 * ICaL * A_CAP_V_SS_F),
        D(ca_m) ~ δCA * (vUni - vNaCa),
        D(adp_i) ~ vCK_mito - V_MITO_V_MYO * vANT + vAm + 0.5 * Jup + A_CAP_V_MYO_F * (IPMCA + INaK),
        D(adp_m) ~ vANT - vSL - vC5,
        D(nadh_m) ~ -vC1 + vIDH + vKGDH + vMDH,
        D(dpsi) ~ inv(CM_MITO) * (vHres - vHu - vANT - vHleak - vNaCa - 2 * vUni - vTrROS),
        D(sox_m) ~ vROS - vTrROS,
        ΣA_i ~ atp_i + adp_i,
        ΣA_m ~ atp_m + adp_m,
        ΣNAD_m ~ nad_m + nadh_m,
    ]

    if bcl > 0
        discrete_events = [collect(tstart:bcl:tend) => [iStim ~ istim], collect(tstart+duty:bcl:tend) => [iStim ~ 0]]
        sys = ODESystem(eqs, t; name, discrete_events)
    else
        sys = ODESystem(eqs, t; name)
    end

    for s in (lccsys, ryrsys, cksys, forcesys, jcasys, iksys, inasys, inaksys, tcassys, etcsys, c5sys, rossys, mitocasys)
        sys = extend(sys, s)
    end

    if simplify
        sys = structural_simplify(sys)
    end

    return sys
end # build_model()

end # Module
