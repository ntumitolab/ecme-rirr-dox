module ECMEDox

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using NaNMath
using LogExpFunctions

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
        sys.sox_i => 1.0222129545264641e-6,
        sys.h2o2_i => 7.55888231084137e-7,
        sys.gssg_i => 0.002220885967007533,
        sys.Q_n => 1.7555277577733173,
        sys.Qdot_n => 0.250155168931016,
        sys.QH2_n => 0.10775990996722686,
        sys.QH2_p => 0.10771653709539612,
        sys.Qdot_p => 0.023269397495713263,
        sys.b1 => 0.19471601197446434,
        sys.b2 => 0.08033504811661059,
        sys.b3 => 0.047309996276701564,
        sys.fes_ox => 0.20066008052559806,
        sys.cytc1_ox => 0.28338319859718647,
        sys.cytc_ox => 0.24903758981942076,
        sys.isoc => 0.05159054318098895,
        sys.akg => 0.051145197718677655,
        sys.scoa => 0.03508849487000582,
        sys.suc => 0.0019107469302081612,
        sys.fum => 0.1751906841603877,
        sys.mal => 0.15856757152954906,
        sys.oaa => 0.011576938766421891,
        sys.m_na => 0.0010899662356754164,
        sys.h_na => 0.9905304527574307,
        sys.j_na => 0.9892354695981638,
        sys.x_k => 0.0012736748778708119,
        sys.ltr_ca => 0.02049022737620247,
        sys.htr_ca => 0.1372356208938611,
        sys.p0 => 0.0036716635421201937,
        sys.p1 => 0.0033538440069086255,
        sys.p2 => 0.006345552748082013,
        sys.p3 => 0.005553300147495315,
        sys.n1 => 0.012559971511989917,
        sys.adp_ic => 0.34296593237569856,
        sys.crp_i => 4.3883068168196715,
        sys.crp_ic => 4.376369356028448,
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
        sys.y_ca_lcc => 0.737396088676386,
        sys.na_i => 15.360882167675042,
        sys.k_i => 142.39854398555977,
        sys.ca_i => 0.00016134091543239482,
        sys.ca_nsr => 0.9944103502674689,
        sys.ca_jsr => 0.9760651184113716,
        sys.ca_ss => 0.000169823015720396,
        sys.ca_m => 0.0002427548693370271,
        sys.adp_i => 0.19176238251620303,
        sys.adp_m => 0.12562731414699957,
        sys.nadh_m => 0.7273229914481529,
        sys.dpsi => 0.14653914595776085,
        sys.sox_m => 7.226067750253397e-9,
        sys.vm => -0.08708987549736572,
    ]
end

function build_model(;name, use_mg=false, simplify=true, bcl=1second, istim=-80μA/cm², tstart=0second, tend = 10second, duty=0.5ms)
    @parameters begin
        iStim = 0μA/cm²     # Stimulation current
        DOX = 0mM           # Doxorubicin concentration
        MT_PROT = 1         # OXPHOS protein content
        ΣA_m = 1.01mM       # Mitochondrial ATP + ADP pool (Gauthier-2013)
        ΣA_i = 8mM          # Cytosolic ATP + ADP pool (Li-2015)
        ΣNAD_m = 10mM       # Mitochondrial NAD + NADH pool # 1.0mM (Gauthier-2013)
        ΣNADP_m = 0.1mM     # Mitochondrial NADP + NADPH pool (Gauthier-2013)

        CM = 1μF/cm²                # Plasma membrane capacitance
        CM_MITO = 1.812E-3mM / mV   # Inner mitochondrial membrane capacitance

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
        KM_CA_CMDN=2.38μM
        ET_CMDN=50μM
        KM_CA_CSQN=0.8mM
        ET_CSQN=5.0mM
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
        sox_m(t) = 0.5μM     # Mitochondrial superoxide
    end

    D = Differential(t)
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
    @unpack iNa, iNsNa, iCaB, iNaB = inasys
    @unpack iPMCA, jUp, iNaCa, jTr, jXfer = jcasys
    @unpack iK, iK1, iKp, iKatp = iksys
    @unpack iNaK = inaksys
    @unpack iCaK, iCaL = lccsys
    @unpack vck_mito = cksys
    @unpack vANT, vATPase, vHu, vHleak = c5sys
    @unpack jRel = ryrsys
    @unpack vAM, jTrpn = forcesys
    @unpack vUni, vNaCa = mitocasys

    eqs = [
        D(vm) * (-CM) ~ iNa + iCaL + iK + iK1 + iKp + iNaCa + iNaK + iNsNa + iPMCA + iCaB + iNaB + iStim + iKatp + iCaK,
        D(na_i) ~ -A_CAP_V_MYO_F * (iNa + iNaB + iNsNa + 3 * (iNaCa + iNaK)),
        D(k_i) ~ -A_CAP_V_MYO_F * (iK + iK1 + iKp + iKatp + iStim - 2 * iNaK + iCaK),
        D(ca_i) ~ β_ca(ca_i, KM_CA_CMDN, ET_CMDN) * (jXfer - jUp - jTrpn - 0.5 * A_CAP_V_MYO_F * (iPMCA + iCaB - 2 * iNaCa) + V_MITO_V_MYO * (vNaCa - vUni)),
        D(ca_nsr) ~ β_ca(ca_nsr, KM_CA_CSQN, ET_CSQN) * (V_MYO * jUp - V_JSR * jTr) / V_NSR,
        D(ca_jsr) ~ β_ca(ca_jsr, KM_CA_CSQN, ET_CSQN) * (jTr - jRel),
        D(ca_ss) ~ β_ca(ca_ss, KM_CA_CMDN, ET_CMDN) * ( (V_JSR * jRel - V_MYO * jXfer)/V_SS - 0.5 * iCaL * A_CAP_V_SS_F),
        D(ca_m) ~ δCA * (vUni - vNaCa),
        D(adp_i) ~ vck_mito - V_MITO_V_MYO * vANT + vAM + 0.5 * jUp + A_CAP_V_MYO_F * (iPMCA + iNaK),
        D(adp_m) ~ vANT - vSL - vATPase,
        D(nadh_m) ~ -vC1 + vIDH + vKGDH + vMDH,
        D(dpsi) ~ inv(CM_MITO) * (vHres - vHu - vANT - vHleak - vNaCa - 2 * vUni - vTrROS),
        D(sox_m) ~ vROS - vTrROS,
        ΣA_i ~ atp_i + adp_i,
        ΣA_m ~ atp_m + adp_m,
        ΣNAD_m ~ nad_m + nadh_m,
    ]

    if bcl > 0
        discrete_events = [collect(tstart:bcl:tend) => [iStim ~ istim], collect(tstart+duty:bcl:tend)=> [iStim ~ 0]]
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
