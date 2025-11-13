module ECMEDox

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using NaNMath

export build_model, build_u0, build_stim_callbacks

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
include("u0.jl")

function build_model(; name=:ecmesys, bcl=0second, istim=-80μAcm⁻², tstart=0second, tend=10second, duty=0.5ms, use_mg=false, simplify=true)
    @parameters begin
        iStim(t) = 0μAcm⁻²          # Stimulation current
        DOX = 0mM                # Doxorubicin concentration
        kdiffO2 = 1000Hz            # Oxygen diffusion rate
        MT_PROT = 1                 # OXPHOS protein content
        ΣA_m = 1.01mM               # Mitochondrial ATP + ADP pool (Gauthier-2013)
        ΣA_i = 8mM                  # Cytosolic ATP + ADP pool (Li-2015)
        ΣNAD_m = 3mM                # Mitochondrial NAD + NADH pool # 1.0mM (Gauthier-2013)
        ΣNADP_m = 0.1mM             # Mitochondrial NADP + NADPH pool (Gauthier-2013)
        iCM = inv(1μFcm⁻²)          # inverse of Plasma membrane capacitance
        iCMito = inv(1.812μM / mV)  # inverse of inner mitochondrial membrane capacitance

        # Cell geometries
        A_CAP = 1.534E-4cm²         # capacitance area
        V_MYO = 25.84pL             # Cytosol volume
        V_MT = 15.89pL              # Mitochondrial (matrix) volume
        V_NSR = 1.4pL               # Network SR volume
        V_JSR = 0.16pL              # Junctional SR volume
        V_SS = 0.495E-3pL           # Dyadic space volume
        V_MITO_V_MYO = V_MT / V_MYO         # Ratio of mitochondrial matrix to cytosolic ion diffusion space
        # Conversion factor from current density (μA/cm²) to ion flux (μM/ms)
        A_CAP_V_MYO_F = A_CAP / (V_MYO * Faraday)
        A_CAP_V_SS_F = A_CAP / (V_SS * Faraday)

        # Calcium buffering
        KM_CA_CMDN = 2.38μM
        ET_CMDN = 50μM
        KM_CA_CSQN = 800μM
        ET_CSQN = 5mM
        δCA = 3E-4          # Mito. Ca buffering factor

        # constant concentrations
        h_i = exp10(-7) * Molar   # Cytosolic proton concentration
        h_m = exp10(-7.6) * Molar # Mitochondrial proton concentration
        nadph_i = 75μM      # Cytosolic NADPH (Gauthier-2013) # 1.0mM (Li-2015)
        pi_m = 8.6512mM     # Inorganic phosphate (Gauthier-2013 and Kembro-2013)
        k_o = 5.4mM         # Extracellular potassium
        na_o = 140mM        # Extracellular sodium
        ca_o = 2mM          # Extracellular calcium
        mg_i = 1mM          # Cytosolic magnesium (Gauthier-2013) # 3.1mM (Li-2015)
        mg_m = 0.4mM        # Mitochondrial magnesium (li-2015)
        pi_i = 3mM          # Cytosolic inorganic phosphate
        O2_o = 6μM       # Extracellular oxygen
    end

    @variables begin
        vm(t) = -87.28mV        # Sarcolemmal membrane potential
        dpsi(t) = 170mV         # Mitochondrial membrane potential
        atp_i(t) # Conserved
        adp_i(t) = 50μM
        atp_m(t) # Conserved
        adp_m(t) = 50μM
        nad_m(t) # Conserved
        nadh_m(t) = 500μM
        na_i(t) = 8.2143mM      # Cytoplasmic Na
        k_i(t) = 150.8mM        # Cytoplasmic K
        ca_i(t) = 5.8653e-5mM   # Cytoplasmic Ca
        ca_jsr(t) = 0.1948mM    # Junctional SR Ca
        ca_nsr(t) = 0.1948mM    # Network SR Ca
        ca_ss(t) = 7.057e-5mM   # Subspace calcium
        ca_m(t) = 2.19e-5mM     # Mitochondrial calcium
        sox_m(t) = 0.1μM        # Mitochondrial superoxide
        O2(t) = 6μM             # Mitochondrial oxygen concentration
    end

    # Subsystems
    @unpack eqs_cicr16, ICaL, ICaK, Jrel, Jtr, Jxfer = get_cicr16_eqs(; vm, ca_ss, ca_i, ca_jsr, ca_nsr, k_i, ca_o, k_o)
    @unpack eqs_ck, vCK_mito, vCK_cyto, vTR_crp = get_ck_eqs(; atp_i, adp_i)
    @unpack eqs_force, vAm, Jtrpn = get_force_eqs(; atp_i, adp_i, ca_i)
    @unpack eqs_jca, IPMCA, Jup, INaCa = get_jca_eqs(; atp_i, adp_i, ca_i, ca_nsr, ca_o, na_i, na_o, vm)
    @unpack eqs_ik, IK1, IK, IKp, IKatp = get_ik_eqs(; na_i, na_o, k_i, k_o, vm, atp_i, adp_i, mg_i)
    @unpack eqs_ina, INa, INsNa, INaB, ICaB = get_ina_eqs(; na_i, na_o, ca_i, ca_o, vm)
    @unpack eqs_nak, INaK = get_inak_eqs(; atp_i, adp_i, vm, na_i, na_o, k_o)
    @unpack eqs_tca, vIDH, vKGDH, vMDH, vSL, suc, fum, oaa = get_tca_eqs(; atp_m, adp_m, nad_m, nadh_m, ca_m, h_m, pi_m, mg_m, use_mg)
    @unpack eqs_etc, vHres, vROS, vSDH, vNADHC1, vO2 = get_etc_eqs(; nad_m, nadh_m, dpsi, sox_m, suc, fum, oaa, DOX, MT_PROT, O2, h_i, h_m)
    @unpack eqs_c5, vANT, vC5, vHu, vHleak = get_c5_eqs(; dpsi, h_i, h_m, atp_i, adp_i, atp_m, adp_m, pi_m, MT_PROT, mg_i, mg_m)
    @unpack eqs_ros, vTrROS, vIMAC, vCAT, vGPX_i, vGR_i, vSOD_i = get_ros_eqs(; dpsi, sox_m, nadph_i, V_MITO_V_MYO)
    @unpack eqs_mitoca, vUni, vNaCa = get_mitoca_eqs(; na_i, ca_m, ca_i, dpsi)

    eqs = [
        ΣA_i ~ atp_i + adp_i,
        ΣA_m ~ atp_m + adp_m,
        ΣNAD_m ~ nad_m + nadh_m,
        D(vm) ~ -iCM * (INa + ICaL + IK + IK1 + IKp + INaCa + INaK + INsNa + IPMCA + ICaB + INaB + iStim + IKatp + ICaK),
        D(na_i) ~ -A_CAP_V_MYO_F * (INa + INaB + INsNa + 3INaCa + 3INaK),
        D(k_i) ~ -A_CAP_V_MYO_F * (IK + IK1 + IKp + IKatp + iStim - 2INaK + ICaK),
        D(ca_i) ~ β_ca(ca_i, KM_CA_CMDN, ET_CMDN) * (Jxfer - Jup - Jtrpn - 0.5 * A_CAP_V_MYO_F * (IPMCA + ICaB - 2INaCa) + V_MITO_V_MYO * (vNaCa - vUni)),
        D(ca_nsr) ~ (V_MYO * Jup - V_JSR * Jtr) / V_NSR,
        D(ca_jsr) ~ β_ca(ca_jsr, KM_CA_CSQN, ET_CSQN) * (Jtr - Jrel),
        D(ca_ss) ~ β_ca(ca_ss, KM_CA_CMDN, ET_CMDN) * ((V_JSR * Jrel - V_MYO * Jxfer) / V_SS - 0.5ICaL * A_CAP_V_SS_F),
        D(ca_m) ~ δCA * (vUni - vNaCa),
        D(adp_i) ~ -V_MITO_V_MYO * vANT + vCK_mito + vAm + 0.5Jup + A_CAP_V_MYO_F * (IPMCA + INaK),
        D(adp_m) ~ vANT - vSL - vC5,
        D(nadh_m) ~ vNADHC1 + vIDH + vKGDH + vMDH,
        D(dpsi) ~ iCMito * (vHres - vHu - vANT - vHleak - vNaCa - 2vUni - vIMAC),
        D(sox_m) ~ vROS - vTrROS,
        D(O2) ~ V_MT / (V_MYO + V_MT) * (-vROS - vO2) + kdiffO2 * (O2_o - O2),
    ]

    eqs = [eqs; eqs_cicr16; eqs_ck; eqs_force; eqs_jca; eqs_ik; eqs_ina; eqs_nak; eqs_tca; eqs_etc; eqs_c5; eqs_ros; eqs_mitoca]

    if bcl > 0 && duty > 0 && istim != 0
        ts0 = collect(tstart:bcl:tend)
        ts1 = ts0 .+ duty

        sstart = ModelingToolkit.SymbolicDiscreteCallback( ts0 => [iStim ~ istim]; discrete_parameters = [iStim], iv = t)
        send = ModelingToolkit.SymbolicDiscreteCallback( ts1 => [iStim ~ 0]; discrete_parameters = [iStim], iv = t)
        sys = System(eqs, t; name, discrete_events = [sstart, send])
    else
        sys = System(eqs, t; name)
    end

    if simplify
        sys = mtkcompile(sys)
    end

    return sys
end # build_model()

end # module ECMEDox
