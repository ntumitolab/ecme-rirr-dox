using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using NaNMath
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms
using Plots

# Adapted from Markevich, 2015
function c1_markevich(; name=:c1sys,
    Q_n=3.0mM, QH2_n=0.3mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)
    @parameters begin
        Em_O2_SOX = -160mV        # O2/Superoxide redox potential
        Em_FMNsq_FMNH = -375mV    # FMN semiquinone/FMNH- redox potential
        Em_NAD = -320mV           # NAD/NADH redox potential
        Em_FMN_FMNH = -340mV      # FMN/FMNH- redox potential
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -400mV       # -213mV
        Em_SQ_QH2_C1 = +600mV
        ET_C1 = 0.4mM  # Activity of complex I
        KI_NAD_C1 = 1mM
        KI_NADH_C1 = 50μM
        KD_NAD_C1 = 25μM
        KD_NADH_C1 = 100μM
        KI_DOX_C1 = 400μM  # DOX inhibition concentration (IC50) on complex I
        K1_C1 = 10Hz / μM
        KEQ1_C1 = 0.1 / μM     # Association constant for Q
        K2_C1 = 4E5Hz / μM
        KEQ2_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        K3_C1 = 2.7E6Hz / μM
        K4_C1 = 1000Hz
        KEQ4_C1 = 20μM
        K5_C1 = 2Hz / μM
        KEQ5_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        K6_C1 = 0.04Hz / μM
        KEQ6_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        I_C1(t) # Conserved
        Q_C1(t) = 0
        SQ_C1(t) = 0
        QH2_C1(t) = 0
        fFMNH_C1(t)  # fraction of fully-reduced FMN
        fFMNsq_C1(t) # fraction of FMN semiquinone
        N2_C1(t)
        N2m_C1(t)
        KEQ3_C1(t)
        vQ_C1(t)
        vQH2_C1(t)
        vNADH_C1(t)
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
        TN_C1(t) # Turnover number
        ROSshunt_C1(t)
    end

    # Q association
    v1 = K1_C1 * (I_C1 * Q_n - Q_C1 / KEQ1_C1)
    # First electron transfer
    v2 = K2_C1 * (Q_C1 * N2m_C1 - SQ_C1 * N2_C1 / KEQ2_C1)
    # Second electron transfer
    fhm = h_m / 1E-7Molar
    v3 = K3_C1 * (SQ_C1 * N2m_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ3_C1)
    # QH2 dissociation
    v4 = K4_C1 * (QH2_C1 - I_C1 * QH2_n / KEQ4_C1)
    # Flavin site ROS generation
    v5 = K5_C1 * ET_C1 * (fFMNH_C1 * O2 - fFMNsq_C1 * sox_m / KEQ5_C1)
    # Quinone site ROS generation
    v6 = K6_C1 * (SQ_C1 * O2 - Q_C1 * sox_m / KEQ6_C1)

    rN2 = exp(iVT * (Em_N2 - Em_NAD)) * (nadh / nad)

    eqs = [
        fFMNH_C1 ~ inv(1 + nad / KD_NAD_C1 + nadh / KI_NADH_C1 + exp(2iVT * (Em_FMN_FMNH - Em_NAD) * (nad / nadh) * (1 + nad / KI_NAD_C1 + nadh / KD_NADH_C1))),
        fFMNsq_C1 ~ exp(iVT * (Em_FMNsq_FMNH - Em_FMN_FMNH)) * fFMNH_C1,
        KEQ3_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_m / h_i)^4,
        N2m_C1 ~ ET_C1 * rN2 / (1 + rN2),
        ET_C1 ~ N2m_C1 + N2_C1,
        # Positive: production; negative: consumption
        vQ_C1 ~ -v1,
        vQH2_C1 ~ v4,
        vROSIf ~ v5,
        vROSIq ~ v6,
        vROS_C1 ~ vROSIf + vROSIq,
        vNADH_C1 ~ -0.5 * (v2 + v3 + v5),
        TN_C1 ~ -vNADH_C1 / ET_C1,
        ROSshunt_C1 ~ 0.5vROS_C1 / (-vNADH_C1-vQH2_C1+vQ_C1),
        ET_C1 ~ I_C1 + Q_C1 + SQ_C1 + QH2_C1,
        D(Q_C1) ~ v1 - v2 + v6,
        D(SQ_C1) ~ v2 - v3 - v6,
        D(QH2_C1) ~ v3 - v4,
    ]
    return ODESystem(eqs, t; name)
end

function c1_gauthier(; name=:c1sys,
    Q_n=3.0mM, QH2_n=0.3mM, nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    C1_INHIB = 1)
    @parameters begin
        ET_C1 = 8.85mM      # Activity of complex I
        dpsi_B_C1 = 50mV   # Phase boundary potential
        # Transition rates
        K12_C1 = 6339.6Hz # pH = 7
        K21_C1 = 5Hz
        K56_C1 = 100Hz
        K65_C1 = 251190Hz # pH = 7
        K61_C1 = 1e7Hz
        K16_C1 = 130Hz
        K23_C1 = 3886.7Hz / sqrt(mM)
        K32_C1 = 9.1295e6Hz
        K34_C1 = 639.1364Hz
        K43_C1 = 3.2882Hz / sqrt(mM)
        K47_C1 = 1.5962E7Hz/mM
        K74_C1 = 65.2227Hz
        K75_C1 = 24.615E3Hz
        K57_C1 = 1.1667E3Hz / sqrt(mM)
        K42_C1 = 6.0318Hz / mM
        Em_O2_SOX = -160mV         # O2/Superoxide redox potential
        Em_FMNH2_FMNH = -375mV     # FMNH/FMNH2 redox potential
    end

    @variables begin
        C1_1(t) = 0
        C1_2(t) # Conserved
        C1_3(t) = 0
        C1_4(t) = 0
        C1_5(t) = 0
        C1_6(t) = 0
        C1_7(t) = 0
        vQC1(t)
        vNADHC1(t)
        vROSC1(t)
        TN_C1(t) # Turnover number
        ROSshunt_C1(t)
    end

    fhi = h_i / 1E-7Molar
    fhm = h_m / 1E-7Molar
    fv = exp(iVT * (dpsi - dpsi_B_C1))
    # State transition rates
    a12 = K12_C1 * fhm^2
    a21 = K21_C1
    a65 = K65_C1 * fhi^2
    a56 = K56_C1
    a61 = K61_C1 / fv
    a16 = K16_C1 * fv
    a23 = K23_C1 * NaNMath.sqrt(nadh)
    a32 = K32_C1
    a34 = K34_C1
    a43 = K43_C1 * NaNMath.sqrt(nad)
    a47 = C1_INHIB * K47_C1 * NaNMath.sqrt(Q_n * h_m)
    a74 = K74_C1
    a57 = C1_INHIB * K57_C1 * NaNMath.sqrt(QH2_n)
    a75 = K75_C1
    a42 = K42_C1 * O2
    a24 = K42_C1 * exp(iVT * (Em_FMNH2_FMNH - Em_O2_SOX)) * sox_m

    v12 = a12 * C1_1 - a21 * C1_2
    v23 = a23 * C1_2 - a32 * C1_3
    v34 = a34 * C1_3 - a43 * C1_4
    v47 = a47 * C1_4 - a74 * C1_7
    v75 = a75 * C1_7 - a57 * C1_5
    v56 = a56 * C1_5 - a65 * C1_6
    v61 = a61 * C1_6 - a16 * C1_1
    v42 = a42 * C1_4 - a24 * C1_2

    eqs = [
        ET_C1 ~ C1_1+C1_2+C1_3+C1_4+C1_5+C1_6+C1_7,
        D(C1_1) ~ v61 - v12,
        D(C1_3) ~ v23 - v34,
        D(C1_4) ~ v34 - v47 - v42,
        D(C1_7) ~ v47 - v75,
        D(C1_5) ~ v75 - v56,
        D(C1_6) ~ v56 - v61,
        vQC1 ~ -0.5v47,
        vROSC1 ~ v42,
        vNADHC1 ~ -v23,
        TN_C1 ~ -vNADHC1 / ET_C1,
        ROSshunt_C1 ~ -vROSC1 / vNADHC1,
    ]
    return ODESystem(eqs, t; name)
end

@parameters begin
    Q_n = 3.0mM
    QH2_n = 0.3mM
    nad = 500μM
    nadh = 500μM
    dpsi = 150mV
end

markevich = c1_markevich(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
gauthier = c1_gauthier(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify

dpsirange = 100mV:10mV:200mV
alter_dpsi = (prob, i, repeat) -> remake(prob, p=[dpsi => dpsirange[i]])

tend = 1000.0ms
prob_m = SteadyStateProblem(markevich, [])
prob_g = SteadyStateProblem(gauthier, [])
alg = DynamicSS(Rodas5P())
# Try
@time sol_m = solve(prob_m, alg)
@time sol_g = solve(prob_g, alg)

# Change in MMP
dpsirange = 100mV:5mV:200mV
alter_dpsi = (prob, i, repeat) -> remake(prob, p=[dpsi => dpsirange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func = alter_dpsi)
eprob_m = EnsembleProblem(prob_m; prob_func = alter_dpsi)
@time sim_g = solve(eprob_g, alg, EnsembleSerial(); trajectories=length(dpsirange))
@time sim_m = solve(eprob_m, alg, EnsembleSerial(); trajectories=length(dpsirange))

# MMP vs NADH turnover
xs = dpsirange
ys_g = map(sim_g) do sol
    sol[gauthier.TN_C1] * 1000
end
ys_m = map(sim_m) do sol
    sol[markevich.TN_C1]
end
plot(xs, [ys_g ys_m], xlabel="MMP (mV)", ylabel="NADH turnover", label=["Gauthier (Hz)" "Markevich (kHz)"])

# MMP vs ROS production
xs = dpsirange
ys_g = map(sim_g) do sol
    sol[gauthier.ROSshunt_C1]
end
ys_m = map(sim_m) do sol
    sol[markevich.ROSshunt_C1] |> abs
end
plot(xs, [ys_g ys_m], xlabel="MMP (mV)", ylabel="ROS production", label=["Gauthier" "Markevich"])

ys_if = map(sim_m) do sol
    sol[markevich.vROSIf] |> abs
end
ys_iq = map(sim_m) do sol
    sol[markevich.vROSIq] |> abs
end
plot(xs, [ys_if ys_iq], xlabel="MMP (mV)", ylabel="ROS production", label=["IF" "IQ"])

sol_m[markevich.TN_C1]
sol_g[gauthier.TN_C1]

sol_m[markevich.vROSIf]
sol_m[markevich.vROSIq]
sol_m[markevich.ROSshunt_C1]
sol_g[gauthier.ROSshunt_C1]
