# # Complex I model
# Comparing Gauthier, Markevich, and simplified complex I models
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using NaNMath
using Plots
using ECMEDox
using ECMEDox: mM, μM, nM, iVT, mV, Molar, Hz, ms

# Gauthier 2012 7-state QSSA model
function c1_gauthier(; name=:c1gauthier,
    Q_n=3.0mM, QH2_n=0.3mM, nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    C1_INHIB=1)
    @parameters begin
        ET_C1 = 8.85mM      ## Activity of complex I
        dpsi_B_C1 = 50mV    ## Phase boundary potential
        ## Transition rates
        K12_C1 = 6339.6Hz ## pH = 7
        K21_C1 = 5Hz
        K56_C1 = 100Hz
        K65_C1 = 251190Hz ## pH = 7
        K61_C1 = 1e7Hz
        K16_C1 = 130Hz
        K23_C1 = 3886.7Hz / sqrt(mM)
        K32_C1 = 9.1295e6Hz
        K34_C1 = 639.1364Hz
        K43_C1 = 3.2882Hz / sqrt(mM)
        K47_C1 = 1.5962E7Hz / mM
        K74_C1 = 65.2227Hz
        K75_C1 = 24.615E3Hz
        K57_C1 = 1.1667E3Hz / sqrt(mM)
        K42_C1 = 6.0318Hz / mM
        Em_O2_SOX = -160mV         ## O2/Superoxide redox potential
        Em_FMNH2_FMNH = -375mV     ## FMNH/FMNH2 redox potential
    end

    @variables begin
        C1_1(t) = 0
        C1_2(t) ## Conserved
        C1_3(t) = 0
        C1_4(t) = 0
        C1_5(t) = 0
        C1_6(t) = 0
        C1_7(t) = 0
        vQ_C1(t)
        vNADH_C1(t)
        vROS_C1(t)
        TN_C1(t) ## Turnover number
    end

    fhi = h_i / 1E-7Molar
    fhm = h_m / 1E-7Molar
    fv = exp(iVT * (dpsi - dpsi_B_C1))
    ## State transition rates
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
        ET_C1 ~ C1_1 + C1_2 + C1_3 + C1_4 + C1_5 + C1_6 + C1_7,
        D(C1_1) ~ v61 - v12,
        D(C1_3) ~ v23 - v34,
        D(C1_4) ~ v34 - v47 - v42,
        D(C1_7) ~ v47 - v75,
        D(C1_5) ~ v75 - v56,
        D(C1_6) ~ v56 - v61,
        vQ_C1 ~ -0.5v47,
        vROS_C1 ~ v42,
        vNADH_C1 ~ -0.5v23,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

# Markevich 2015 mass action model
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4426091/
function c1_markevich_full(; name=:c1markevich_full,
    Q_n=1.8mM, QH2_n=0.2mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)
    @parameters begin
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -80mV
        Em_N1a = -370mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## ~800mV (?) in Markevich, 2015
        ET_C1 = 17μM              ## Activity of complex I
        ## DOX IC50 on complex I
        KI_DOX_C1 = 400μM
        kf1_C1 = 83Hz / μM
        KEQ1_C1 = 0.01 / μM
        kf2_C1 = 1.44e12Hz
        kf3_C1 = 1e6Hz
        KEQ3_C1 = 25μM
        KEQ2_C1 = exp(2iVT * (Em_FMN_FMNH - Em_NAD)) / KEQ1_C1 / KEQ3_C1
        kf4_C1 = 1Hz / μM
        KEQ4_C1 = 0.001 / μM
        kf5_C1 = 2Hz / μM
        KEQ5_C1 = 0.02 / μM
        kf6_C1 = 5e8Hz / μM
        KEQ6_C1 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
        kf7_C1 = 10000Hz / μM
        KEQ7_C1 = exp(iVT * (Em_N2 - Em_N3))
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = 0.1 / μM         ## Association constant for Q
        kf9_C1 = 4E5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kf10_C1 = 2e6Hz / μM
        KEQ10_C1 = exp(iVT * (Em_N1a - Em_FMN_FMNsq))
        KEQ10B_C1 = exp(iVT * (Em_FMNsq_FMNH - Em_N1a))
        kf11_C1 = 1e9Hz / μM
        KEQ11_C1 = exp(iVT * (Em_N3 - Em_FMN_FMNsq))
        kf13_C1 = 2.7e6Hz / μM
        kf14_C1 = 1000Hz
        KEQ14_C1 = 20μM             ## Dissociation constant for QH2
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        KEQ17_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        Iq_C1(t) ## Conserved
        Q_C1(t) = 0
        SQ_C1(t) = 0
        QH2_C1(t) = 0
        FMN(t) ## Conserved
        FMN_NADH(t) = 0
        FMNH_NAD(t) = 0
        FMN_NAD(t) = 0
        FMNH_NADH(t) = 0
        FMNH(t) = 0
        FMNsq(t) = 0
        N2_C1(t)
        N2r_C1(t) = 0
        N3_C1(t)
        N3r_C1(t) = 0
        N1a_C1(t)
        N1ar_C1(t) = 0
        KEQ13_C1(t)
        vQ_C1(t)
        vQH2_C1(t)
        vNADH_C1(t)
        vNAD_C1(t)
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
        TN_C1(t) ## NADH turnover number
        vHres_C1(t)
    end

    fhm = h_m / 1E-7Molar
    C1_INHIB = (1 - ROTENONE_BLOCK)
    ## NADH + FMN = FMN.NADH
    v1 = kf1_C1 * (nadh * FMN - FMN_NADH / KEQ1_C1)
    ## FMN.NADH = FMNH−.NAD+
    v2 = kf2_C1 * (FMN_NADH - FMNH_NAD / KEQ2_C1)
    ## FMNH−.NAD+ = FMNH− + NAD+
    v3 = kf3_C1 * (FMNH_NAD - FMNH * nad / KEQ3_C1)
    ## FMN + NAD+ = FMN.NAD+
    v4 = kf4_C1 * (nad * FMN - FMN_NAD / KEQ4_C1)
    ## FMNH− + NADH = FMNH−.NADH
    v5 = kf5_C1 * (FMNH * nadh - FMNH_NADH / KEQ5_C1)
    ## FMNH− + N3 = FMNHsq + N3−
    v6 = kf6_C1 * (FMNH * N3_C1 - FMNsq * N3r_C1 / KEQ6_C1)
    ## N3− + N2 = N3 + N2−
    v7 = kf7_C1 * (N3r_C1 * N2_C1 - N3_C1 * N2r_C1 / KEQ7_C1)
    ## Q association
    q = Q_n * C1_INHIB
    v8 = kf8_C1 * (Iq_C1 * q - Q_C1 / KEQ8_C1)
    ## CI.Q + N2− = CIQsq + N2
    v9 = kf9_C1 * (Q_C1 * N2r_C1 - SQ_C1 * N2_C1 / KEQ9_C1)
    ## FMNHsq + N1a = FMN + N1a− + Hi+
    v10 = kf10_C1 * (FMNsq * N1a_C1 - FMN * N1ar_C1 * fhm / KEQ10_C1)
    ## FMNsq + N1a− = FMNH- + N1a
    v10b = kf10_C1 * (FMNsq * N1ar_C1 - FMNH * N1a_C1 / KEQ10B_C1)
    ## FMNHsq + N3 = FMN + N3− + Hi+
    v11 = kf11_C1 * (FMNsq * N3_C1 - FMN * N3r_C1 * fhm / KEQ11_C1)
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = kf13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ13_C1)
    ## QH2 dissociation
    qh2 = QH2_n * C1_INHIB
    v14 = kf14_C1 * (QH2_C1 - Iq_C1 * qh2 / KEQ14_C1)
    ## Flavin site ROS generation
    v16 = kf16_C1 * (FMNH * O2 - FMNsq * sox_m / KEQ16_C1)
    ## Quinone site ROS generation
    v17 = kf17_C1 * (SQ_C1 * O2 - Q_C1 * sox_m / KEQ17_C1)

    eqs = [
        ET_C1 ~ N2r_C1 + N2_C1,
        ET_C1 ~ N3r_C1 + N3_C1,
        ET_C1 ~ N1ar_C1 + N1a_C1,
        ET_C1 ~ Iq_C1 + Q_C1 + SQ_C1 + QH2_C1,
        ET_C1 ~ FMN + FMN_NADH + FMNH_NAD + FMN_NAD + FMNH_NADH + FMNH + FMNsq,
        KEQ13_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_m / h_i)^4,
        D(FMN_NADH) ~ v1 - v2,
        D(FMNH_NAD) ~ v2 - v3,
        D(FMN_NAD) ~ v4,
        D(FMNH_NADH) ~ v5,
        D(FMNH) ~ v3 - v5 - v16 - v6 + v10b,
        D(FMNsq) ~ v6 + v16 - v10 - v11 - v10b,
        D(N1ar_C1) ~ v10 - v10b,
        D(N3r_C1) ~ v6 + v11 - v7 - v12,
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        D(Q_C1) ~ v8 - v9 + v17,
        D(SQ_C1) ~ v9 - v17 - v13,
        D(QH2_C1) ~ v13 - v14,

        ## Positive: production; negative: consumption
        vNADH_C1 ~ -(v1 + v5),
        vNAD_C1 ~ v3 - v4,
        vQ_C1 ~ -v8,
        vQH2_C1 ~ v14,
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROS_C1 ~ vROSIf + vROSIq,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

# Q-site complex I model
function c1q(; name=:c1q,
    Q_n=1800μM, QH2_n=200μM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)

    @parameters begin
        ET_C1 = 1μM                ## Activity of complex I
        Em_O2_SOX = -160mV          ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV       ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV      ## FMN semiquinone/FMNH- redox potential
        ## FMN/FMNH- avg redox potential
        Em_FMN_FMNH = (Em_FMN_FMNsq + Em_FMNsq_FMNH)/2
        Em_NAD = -320mV             ## NAD/NADH avg redox potential
        Em_Q_SQ_C1 = -300mV         ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV       ## 800mV in Markevich, 2015
        Em_N3 = -250mV
        Em_N2 = -80mV
        KI_NADH_C1 = 50μM
        KD_NADH_C1 = 100μM
        KI_NAD_C1 = 1000μM
        KD_NAD_C1 = 25μM
        ## NADH + FMN = NAD+ + FMNH-
        KEQ_NADH_FMN = exp(2iVT * (Em_FMN_FMNH - Em_NAD))
        ## 2FMNsq = (N1a) = FMN + FMNH- + H+
        rKEQ_FMNsq_Dis = exp(-iVT * (Em_FMNsq_FMNH - Em_FMN_FMNsq))
        ## FMNH- + N2 = FMNsq + N2-
        KEQ_FMNH_N2 = exp(iVT * (Em_N2 - Em_FMNsq_FMNH))
        ## N2r + Q = N2 + SQ
        kf7_C1 = 10000Hz
        rKEQ_N2r_Q = exp(-iVT * (Em_Q_SQ_C1 - Em_N2))
        ## N2r + SQ = N2 + QH2
        kf13_C1 = 2.7e6Hz
        ## Q binding and QH2 unbinding
        kf14_C1 = 10Hz
        rKD_Q_C1 = inv(10μM)
        rKD_QH2_C1 = inv(20μM)
        ## SOX production from IF site
        kf16_C1 = 0.001Hz / μM
        rKEQ16_C1 = exp(-iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        ## SOX production from IQ site
        kf17_C1 = 0.001Hz / μM / 20
        rKEQ17_C1 = exp(-iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        ## Flavin site
        FMN(t)
        FMN_NAD(t)
        FMNsq(t)
        FMNH(t)
        FMNH_NADH(t)
        FMN_NADH(t)
        FMNH_NAD(t)
        N2_C1(t)
        N2r_C1(t)
        ## Quinone site
        Q_C1(t)
        SQ_C1(t) = 0
        QH2_C1(t) = 0
        rKEQ_N2r_SQ(t)
        ## Reaction rates
        vQ_C1(t)
        vROS_C1(t)
        vROSIf(t)
        vROSIq(t)
        vNADH_C1(t)
        TN_C1(t)
    end

    ## Mitochondrial pH factor
    fhm = h_m * inv(1E-7Molar)
    ## Weights in the flavin site
    wFMN = 1
    wFMN_NAD = wFMN * nad / KI_NAD_C1
    wFMN_NADH = wFMN * nadh / KD_NADH_C1
    wFMNH = wFMN * (nadh / nad) * KEQ_NADH_FMN
    wFMNH_NAD = wFMNH * nad / KD_NAD_C1
    wFMNH_NADH = wFMNH * nadh / KI_NADH_C1
    wFMNsq = NaNMath.sqrt(wFMN * wFMNH * rKEQ_FMNsq_Dis * fhm)
    denf = wFMN + wFMN_NAD + wFMNH + wFMNH_NADH + wFMNsq + wFMN_NADH + wFMNH_NAD
    ## First electron transfer
    v7 = kf7_C1 * (N2r_C1 * Q_C1 - N2_C1 * SQ_C1 * rKEQ_N2r_Q)
    ## Second electron transfer
    v13 = kf13_C1 * (N2r_C1 * SQ_C1 * fhm^2 - N2_C1 * QH2_C1 * rKEQ_N2r_SQ)
    ## Q binding and QH2 unbinding
    q = Q_n * rKD_Q_C1
    qh2 = QH2_n * rKD_QH2_C1
    v14 = kf14_C1 * (QH2_C1 * q - Q_C1 * qh2)
    ## Flavin site ROS generation
    v16 = kf16_C1 * (FMNH * O2 - FMNsq * sox_m * rKEQ16_C1)
    ## Quinone site ROS generation
    v17 = kf17_C1 * (SQ_C1 * O2 - Q_C1 * sox_m * rKEQ17_C1)

    eqs = [
        rKEQ_N2r_SQ ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
        FMN ~ wFMN * ET_C1 / denf,
        FMN_NAD ~ wFMN_NAD * ET_C1 / denf,
        FMNH ~ wFMNH * ET_C1 / denf,
        FMNsq ~ wFMNsq * ET_C1 / denf,
        FMNH_NADH ~ wFMNH_NADH * ET_C1 / denf,
        FMN_NADH ~ wFMN_NADH * ET_C1 / denf,
        FMNH_NAD ~ wFMNH_NAD * ET_C1 / denf,
        N2_C1 ~ FMNsq / (FMNsq + FMNH * KEQ_FMNH_N2),
        N2r_C1 ~ 1 - N2_C1,
        ET_C1 ~ Q_C1 + SQ_C1 + QH2_C1,
        D(SQ_C1) ~ v7 - v13 - v17,
        D(QH2_C1) ~ v13 - v14,
        vQ_C1 ~ -v14,
        vNADH_C1 ~ -0.5 * (v7 + v13 + v16),
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROS_C1 ~ vROSIf + vROSIq,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

#---
@parameters begin
    Q_n = 1.8mM
    QH2_n = 0.2mM
    nad = 500μM
    nadh = 500μM
    dpsi = 150mV
end

#---
qsys = c1q(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
markevich = c1_markevich_full(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
gauthier = c1_gauthier(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify

prob_q = SteadyStateProblem(qsys, [
    qsys.ET_C1 => 1μM,
    qsys.kf7_C1 => 10000Hz,
    qsys.kf13_C1 => 2.7e6Hz,
    qsys.kf14_C1 => 10Hz,
    qsys.kf16_C1 => 0.020Hz / μM,
    qsys.kf17_C1 => 0.017Hz / μM / 20,
])
prob_m = SteadyStateProblem(markevich, [
    markevich.ET_C1 => 17μM,
    markevich.kf16_C1 => 0.001Hz / μM,
    markevich.kf17_C1 => 0.001Hz / μM / 20,
])
prob_g = SteadyStateProblem(gauthier, [])
alg = DynamicSS(Rodas5P())
ealg = EnsembleThreads()

# ## Varying MMP
dpsirange = 100mV:5mV:200mV
alter_dpsi = (prob, i, repeat) -> begin
    prob.ps[dpsi] = dpsirange[i]
    prob
end

eprob_q = EnsembleProblem(prob_q; prob_func=alter_dpsi)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_dpsi)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi)
@time sim_q = solve(eprob_q, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# Helper function
extract(sim, k) = map(s -> s[k], sim)
# MMP vs NADH turnover
# markevich model has a steeper dependence
xs = dpsirange
ys = hcat(extract(sim_g, gauthier.vNADH_C1), extract(sim_m, markevich.vNADH_C1), extract(sim_q, qsys.vNADH_C1))

plot(xs, ys, xlabel="MMP (mV)", ylabel="NADH rate (μM/ms)", label=["Gauthier" "Markevich" "IQ"])

# Turnover rate (Hz)
extract(sim_q, qsys.TN_C1) .* 1000

# IF redox potential (mV)
extract(sim_q, -80 + 26.7 * log(qsys.N2_C1 / qsys.N2r_C1))
#---
ys = stack(extract.(Ref(sim_q), [qsys.Q_C1, qsys.SQ_C1, qsys.QH2_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1"], legend=:left)

#---
ys = stack(extract.(Ref(sim_q), [qsys.FMN, qsys.FMNsq, qsys.FMNH, qsys.FMN_NAD, qsys.FMNH_NADH, qsys.FMN_NADH, qsys.FMNH_NAD]), dims=2)
pl1 = plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH" "FMN_NADH" "FMNH_NAD"], title="IQ model", legend=:right)

ys = stack(extract.(Ref(sim_m), [markevich.FMN, markevich.FMNsq, markevich.FMNH, markevich.FMN_NAD, markevich.FMNH_NADH, markevich.FMN_NADH, markevich.FMNH_NAD]), dims=2)
pl2 = plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH" "FMN_NADH" "FMNH_NAD"], title="M model", legend=:right)

plot(pl1, pl2)

# MMP vs ROS production
xs = dpsirange
ys_g = extract(sim_g, gauthier.vROS_C1) .* 1000
ys_m = extract(sim_m, markevich.vROS_C1) .* 1000
ys_q = extract(sim_q, qsys.vROS_C1) .* 1000
plot(xs, [ys_g ys_m ys_q], xlabel="MMP (mV)", ylabel="ROS production (μM/s)", label=["Gauthier" "Markevich" "IQ"])

# ## Varying NADH
nadhrange = 10μM:10μM:990μM
alter_nadh = (prob, i, repeat) -> begin
    prob.ps[nadh] = nadhrange[i]
    prob.ps[nad] = 1000μM - prob.ps[nadh]
    prob
end

eprob_q = EnsembleProblem(prob_q; prob_func=alter_nadh)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_nadh)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_nadh)
@time sim_q = solve(eprob_q, alg, ealg; trajectories=length(nadhrange), abstol=1e-8, reltol=1e-8)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(nadhrange), abstol=1e-8, reltol=1e-8)
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(nadhrange), abstol=1e-8, reltol=1e-8)

# NADH vs turnover
xs = nadhrange
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_q = extract(sim_q, qsys.vNADH_C1)

plot(xs, [ys_g ys_m ys_q], xlabel="NADH (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich" "IQ"])

# NADH vs ROS production
xs = nadhrange
ys = [extract(sim_g, gauthier.vROS_C1) extract(sim_m, markevich.vROS_C1) extract(sim_q, qsys.vROS_C1)]

plot(xs, ys, xlabel="NADH (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "IQ"])

#---
ys = stack(extract.(Ref(sim_m), [markevich.Q_C1, markevich.SQ_C1, markevich.QH2_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1"], legend=:left, title="M model")

#---
ys = stack(extract.(Ref(sim_m), [markevich.FMN, markevich.FMNsq, markevich.FMNH, markevich.FMN_NAD, qsys.FMNH_NADH]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH"], legend=:topleft, title="M model")

#---
ys = stack(extract.(Ref(sim_q), [qsys.Q_C1, qsys.SQ_C1, qsys.QH2_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1"], title="Iq model")

#---
ys = stack(extract.(Ref(sim_q), [qsys.FMN, qsys.FMNsq, qsys.FMNH, qsys.FMN_NAD, qsys.FMNH_NADH]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH"])

# ## Varying Q
qh2range = 10μM:10μM:1990μM
alter_qh2 = (prob, i, repeat) -> begin
    prob.ps[QH2_n] = qh2range[i]
    prob.ps[Q_n] = 2000μM - prob.ps[QH2_n]
    prob
end

eprob_q = EnsembleProblem(prob_q; prob_func=alter_qh2)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_qh2)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_qh2)
@time sim_q = solve(eprob_q, alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# QH2 vs NADH turnover
xs = qh2range
ys = [extract(sim_g, gauthier.vNADH_C1) extract(sim_m, markevich.vNADH_C1) extract(sim_q, qsys.vNADH_C1)]

plot(xs, ys, xlabel="QH2 (μM)", ylabel="NADH rate (μM/ms)", label=["Gauthier" "Markevich" "IQ"])

# QH2 vs Q turnover
xs = qh2range
ys = [extract(sim_g, gauthier.vQ_C1) extract(sim_m, markevich.vQ_C1) extract(sim_q, qsys.vQ_C1)]

plot(xs, ys, xlabel="QH2 (μM)", ylabel="Q rate (μM/ms)", label=["Gauthier" "Markevich" "IQ"])

# QH2 vs ROS production
# Gauthier model produces a lot of SOX on high QH2
xs = qh2range
ys = [extract(sim_g, gauthier.vROS_C1) extract(sim_m, markevich.vROS_C1) extract(sim_q, qsys.vROS_C1)]
plot(xs, ys, xlabel="QH2 (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "IQ"])

#---
ys = stack(extract.(Ref(sim_q), [qsys.Q_C1, qsys.SQ_C1, qsys.QH2_C1]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1"], legend=:left)

#---
@unpack C1_1, C1_2, C1_3, C1_4, C1_5, C1_6, C1_7 = gauthier
ys = stack(extract.(Ref(sim_g), [C1_3, C1_4, C1_6]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc (μM)", label=["C1_3" "C1_4" "C1_6"], legend=:right, lw=1.5)
