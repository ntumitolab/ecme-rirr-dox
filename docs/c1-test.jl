# # Complex I model
# Comparing Gauthier, Markevich, and simplified complex I models
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using NaNMath
using Plots
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms

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
    q = Q_n * (1 - ROTENONE_BLOCK)
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
    qh2 = QH2_n * (1 - ROTENONE_BLOCK)
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
    Q_n=1.8mM, QH2_n=0.2mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)

    @parameters begin
        ET_C1 = 17μM                ## Activity of complex I
        Em_O2_SOX = -160mV          ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV       ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV      ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV        ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV             ## NAD/NADH avg redox potential
        Em_N3 = -250mV              ## FeS N3 redox potential
        Em_N2 = -80mV               ## FeS N2 redox potential
        Em_Q_SQ_C1 = -300mV         ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV       ## 800mV in Markevich, 2015
        KI_NADH_C1 = 50μM
        KI_NAD_C1 = 1000μM
        ## NADH + FMN = NAD+ + FMNH-
        KEQ2_C1 = exp(2iVT * (Em_FMN_FMNH - Em_NAD))
        ## 0.5NADH + N3 = 0.5NAD + N3r + 0.5H+
        KEQ_NADH_N3 = exp(iVT * (Em_N3 - Em_NAD))
        rKEQ6_C1 = exp(-iVT * (Em_N3 - Em_FMNsq_FMNH))
        kf7_C1 = 10000Hz / μM
        rKEQ7_C1 = exp(-iVT * (Em_N2 - Em_N3))
        kf8_C1 = 10Hz / μM
        rKEQ8_C1 = 10μM         ## Dissociation constant for Q
        kf9_C1 = 4E5Hz / μM
        rKEQ9_C1 = exp(-iVT * (Em_Q_SQ_C1 - Em_N2))
        rKEQ11_C1 = exp(-iVT * (Em_N3 - Em_FMN_FMNsq))
        kf13_C1 = 2.7e6Hz / μM
        kf14_C1 = 1000Hz
        rKEQ14_C1 = inv(20μM)       ## Dissociation constant for QH2
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        rKEQ16_C1 = exp(-iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        rKEQ17_C1 = exp(-iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        ## Flavin site
        FMN(t)
        FMN_NAD(t)
        FMNsq(t)
        FMNH(t)
        FMNH_NADH(t)
        ## Quinone site
        Iq_C1(t)
        Q_C1(t)
        SQ_C1(t)
        QH2_C1(t)
        wIq(t)
        wIqQ(t)
        wIqSQ(t)
        wIqQH2(t)
        rKEQ13_C1(t)
        ## FeS clusters
        N2_C1(t) ## Conserved
        N2r_C1(t) = 0
        N3_C1(t)
        N3r_C1(t)
        ## N3r/N3 ratio from reaction 0.5NADH + N3 = 0.5NAD + N3r + 0.5H+
        rN3_C1(t)
        ## Reaction rates
        vQ_C1(t)
        vROS_C1(t)
        vROSIf(t)
        vROSIq(t)
        vNADH_C1(t)
        TN_C1(t)
    end

    ## Mitochondrial pH
    fhm = h_m * inv(1E-7Molar)
    ## N3− + N2 = N3 + N2−
    v7 = kf7_C1 * (N3r_C1 * N2_C1 - N3_C1 * N2r_C1 * rKEQ7_C1)
    ## Q association
    q = Q_n * (1 - ROTENONE_BLOCK)
    v8 = kf8_C1 * (Iq_C1 * q - Q_C1 * rKEQ8_C1)
    ## CI.Q + N2− = CIQsq + N2
    v9 = kf9_C1 * (Q_C1 * N2r_C1 - SQ_C1 * N2_C1 * rKEQ9_C1)
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = kf13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 * rKEQ13_C1)
    ## QH2 dissociation
    qh2 = QH2_n * (1 - ROTENONE_BLOCK)
    v14 = kf14_C1 * (QH2_C1 - Iq_C1 * qh2 * rKEQ14_C1)
    ## Flavin site ROS generation
    v16 = kf16_C1 * (FMNH * O2 - FMNsq * sox_m * rKEQ16_C1)
    ## Quinone site ROS generation
    v17 = kf17_C1 * (SQ_C1 * O2 - Q_C1 * sox_m * rKEQ17_C1)
    ## FMN + NADH = FMNH- + NAD+
    rFMNH_FMN = (nadh / nad) * KEQ2_C1
    ## FMNHsq + N3 = FMN + N3− + Hi+
    rFMNHsq_FMN = rN3_C1 * fhm * rKEQ11_C1
    ## Weights in the flavin site
    denf = 1 + rFMNH_FMN + rFMNHsq_FMN
    fFMN = KI_NAD_C1 / (nad + KI_NAD_C1)
    fFMNH = KI_NADH_C1 / (nadh + KI_NADH_C1)

    ## State transition rates in the quinone site
    ## 1 = Iq 2 = IqQ, 3 = IqSQ, 4 = IqQH2
    b12 = kf8_C1 * Q_n
    b21 = kf8_C1 * rKEQ8_C1
    b23 = kf9_C1 * N2r_C1 + kf17_C1 * rKEQ17_C1 * sox_m
    b32 = kf9_C1 * rKEQ9_C1 * N2_C1 + kf17_C1 * O2
    b34 = kf13_C1 * N2r_C1 * fhm^2
    b43 = kf13_C1 * rKEQ13_C1 * N2_C1
    b41 = kf14_C1
    b14 = kf14_C1 * rKEQ14_C1 * QH2_n
    qDen = wIq + wIqQ + wIqSQ + wIqQH2
    qC1 = ET_C1 / qDen

    eqs = [
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        rKEQ13_C1 ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
        wIq ~ b21*b32*b41 + b21*b32*b43 + b21*b34*b41 + b23*b34*b41,
        wIqQ ~ b12*b32*b41 + b12*b32*b43 + b12*b34*b41 + b14*b32*b43,
        wIqSQ ~ b12*b23*b41 + b12*b23*b43 + b14*b21*b43 + b14*b23*b43,
        wIqQH2 ~ b12*b23*b34 + b14*b21*b32 + b14*b21*b34 + b14*b23*b34,
        Iq_C1 ~ wIq * qC1,
        Q_C1 ~ wIqQ * qC1,
        SQ_C1 ~ wIqSQ * qC1,
        QH2_C1 ~ wIqQH2 * qC1,
        ET_C1 ~ N2_C1 + N2r_C1,
        ET_C1 ~ N3_C1 + N3r_C1,
        N3_C1 ~ ET_C1 / (1 + rN3_C1),
        rN3_C1 ~ KEQ_NADH_N3 * NaNMath.sqrt(nadh / (nad * fhm)),
        FMN ~ ET_C1 / denf * fFMN,
        FMN_NAD ~ ET_C1 / denf * (1 - fFMN),
        FMNsq ~ ET_C1 / denf * rFMNHsq_FMN,
        FMNH ~ ET_C1 / denf * rFMNH_FMN * fFMNH,
        FMNH_NADH ~ ET_C1 / denf * rFMNH_FMN * (1 - fFMNH),
        vQ_C1 ~ -v8,
        vNADH_C1 ~ -0.5 * (v7 + v12 + v16),
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

prob_q = SteadyStateProblem(qsys, [qsys.ET_C1 => 17μM, qsys.kf16_C1 => 0.001Hz / μM, qsys.kf17_C1 => 0.001Hz / μM / 20])
prob_m = SteadyStateProblem(markevich, [markevich.ET_C1 => 17μM, markevich.kf16_C1 => 0.001Hz / μM, markevich.kf17_C1 => 0.001Hz / μM / 20])
prob_g = SteadyStateProblem(gauthier, [])
alg = DynamicSS(Rodas5P())
ealg = EnsembleSerial()

# ## Varying MMP
dpsirange = 100mV:5mV:200mV
alter_dpsi = (prob, i, repeat) -> remake(prob, p=[dpsi => dpsirange[i]])

eprob_q = EnsembleProblem(prob_q; prob_func=alter_dpsi, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_dpsi, safetycopy=false)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi, safetycopy=false)
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

# MMP vs Q turnover
xs = dpsirange
ys = hcat(extract(sim_g, gauthier.vQ_C1), extract(sim_m, markevich.vQ_C1), extract(sim_q, qsys.vQ_C1))

plot(xs, ys, xlabel="MMP (mV)", ylabel="Q rate (μM/ms)", label=["Gauthier" "Markevich" "IQ"])

#---
ys = stack(extract.(Ref(sim_q), [qsys.Q_C1, qsys.SQ_C1, qsys.QH2_C1, qsys.N2r_C1, qsys.N3r_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "N2r_C1" "N3r_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_q), [qsys.FMN, qsys.FMNsq, qsys.FMNH, qsys.FMN_NAD, qsys.FMNH_NADH]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH"], legend=:right)

# MMP vs ROS production
xs = dpsirange
ys_g = extract(sim_g, gauthier.vROS_C1) .* 1000
ys_m = extract(sim_m, markevich.vROS_C1) .* 1000
ys_q = extract(sim_q, qsys.vROS_C1) .* 1000
plot(xs, [ys_g ys_m ys_q], xlabel="MMP (mV)", ylabel="ROS production (μM/s)", label=["Gauthier" "Markevich" "IQ"])

# ## Varying NADH
nadhrange = 10μM:10μM:990μM
nadrange = 1000μM .- nadhrange
alter_nadh = (prob, i, repeat) -> remake(prob, p=[nadh => nadhrange[i], nad => nadrange[i]])

eprob_q = EnsembleProblem(prob_q; prob_func=alter_nadh, safetycopy=false)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_nadh, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_nadh, safetycopy=false)
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
# Q site model is more sensitive
xs = nadhrange
ys = [extract(sim_g, gauthier.vROS_C1) extract(sim_m, markevich.vROS_C1) extract(sim_q, qsys.vROS_C1)]

plot(xs, ys, xlabel="NADH (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "IQ"])

#---
ys = stack(extract.(Ref(sim_q), [qsys.Q_C1, qsys.SQ_C1, qsys.QH2_C1, qsys.N2r_C1, qsys.N3r_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "N2r_C1" "N3r_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_q), [qsys.FMN, qsys.FMNsq, qsys.FMNH, qsys.FMN_NAD, qsys.FMNH_NADH]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH"], legend=:right)

# ## Varying Q
qh2range = 10μM:10μM:1990μM
qrange = 2000μM .- qh2range
alter_qh2 = (prob, i, repeat) -> remake(prob, p=[QH2_n => qh2range[i], Q_n => qrange[i]])

eprob_q = EnsembleProblem(prob_q; prob_func=alter_qh2, safetycopy=false)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_qh2, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_qh2, safetycopy=false)
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
@unpack C1_1, C1_2, C1_3, C1_4, C1_5, C1_6, C1_7 = gauthier
ys = stack(extract.(Ref(sim_g), [C1_3, C1_4, C1_6]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc (μM)", label=["C1_3" "C1_4" "C1_6"], legend=:right, lw=1.5)
