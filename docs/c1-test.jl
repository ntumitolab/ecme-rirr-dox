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

# From Markevich, 2015
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
        KEQ2_C1 = exp(iVT * (Em_FMN_FMNH - Em_NAD)) / KEQ1_C1 / KEQ3_C1
        kf4_C1 = 1Hz / μM
        KEQ4_C1 = 0.001 / μM
        kf5_C1 = 2Hz / μM
        KEQ5_C1 = 0.02 / μM
        kf6_C1 = 5e8Hz / μM
        KEQ6_C1 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
        kf7_C1 = 1E4Hz / μM
        KEQ7_C1 = exp(iVT * (Em_N2 - Em_N3))
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = 0.1 / μM         ## Association constant for Q
        kf9_C1 = 4E5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kf10_C1 = 2e6Hz / μM
        KEQ10_C1 = exp(iVT * (Em_N1a - Em_FMN_FMNsq))
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
    v8 = kf8_C1 * (Iq_C1 * Q_n - Q_C1 / KEQ8_C1)
    ## CI.Q + N2− = CIQsq + N2
    v9 = kf9_C1 * (Q_C1 * N2r_C1 - SQ_C1 * N2_C1 / KEQ9_C1)
    ## FMNHsq + N1a = FMN + N1a− + Hi+
    v10 = kf10_C1 * (FMNsq * N1a_C1 - FMN * N1ar_C1 * fhm / KEQ10_C1)
    ## FMNHsq + N3 = FMN + N3− + Hi+
    v11 = kf11_C1 * (FMNsq * N3_C1 - FMN * N3r_C1 * fhm / KEQ11_C1)
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = kf13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ13_C1)
    ## QH2 dissociation
    v14 = kf14_C1 * (QH2_C1 - Iq_C1 * QH2_n / KEQ14_C1)
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
        D(FMNH) ~ v3 - v5 - v16 - v6,
        D(FMNsq) ~ v6 + v16 - v10 - v11,
        D(N1ar_C1) ~ v10,
        D(N3r_C1) ~ v6 + v11 - v7 - v12,
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        D(Q_C1) ~ v8 - v9 + v17,
        D(SQ_C1) ~ v9 - v17 - v13,
        D(QH2_C1) ~ v13 - v14,

        ## Positive: production; negative: consumption
        vNADH_C1 ~ -v1,
        vNAD_C1 ~ v3,
        vQ_C1 ~ -v8,
        vQH2_C1 ~ v14,
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROS_C1 ~ vROSIf + vROSIq,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

# Simplified Markevich model
# QSSA for the catalytic cycle in the flavin (FMN) site
function c1_birb(; name=:c1birb,
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
        Em_SQ_QH2_C1 = +500mV     ## 800mV (?) in Markevich, 2015
        ET_C1 = 17μM              ## Activity of complex I
        ## DOX IC50 on complex I
        KI_DOX_C1 = 400μM
        kf1_C1 = 83Hz / μM
        KEQ1_C1 = 0.01 / μM
        kr1_C1 = kf1_C1 / KEQ1_C1
        ## kf2_C1 = 1.44e12Hz # Rapid equlibrium
        kf3_C1 = 1e6Hz
        KEQ3_C1 = 25μM
        kr3_C1 = kf3_C1 / KEQ3_C1
        KEQ2_C1 = exp(iVT * (Em_FMN_FMNH - Em_NAD)) / KEQ1_C1 / KEQ3_C1
        kf6_C1 = 5e8Hz / μM
        KEQ6_C1 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
        kr6_C1 = kf6_C1 / KEQ6_C1
        kf7_C1 = 1E4Hz / μM
        KEQ7_C1 = exp(iVT * (Em_N2 - Em_N3))
        kr7_C1 = kf7_C1 / KEQ7_C1
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = 0.1 / μM         ## Association constant for Q
        kr8_C1 = kf8_C1 / KEQ8_C1
        kf9_C1 = 4E5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kr9_C1 = kf9_C1 / KEQ9_C1
        kf10_C1 = 2e6Hz / μM
        KEQ10_C1 = exp(iVT * (Em_N1a - Em_FMN_FMNsq))
        kr10_C1 = kf10_C1 / KEQ10_C1
        kf11_C1 = 1e9Hz / μM
        KEQ11_C1 = exp(iVT * (Em_N3 - Em_FMN_FMNsq))
        kr11_C1 = kf11_C1 / KEQ11_C1
        kf13_C1 = 2.7e6Hz / μM
        kf14_C1 = 1000Hz
        KEQ14_C1 = 20μM             ## Dissociation constant for QH2
        kr14_C1 = kf14_C1 / KEQ14_C1
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kr16_C1 = kf16_C1 / KEQ16_C1
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        KEQ17_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
        kr17_C1 = kf17_C1 / KEQ17_C1
        KI_NAD_C1 = 1mM
        KI_NADH_C1 = 50μM
    end

    @variables begin
        N2_C1(t)
        N2r_C1(t) = 0
        N3_C1(t)
        N3r_C1(t) = 0
        N1a_C1(t)
        N1ar_C1(t) = 0
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
        vHres_C1(t)
    end

    @variables begin
        wFMN_FMNNAD(t)
        wFMNNADH_FMNHNAD(t)
        wFMNH_FMNHNADH(t)
        wFMNsq(t)
        FMN(t)
        FMN_NADH(t)
        FMNH_NAD(t)
        FMN_NAD(t)
        FMNH_NADH(t)
        FMNH(t)
        FMNsq(t)
        TN_C1f(t) ## NADH turnover number
        TN_C1(t)  ## Flavin site turnover number
        vNADH_C1(t)
        vNAD_C1(t)
    end

    fhm = h_m * inv(1E-7Molar)
    fFMN = KI_NAD_C1 / (nad + KI_NAD_C1)
    fFMNH = KI_NADH_C1 / (nadh + KI_NADH_C1)
    fFMN_NADH = 1 / (1 + KEQ2_C1)
    fFMNH_NAD = KEQ2_C1 / (1 + KEQ2_C1)
    fDen = wFMN_FMNNAD + wFMNNADH_FMNHNAD + wFMNH_FMNHNADH + wFMNsq
    fC1 = ET_C1 / fDen

    ## State transition rates in the flavin site
    ## 1 = FMN + FMN_NAD, 2 = FMN_NADH + FMNH_NAD, 3 = FMNH + FMNH_NADH, 4 = FMNsq
    a12 = kf1_C1 * nadh * fFMN
    a21 = kr1_C1 * fFMN_NADH
    a23 = kf3_C1 * fFMNH_NAD
    a32 = kr3_C1 * nad * fFMNH
    a34 = (kf6_C1 * N3_C1 + kf16_C1 * O2) * fFMNH
    a43 = kr6_C1 * N3r_C1 + kr16_C1 * sox_m
    a41 = kf11_C1 * N3_C1 + kf10_C1 * N1a_C1
    a14 = (kr11_C1 * N3r_C1 + kr10_C1 * N1ar_C1) * fhm * fFMN

    eqsf = [
        wFMN_FMNNAD ~ a21*a32*a41 + a21*a32*a43 + a21*a34*a41 + a23*a34*a41,
        wFMNNADH_FMNHNAD ~ a12*a32*a41 + a12*a32*a43 + a12*a34*a41 + a14*a32*a43,
        wFMNH_FMNHNADH ~ a12*a23*a41 + a12*a23*a43 + a14*a21*a43 + a14*a23*a43,
        wFMNsq ~ a12*a23*a34 + a14*a21*a32 + a14*a21*a34 + a14*a23*a34,
        FMN ~ fC1 * wFMN_FMNNAD * fFMN,
        FMN_NAD ~ fC1 * wFMN_FMNNAD * (1 - fFMN),
        FMN_NADH ~ fC1 * wFMNNADH_FMNHNAD * fFMN_NADH,
        FMNH_NAD ~ fC1 * wFMNNADH_FMNHNAD * fFMNH_NAD,
        FMNH_NADH ~ fC1 * wFMNH_FMNHNADH * (1 - fFMNH),
        FMNH ~ fC1 * wFMNH_FMNHNADH * fFMNH,
        FMNsq ~ fC1 * wFMNsq,
        TN_C1 ~ (a12 * a23 * a34 * a41 - a14 * a43 * a32 * a21) / fDen,
    ]

    ## State transition rates in the quinone site
    ## 1 = Iq 2 = IqQ, 3 = IqSQ, 4 = IqQH2
    @variables begin
        Iq_C1(t)
        Q_C1(t)
        SQ_C1(t)
        QH2_C1(t)
        wIq(t)
        wIqQ(t)
        wIqSQ(t)
        wIqQH2(t)
        vQ_C1(t)
        vQH2_C1(t)
        KEQ13_C1(t)
        vHres_C1(t)
    end

    b12 = kf8_C1 * Q_n
    b21 = kr8_C1
    b23 = kf9_C1 * N2r_C1 + kr17_C1 * sox_m
    b32 = kr9_C1 * N2_C1 + kf17_C1 * O2
    b34 = kf13_C1 * N2r_C1 * fhm^2
    b43 = kf13_C1 / KEQ13_C1 * N2_C1
    b41 = kf14_C1
    b14 = kr14_C1 * QH2_n
    qDen = wIq + wIqQ + wIqSQ + wIqQH2
    qC1 = ET_C1 / qDen

    eqsq = [
        KEQ13_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_m / h_i)^4,
        wIq ~ b21*b32*b41 + b21*b32*b43 + b21*b34*b41 + b23*b34*b41,
        wIqQ ~ b12*b32*b41 + b12*b32*b43 + b12*b34*b41 + b14*b32*b43,
        wIqSQ ~ b12*b23*b41 + b12*b23*b43 + b14*b21*b43 + b14*b23*b43,
        wIqQH2 ~ b12*b23*b34 + b14*b21*b32 + b14*b21*b34 + b14*b23*b34,
        Iq_C1 ~ wIq * qC1,
        Q_C1 ~ wIqQ * qC1,
        SQ_C1 ~ wIqSQ * qC1,
        QH2_C1 ~ wIqQH2 * qC1,
        vQ_C1 ~ kr8_C1 * Q_C1 - kf8_C1 * Q_n * Iq_C1,
        vQH2_C1 ~ kf14_C1 * QH2_C1 - kr14_C1 * QH2_n * Iq_C1,
    ]

    ## NADH + FMN = FMN.NADH
    v1 = kf1_C1 * nadh * FMN - kr1_C1 * FMN_NADH
    ## FMNH−.NAD+ = FMNH− + NAD+
    v3 = kf3_C1 * FMNH_NAD - kr3_C1 * FMNH * nad
    ## FMNH− + N3 = FMNHsq + N3−
    v6 = kf6_C1 * FMNH * N3_C1 - kr6_C1 * FMNsq * N3r_C1
    ## N3− + N2 = N3 + N2−
    v7 = kf7_C1 * N3r_C1 * N2_C1 - kr7_C1 * N3_C1 * N2r_C1
    ## First electron transfer CI.Q + N2− = CIQsq + N2
    v9 = kf9_C1 * Q_C1 * N2r_C1 - kr9_C1 * SQ_C1 * N2_C1
    ## FMNHsq + N1a = FMN + N1a− + Hi+
    v10 = kf10_C1 * FMNsq * N1a_C1 - kr10_C1 * FMN * N1ar_C1 * fhm
    ## FMNHsq + N3 = FMN + N3− + Hi+
    v11 = kf11_C1 * FMNsq * N3_C1 - kr11_C1 * FMN * N3r_C1 * fhm
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = kf13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ13_C1)
    ## Flavin site ROS generation
    v16 = kf16_C1 * FMNH * O2 - kr16_C1 * FMNsq * sox_m
    ## Quinone site ROS generation
    v17 = kf17_C1 * SQ_C1 * O2 - kr17_C1 * Q_C1 * sox_m

    eqs = [
        ET_C1 ~ N2r_C1 + N2_C1,
        ET_C1 ~ N3r_C1 + N3_C1,
        ET_C1 ~ N1ar_C1 + N1a_C1,
        D(N1ar_C1) ~ v10,
        D(N3r_C1) ~ v6 + v11 - v7 - v12,
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROS_C1 ~ vROSIf + vROSIq,
        vHres_C1 ~ 4v13,
        vNADH_C1 ~ -v1,
        vNAD_C1 ~ v3,
    ]
    return ODESystem([eqsf; eqsq; eqs], t; name)
end

# From Gauthier 2012
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
        vNADH_C1 ~ -v23,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

@parameters begin
    Q_n = 1.8mM
    QH2_n = 0.2mM
    nad = 500μM
    nadh = 500μM
    dpsi = 150mV
end

#---
birb = c1_birb(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
markevich = c1_markevich_full(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
gauthier = c1_gauthier(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify

prob_m = SteadyStateProblem(markevich, [markevich.ET_C1 => 17μM, markevich.kf16_C1 => 0.001Hz / μM, markevich.kf17_C1 => 0.001Hz / μM / 20])
prob_b = SteadyStateProblem(birb, [birb.ET_C1 => 17μM, birb.kf16_C1 => 0.001Hz / μM, birb.kf17_C1 => 0.001Hz / μM / 20])
prob_g = SteadyStateProblem(gauthier, [])
alg = DynamicSS(TRBDF2())
ealg = EnsembleSerial()

# ## Varying MMP
dpsirange = 100mV:5mV:200mV
alter_dpsi = (prob, i, repeat) -> remake(prob, p=[dpsi => dpsirange[i]])

eprob_b = EnsembleProblem(prob_b; prob_func=alter_dpsi, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_dpsi, safetycopy=false)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi, safetycopy=false)
@time sim_b = solve(eprob_b, alg, ealg; trajectories=length(dpsirange))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(dpsirange))
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(dpsirange))

# Helper function
extract(sim, k) = map(s -> s[k], sim)
# MMP vs NADH turnover
# markevich model has a steeper dependence
xs = dpsirange
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_b = extract(sim_b, birb.vNADH_C1)

plot(xs, [ys_g ys_m ys_b], xlabel="MMP (mV)", ylabel="NADH rate (μM/ms)", label=["Gauthier" "Markevich" "Birb"])

# MMP vs ROS production
xs = dpsirange
ys_g = extract(sim_g, gauthier.vROS_C1) .* 1000
ys_m = extract(sim_m, markevich.vROS_C1) .* 1000
ys_b = extract(sim_b, birb.vROS_C1) .* 1000
plot(xs, [ys_g ys_m ys_b], xlabel="MMP (mV)", ylabel="ROS production (μM/s)", label=["Gauthier" "Markevich" "Birb"])

#---
ys_if = extract(sim_b, birb.vROSIf)
ys_iq = extract(sim_b, birb.vROSIq)
plot(xs, [ys_if ys_iq], title="Birb model", xlabel="MMP (mV)", ylabel="ROS production", label=["IF" "IQ"])

# Inside the ODE system
ys = stack(extract.(Ref(sim_b), [birb.Q_C1, birb.SQ_C1, birb.QH2_C1, birb.Iq_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "Iq_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_b), [birb.N1ar_C1, birb.N3r_C1, birb.N2r_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
_em(em7, xr, xo) = em7 - 26.7 * log(xr / xo)
ys = stack(extract.(Ref(sim_b), [ _em(-370, birb.N1ar_C1, birb.N1a_C1), _em(-250, birb.N3r_C1, birb.N3_C1), _em(-80, birb.N2r_C1, birb.N2_C1)]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Redox potential (mV)", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
@unpack FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq = birb
ys = stack(extract.(Ref(sim_b), [FMNH, FMNsq]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Conc (μM)", label=["FMNH" "FMNsq"], legend=:right, lw=1.5)

# ## Varying NADH
nadhrange = 10μM:10μM:990μM
nadrange = 1000μM .- nadhrange
alter_nadh = (prob, i, repeat) -> remake(prob, p=[nadh => nadhrange[i], nad => nadrange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_nadh, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_nadh, safetycopy=false)
eprob_b = EnsembleProblem(prob_b; prob_func=alter_nadh, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(nadhrange))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(nadhrange))
@time sim_b = solve(eprob_b, alg, ealg; trajectories=length(nadhrange))

# NADH vs turnover
xs = nadhrange
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_b = extract(sim_b, birb.vNADH_C1)

plot(xs, [ys_g ys_m ys_b], xlabel="NADH (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich" "Birb"])

# NADH vs ROS production
xs = nadhrange
ys_g = extract(sim_g, gauthier.vROS_C1)
ys_m = extract(sim_m, markevich.vROS_C1)
ys_b = extract(sim_b, birb.vROS_C1)

plot(xs, [ys_g ys_m ys_b], xlabel="NADH (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "Birb"])

#---
ys_if = extract(sim_b, birb.vROSIf)
ys_iq = extract(sim_b, birb.vROSIq)
plot(xs, [ys_if ys_iq], title="Birb model", xlabel="NADH (μM)", ylabel="ROS production", label=["IF" "IQ"])

# Inside the ODE system
ys = stack(extract.(Ref(sim_b), [birb.Q_C1, birb.SQ_C1, birb.QH2_C1, birb.Iq_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Fraction", label=["Q_C1" "SQ_C1" "QH2_C1" "Iq_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_b), [birb.N1ar_C1, birb.N3r_C1, birb.N2r_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Fraction", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
ys = extract(sim_b, birb.FMNsq / birb.FMN)
plot(xs, ys, xlabel="NADH (μM)", ylabel="FMNsq / FMN", label=false)

#---
@unpack FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq = birb
ys = stack(extract.(Ref(sim_b), [FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Conc (μM)", label=["FMN" "FMN_NADH" "FMNH_NAD" "FMN_NAD" "FMNH_NADH" "FMNH" "FMNsq"], legend=:right, lw=1.5)

# ## Varying Q
qh2range = 10μM:10μM:1990μM
qrange = 2000μM .- qh2range
alter_qh2 = (prob, i, repeat) -> remake(prob, p=[QH2_n => qh2range[i], Q_n => qrange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_qh2, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_qh2, safetycopy=false)
eprob_b = EnsembleProblem(prob_b; prob_func=alter_qh2, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(qh2range))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(qh2range))
@time sim_b = solve(eprob_b, alg, ealg; trajectories=length(qh2range))

# QH2 vs NADH turnover
xs = qh2range
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_b = extract(sim_b, birb.vNADH_C1)

plot(xs, [ys_g ys_m ys_b], xlabel="QH2 (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich" "Birb"])

# QH2 vs ROS production
# Gauthier model produces a lot of SOX on high QH2
xs = qh2range
ys_g = extract(sim_g, gauthier.vROS_C1)
ys_m = extract(sim_m, markevich.vROS_C1)
ys_b = extract(sim_b, birb.vROS_C1)
plot(xs, [ys_g ys_m ys_b], xlabel="QH2 (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "Birb"])

#---
@unpack C1_1, C1_2, C1_3, C1_4, C1_5, C1_6, C1_7 = gauthier
ys = stack(extract.(Ref(sim_g), [C1_3, C1_4, C1_6]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc (μM)", label=["C1_3" "C1_4" "C1_6"], legend=:right, lw=1.5)

#---
ys_if = extract(sim_b, birb.vROSIf)
ys_iq = extract(sim_b, birb.vROSIq)
plot(xs, [ys_if ys_iq], title="Birb model", xlabel="QH2 (μM)", ylabel="ROS production", label=["IF" "IQ"])

# Inside the ODE system
ys = stack(extract.(Ref(sim_b), [birb.Q_C1, birb.SQ_C1, birb.QH2_C1, birb.Iq_C1]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc. (μM)", label=["Q_C1" "SQ_C1" "QH2_C1" "Iq_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_b), [birb.N1ar_C1, birb.N3r_C1, birb.N2r_C1]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc. (μM)", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
@unpack FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq = birb
ys = stack(extract.(Ref(sim_b), [FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc (μM)", label=["FMN" "FMN_NADH" "FMNH_NAD" "FMN_NAD" "FMNH_NADH" "FMNH" "FMNsq"], legend=:right, lw=1.5)
