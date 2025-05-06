using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using NaNMath
using Plots
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms

# Simplified Markevich model
# QSSA for the catalytic cycles in the flavin and the quinone sites
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
