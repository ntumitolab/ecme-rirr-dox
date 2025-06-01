using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using NaNMath

# Markevich 2015 mass action model
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4426091/
function c1_markevich_full(; name=:c1markevich_full,
    Q_n=1.8mM, QH2_n=0.2mM,
    nad_m=500μM, nadh_m=500μM,
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
    v1 = kf1_C1 * (nadh_m * FMN - FMN_NADH / KEQ1_C1)
    ## FMN.NADH = FMNH−.NAD+
    v2 = kf2_C1 * (FMN_NADH - FMNH_NAD / KEQ2_C1)
    ## FMNH−.NAD+ = FMNH− + NAD+
    v3 = kf3_C1 * (FMNH_NAD - FMNH * nad_m / KEQ3_C1)
    ## FMN + NAD+ = FMN.NAD+
    v4 = kf4_C1 * (nad_m * FMN - FMN_NAD / KEQ4_C1)
    ## FMNH− + NADH = FMNH−.NADH
    v5 = kf5_C1 * (FMNH * nadh_m - FMNH_NADH / KEQ5_C1)
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
    return System(eqs, t; name)
end

# Simplified Markevich 2015 model
# Rapid equlibrium at the flavin site and QSSA in the quinone site
function c1_markevich_s(; name=:c1markevich_s,
    Q_n=1.8mM, QH2_n=0.2mM,
    nad_m=500μM, nadh_m=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0, MT_PROT=1)

    @parameters begin
        ET_C1 = 17μM              ## Activity of complex I
        KI_DOX_C1 = 400μM         ## DOX IC50 on complex I
        K_RC_DOX = 1000 / 15mM    ## DOX redox cycling constant
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## 800mV (?) in Markevich, 2015
        KI_NADH_C1 = 50μM
        KD_NADH_C1 = 100μM
        KI_NAD_C1 = 1000μM
        KD_NAD_C1 = 25μM
        KEQ_NADH_FMN = exp(2iVT * (Em_FMN_FMNH - Em_NAD))
        KEQ6_C1 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
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
        ## Reaction rates
        TNC1(t)  ## Flavin site turnover number
        vQC1(t)
        vQH2C1(t)
        vNADHC1(t)
        vNADC1(t)
        vROSIf(t)
        vROSIq(t)
        vROSC1(t)
        vHresC1(t)
    end

    C1_CONC = ET_C1 * MT_PROT
    ## complex I inhibition by DOX and rotenone
    C1_INHIB = hil(KI_DOX_C1, DOX, 3) * (1 - ROTENONE_BLOCK)
    ## Electron leak scaling factor from complex I
    E_LEAK_C1 = 1 + K_RC_DOX * DOX
    fhm = h_m * inv(1E-7Molar)
    ## Weights in the flavin site
    wFMN = 1
    wFMN_NAD = wFMN * nad_m / KI_NAD_C1
    wFMN_NADH = wFMN * nadh_m / KD_NADH_C1
    wFMNH = wFMN * (nadh_m / nad_m) * KEQ_NADH_FMN
    wFMNH_NAD = wFMNH * nad_m / KD_NAD_C1
    wFMNH_NADH = wFMNH * nadh_m / KI_NADH_C1
    wFMNsq = NaNMath.sqrt(wFMN * wFMNH * rKEQ_FMNsq_Dis * fhm)
    denf = wFMN + wFMN_NAD + wFMNH + wFMNH_NADH + wFMNsq + wFMN_NADH + wFMNH_NAD
    fC1 = C1_CONC / denf

    ## N3− + N2 = N3 + N2−
    v7 = kf7_C1 * (N3r_C1 * N2_C1 - N3_C1 * N2r_C1 * rKEQ7_C1)
    ## Q association
    q = Q_n * C1_INHIB
    v8 = kf8_C1 * (Iq_C1 * q - Q_C1 * rKEQ8_C1)
    ## CI.Q + N2− = CIQsq + N2
    v9 = kf9_C1 * (Q_C1 * N2r_C1 - SQ_C1 * N2_C1 * rKEQ9_C1)
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = kf13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 * rKEQ13_C1)
    ## QH2 dissociation
    qh2 = QH2_n * C1_INHIB
    v14 = kf14_C1 * (QH2_C1 - Iq_C1 * qh2 * rKEQ14_C1)
    ## Flavin site ROS generation
    v16 = kf16_C1 * E_LEAK_C1 * (FMNH * O2 - FMNsq * sox_m * rKEQ16_C1)
    ## Quinone site ROS generation
    v17 = kf17_C1 * (SQ_C1 * O2 - Q_C1 * sox_m * rKEQ17_C1)

    ## State transition rates in the quinone site
    ## 1 = Iq 2 = IqQ, 3 = IqSQ, 4 = IqQH2
    b12 = kf8_C1 * q
    b21 = kf8_C1 * rKEQ8_C1
    b23 = kf9_C1 * N2r_C1 + kf17_C1 * rKEQ17_C1 * sox_m
    b32 = kf9_C1 * rKEQ9_C1 * N2_C1 + kf17_C1 * O2
    b34 = kf13_C1 * N2r_C1 * fhm^2
    b43 = kf13_C1 * rKEQ13_C1 * N2_C1
    b41 = kf14_C1
    b14 = kf14_C1 * rKEQ14_C1 * qh2
    qDen = wIq + wIqQ + wIqSQ + wIqQH2
    qC1 = C1_CONC / qDen

    eqs = [
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        FMN ~ wFMN * fC1,
        FMN_NADH ~ wFMN_NADH * fC1,
        FMNH_NAD ~ wFMNH_NAD * fC1,
        FMN_NAD ~ wFMN_NAD * fC1,
        FMNH_NADH ~ wFMNH_NADH * fC1,
        FMNH ~ wFMNH * fC1,
        FMNsq ~ wFMNsq * fC1,
        ET_C1 ~ N2_C1 + N2r_C1,
        ET_C1 ~ N3_C1 + N3r_C1,
        N3_C1 ~ FMN * (FMN + FMNsq * KEQ6_C1),
        rKEQ13_C1 ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
        wIq ~ b21 * b32 * b41 + b21 * b32 * b43 + b21 * b34 * b41 + b23 * b34 * b41,
        wIqQ ~ b12 * b32 * b41 + b12 * b32 * b43 + b12 * b34 * b41 + b14 * b32 * b43,
        wIqSQ ~ b12 * b23 * b41 + b12 * b23 * b43 + b14 * b21 * b43 + b14 * b23 * b43,
        wIqQH2 ~ b12 * b23 * b34 + b14 * b21 * b32 + b14 * b21 * b34 + b14 * b23 * b34,
        Iq_C1 ~ wIq * qC1,
        Q_C1 ~ wIqQ * qC1,
        SQ_C1 ~ wIqSQ * qC1,
        QH2_C1 ~ wIqQH2 * qC1,
        vQC1 ~ -v8,
        vNADHC1 ~ -0.5 * (v7 + v12 + v16),
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROSC1 ~ vROSIf + vROSIq,
        vQH2C1 ~ v14,
        vHresC1 ~ 4 * v13,
        vNADC1 ~ -vNADHC1,
        TNC1 ~ vNADC1 / C1_CONC,
    ]
    return System(eqs, t; name)
end
