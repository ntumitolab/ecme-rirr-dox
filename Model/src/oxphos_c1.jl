_redoxcycling_c1(DOX=0μM, K_RC_DOX=1000 / 15mM) = 1 + K_RC_DOX * DOX
_inhibit_c1(DOX=0μM, KI_DOX_C1=400μM) = hil(KI_DOX_C1, DOX, 3)

"""
Complex I using a simplified Markevich model

- Rapid equlibrium in the flavin site
- QSSA for the catalytic cycle in the quinone site
"""
function get_c1_eqs(;
    Q_n,                    ## Ubiquinone concentration
    QH2_n,                  ## Ubiquinol concentration
    nad_m,                  ## NAD concentration
    nadh_m,                 ## NADH concentration
    dpsi,                   ## Mitochondrial membrane potential
    sox_m,                  ## Superoxide concentration in the mitochondrial matrix
    O2=6μM,                     ## Oxygen concentration
    h_i=exp10(-7) * Molar,      ## IMS proton concentration
    h_m=exp10(-7.6) * Molar,    ## Matrix proton concentration
    REDOX_CYCLING=_redoxcycling_c1(),
    C1_INHIB_DOX=_inhibit_c1(), ## Complex I inhibition by DOX and rotenone
    ROTENONE_BLOCK=0,
    C1_CONC = 17μM                ## Activity of complex I
    )

    @parameters begin
        # ET_C1 = 17μM            ## Activity of complex I
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -150mV             ## -150mV in B. taurus mitochondrial complex I; -80 mV in bacteria
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## +800mV in Markevich, 2015
        KI_NADH_C1 = 50μM
        KD_NADH_C1 = 100μM
        KI_NAD_C1 = 1000μM
        KD_NAD_C1 = 25μM
        ## NADH + FMN = NAD+ + FMNH-
        KEQ_NADH_FMN = exp(2iVT * (Em_FMN_FMNH - Em_NAD))
        ## 2FMNsq = (N1a) = FMN + FMNH- + H+
        rKEQ_FMNsq_Dis = exp(-iVT * (Em_FMNsq_FMNH - Em_FMN_FMNsq))
        ## FMNH- + N3 = FMNsq + N3-
        KEQ_FMNH_N3 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
        ## N3- + N2 = N3 + N2-
        kf7_C1 = 10000Hz / μM
        KEQ7_C1 = exp(iVT * (Em_N2 - Em_N3))
        kr7_C1 = kf7_C1 / KEQ7_C1
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = inv(10μM)
        kr8_C1 = kf8_C1 / KEQ8_C1
        kf9_C1 = 4E5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kr9_C1 = kf9_C1 / KEQ9_C1
        kf13_C1 = 2.7e6Hz / μM
        kf14_C1 = 1000Hz
        KEQ14_C1 = 20μM
        kr14_C1 = kf14_C1 / KEQ14_C1
        ## SOX production from If site
        kf16_C1 = 2Hz / μM
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kr16_C1 = kf16_C1 / KEQ16_C1
        ## SOX production from Iq site
        kf17_C1 = 0.02Hz / μM
        KEQ17_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
        kr17_C1 = kf17_C1 / KEQ17_C1
    end

    sts = @variables N2r_C1(t) = 0

    @variables begin
        ## Flavin site
        FMN(t)
        FMN_NAD(t)
        FMNsq(t)
        FMNH(t)
        FMNH_NADH(t)
        FMN_NADH(t)
        FMNH_NAD(t)
        ## FeS cluster
        N3_C1(t)
        N3r_C1(t)
        N2_C1(t)
        ## Quinone site
        C1(t)
        Q_C1(t)
        SQ_C1(t)
        QH2_C1(t)
        rKEQ_N2r_SQ(t)
        ## Reaction rates
        TNC1(t)  ## Turn over rates
        vQC1(t)
        vQH2C1(t)
        vNADHC1(t)
        vNADC1(t)
        vROSIf(t)
        vROSIq(t)
        vROSC1(t)
        vHresC1(t)
    end

    C1_INHIB = C1_INHIB_DOX * (1 - ROTENONE_BLOCK)
    ## Mitochondrial pH factor
    fhm = h_m * inv(1E-7Molar)
    ## Flavin site in rapid equilibrium
    ## Weights in the flavin site
    wFMN = 1
    wFMN_NAD = wFMN * nad_m / KI_NAD_C1
    wFMN_NADH = wFMN * nadh_m / KD_NADH_C1
    wFMNH = wFMN * (nadh_m / nad_m) * KEQ_NADH_FMN
    wFMNH_NAD = wFMNH * nad_m / KD_NAD_C1
    wFMNH_NADH = wFMNH * nadh_m / KI_NADH_C1
    wFMNsq = NaNMath.sqrt(wFMN * wFMNH * rKEQ_FMNsq_Dis * fhm)
    fDen = wFMN + wFMN_NAD + wFMNH + wFMNH_NADH + wFMNsq + wFMN_NADH + wFMNH_NAD
    fC1 = C1_CONC / fDen
    ## FMNH + O2 = FMNsq + sox
    v16 = REDOX_CYCLING * (kf16_C1 * FMNH * O2 - kr16_C1 * FMNsq * sox_m)

    ## N3− + N2 = N3 + N2−
    v7 = kf7_C1 * N3r_C1 * N2_C1 - kr7_C1 * N3_C1 * N2r_C1
    v12 = v7

    ## Quinone site state transition rates
    ## C1 + Q = Q_C1
    b12 = kf8_C1 * Q_n * C1_INHIB
    b21 = kr8_C1
    v8 = b12 * C1 - b21 * Q_C1
    ## Q_C1 + N2r = SQ_C1 + N2
    b23a = kf9_C1 * N2r_C1
    b32a = kr9_C1 * N2_C1
    v9 = b23a * Q_C1 - b32a * SQ_C1
    ## C1_SQ + N2r + 6Hm = C1_QH2 + N2 + 4Hi
    b34 = kf13_C1 * N2r_C1 * fhm^2
    b43 = kf13_C1 * rKEQ_N2r_SQ * N2_C1
    v13 = b34 * SQ_C1 - b43 * QH2_C1
    ## C1_QH2 = C1 + QH2
    b41 = kf14_C1
    b14 = kr14_C1 * QH2_n * C1_INHIB
    v14 = b41 * QH2_C1 - b14 * C1
    ## C1_SQ + O2 = C1_Q + sox
    b32b = kf17_C1 * O2
    b23b = kr17_C1 * sox_m
    v17 = b32b * SQ_C1 - b23b * Q_C1
    b23 = b23a + b23b
    b32 = b32a + b32b
    ## KA pattern
    wC1 = b21 * b32 * b41 + b21 * b32 * b43 + b21 * b34 * b41 + b23 * b34 * b41
    wC1_Q = b12 * b32 * b41 + b12 * b32 * b43 + b12 * b34 * b41 + b14 * b32 * b43
    wC1_SQ = b12 * b23 * b41 + b12 * b23 * b43 + b14 * b21 * b43 + b14 * b23 * b43
    wC1_QH2 = b12 * b23 * b34 + b14 * b21 * b32 + b14 * b21 * b34 + b14 * b23 * b34
    qDen = wC1 + wC1_Q + wC1_SQ + wC1_QH2
    qC1 = C1_CONC / qDen

    eqs_c1 = [
        rKEQ_N2r_SQ ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
        FMN ~ wFMN * fC1,
        FMN_NAD ~ wFMN_NAD * fC1,
        FMNH ~ wFMNH * fC1,
        FMNsq ~ wFMNsq * fC1,
        FMNH_NADH ~ wFMNH_NADH * fC1,
        FMN_NADH ~ wFMN_NADH * fC1,
        FMNH_NAD ~ wFMNH_NAD * fC1,
        N3_C1 ~ C1_CONC * FMNsq / (FMNsq + FMNH * KEQ_FMNH_N3),
        C1_CONC ~ N3_C1 + N3r_C1,
        C1_CONC ~ N2_C1 + N2r_C1,
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        C1 ~ wC1 * qC1,
        Q_C1 ~ wC1_Q * qC1,
        SQ_C1 ~ wC1_SQ * qC1,
        QH2_C1 ~ wC1_QH2 * qC1,
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

    return (; eqs_c1, vQC1, vNADHC1, vROSIf, vROSIq, vROSC1, vQH2C1, vHresC1, vNADC1, TNC1)
end
