# Five state complex I model
function c1_5(; name=:c1f,
    Q_n=1800μM, QH2_n=200μM,
    nad=2500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.01μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0, MT_PROT=1)

    @parameters begin
        ET_C1 = 17μM                ## Activity of complex I
        KI_DOX_C1 = 400μM           ## DOX IC50 on complex I
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -150mV            ## B. taurus ISC N2 redox potential
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV
        KI_NADH_C1 = 50μM
        KD_NADH_C1 = 100μM
        KI_NAD_C1 = 1000μM
        KD_NAD_C1 = 25μM
        ## NADH + FMN = NAD+ + FMNH-
        KEQ_NADH_FMN = exp(2iVT * (Em_FMN_FMNH - Em_NAD))
        kf1_C1 = 83Hz / μM
        KEQ1_C1 = 0.01 / μM
        kb1_C1 = kf1_C1 / KEQ1_C1
        kf3_C1 = 1e6Hz
        KEQ3_C1 = KEQ_NADH_FMN / KEQ1_C1
        kb3_C1 = kf3_C1 / KEQ3_C1
        ## FMNsq + N2 = FMN + N2-
        KEQ_FMNsq_N2 = exp(iVT * (Em_N2 - Em_FMN_FMNsq))
        ## FMNH + N2 = FMNsq + N2-
        KEQ_FMNH_N2 = exp(iVT * (Em_N2 - Em_FMNsq_FMNH))
        ## N2- + Q = N2 + SQ
        KEQ_N2r_Q = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        ## FMNsq + Q = FMN + SQ
        KEQ_FMNsq_Q = exp(iVT * (Em_Q_SQ_C1 - Em_FMN_FMNsq))
        ## FMNH + Q = FMNsq + SQ
        KEQ_FMNH_Q = exp(iVT * (Em_Q_SQ_C1 - Em_FMNsq_FMNH))
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = inv(10μM)
        kb8_C1 = kf8_C1 / KEQ8_C1
        kf9_C1 = 10000Hz
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kb9_C1 = kf9_C1 / KEQ9_C1
        kf13_C1 = 2.7e6Hz
        kf14_C1 = 1000Hz
        KEQ14_C1 = 20μM
        kb14_C1 = kf14_C1 / KEQ14_C1
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kb16_C1 = kf16_C1 / KEQ16_C1
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        KEQ17_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
        kb17_C1 = kf17_C1 / KEQ17_C1
    end

    C1_CONC = ET_C1 * MT_PROT

    @variables begin
        C1_0(t) = C1_CONC
        FMN_N2_Q(t)
        C1_1(t) = 0
        FMNsq_N2_Q(t)
        FMN_N2r_Q(t)
        FMN_N2_SQ(t)
        C1_2(t) = 0
        FMNH_N2_Q(t)
        FMNsq_N2r_Q(t)
        FMNsq_N2_SQ(t)
        FMN_N2r_SQ(t)
        C1_3(t) = 0
        FMNH_N2r_Q(t)
        FMNsq_N2r_SQ(t)
        FMNH_N2_SQ(t)
        C1_4(t) ## Conserved
        FMNH_N2r_SQ(t)
        rKEQ_N2r_SQ(t)
        ## Reaction rates
        vQC1(t)
        vQH2C1(t)
        vROSC1(t)
        vROSIf(t)
        vROSIq(t)
        vNADHC1(t)
        vNADC1(t)
        TNC1(t)
        vHresC1(t)
    end

    ## QSSA for the A <--> B <--> C reaction
    _kfkb(kf1, kb1, kf2, kb2) = (kf1 * kf2 / (kf2 + kb1), kb2 * kb1 / (kf2 + kb1))

    C1_INHIB = (1 - ROTENONE_BLOCK) / (1 + (DOX / KI_DOX_C1)^3)
    ## Mitochondrial pH factor
    fhm = h_m * inv(1E-7Molar)

    ## Population of each state
    wFMNsq_N2_Q = 1 * fhm
    wFMN_N2r_Q = KEQ_FMNsq_N2
    wFMN_N2_SQ = KEQ_FMNsq_Q
    w1 = wFMNsq_N2_Q + wFMN_N2r_Q + wFMN_N2_SQ

    wFMNH_N2_Q = 1
    wFMNsq_N2r_Q = KEQ_FMNH_N2
    wFMNsq_N2_SQ = KEQ_FMNH_Q
    wFMN_N2r_SQ = wFMNsq_N2r_Q * KEQ_FMNsq_Q / fhm
    w2 = wFMNH_N2_Q + wFMNsq_N2r_Q + wFMNsq_N2_SQ + wFMN_N2r_SQ

    wFMNH_N2r_Q = 1
    wFMNsq_N2r_SQ = KEQ_FMNH_Q
    wFMNH_N2_SQ = KEQ_N2r_Q
    w3 = wFMNH_N2r_Q + wFMNsq_N2r_SQ + wFMNH_N2_SQ

    ## I(n) + NADH = I(n+2)H- + NAD+
    kfn, kbn = _kfkb(kf1_C1 * nadh, kb1_C1, kf3_C1, kb3_C1 * nad)
    v02_nadh = kfn * FMN_N2_Q - kbn * FMNH_N2_Q
    v13_nadh = kfn * (FMN_N2r_Q + FMN_N2_SQ) - kbn * (FMNH_N2r_Q + FMNH_N2_SQ)
    v24_nadh = kfn * FMN_N2r_SQ - kbn * FMNH_N2r_SQ

    ## I(n)-SQ + 6Hm = I(n-2) + QH2 + 4Hi; I(n-2) + Q = I(n-2)-Q
    kfq, kbq = _kfkb(kf13_C1 * fhm^2, kf13_C1 * rKEQ_N2r_SQ / KEQ14_C1 * QH2_n, kf8_C1 * Q_n, kb8_C1)
    v20_q = kfq * FMN_N2r_SQ - kbq * FMN_N2_Q
    v31_q = kfq * FMNsq_N2r_SQ - kbq * FMNsq_N2_Q
    v42_q = kfq * FMNH_N2r_SQ - kbq * FMNH_N2_Q

    ## Flavin site ROS production
    v21_f = kf16_C1 * FMNH_N2_Q * O2 - kb16_C1 * FMNsq_N2_Q * sox_m
    v32_f = kf16_C1 * (FMNH_N2r_Q + FMNH_N2_SQ) * O2 - kb16_C1 * (FMNsq_N2r_Q + FMNsq_N2_SQ) * sox_m
    v43_f = kf16_C1 * FMNH_N2r_SQ * O2 - kb16_C1 * FMNsq_N2r_SQ * sox_m
    v16 = v21_f + v32_f + v43_f
    ## Quinone site ROS production
    v10_q = kf17_C1 * FMN_N2_SQ * O2 - kb17_C1 * FMN_N2_Q * sox_m
    v21_q = kf17_C1 * (FMNsq_N2_SQ + FMN_N2r_SQ) * O2 - kb17_C1 * (FMNsq_N2_Q + FMN_N2r_Q) * sox_m
    v32_q = kf17_C1 * (FMNsq_N2r_SQ + FMNH_N2_SQ) * O2 - kb17_C1 * (FMNH_N2_Q + FMNsq_N2r_Q) * sox_m
    v43_q = kf17_C1 * FMNH_N2r_SQ * O2 - kb17_C1 * FMNH_N2r_Q * sox_m
    v17 = v10_q + v21_q + v32_q + v43_q

    eqs = [
        C1_CONC ~ C1_0 + C1_1 + C1_2 + C1_3 + C1_4,
        rKEQ_N2r_SQ ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
        FMN_N2_Q ~ C1_0,
        FMNsq_N2_Q ~ C1_1 / w1 * wFMNsq_N2_Q,
        FMN_N2r_Q ~ C1_1 / w1 * wFMN_N2r_Q,
        FMN_N2_SQ ~ C1_1 / w1 * wFMN_N2_SQ,
        FMNH_N2_Q ~ C1_2 / w2 * wFMNH_N2_Q,
        FMNsq_N2r_Q ~ C1_2 / w2 * wFMNsq_N2r_Q,
        FMNsq_N2_SQ ~ C1_2 / w2 * wFMNsq_N2_SQ,
        FMN_N2r_SQ ~ C1_2 / w2 * wFMN_N2r_SQ,
        FMNH_N2r_Q ~ C1_3 / w3 * wFMNH_N2r_Q,
        FMNsq_N2r_SQ ~ C1_3 / w3 * wFMNsq_N2r_SQ,
        FMNH_N2_SQ ~ C1_3 / w3 * wFMNH_N2_SQ,
        FMNH_N2r_SQ ~ C1_4,
        D(C1_0) ~ -v02_nadh + v20_q + v10_q,
        D(C1_1) ~ -v13_nadh + v31_q + v21_f - v10_q + v21_q,
        D(C1_2) ~ v02_nadh - v24_nadh - v20_q + v42_q - v21_f + v32_f - v21_q + v32_q,
        D(C1_3) ~ v13_nadh - v31_q - v32_f + v43_f - v32_q + v43_q,
        vNADHC1 ~ -(v02_nadh + v13_nadh + v24_nadh),
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROSC1 ~ vROSIf + vROSIq,
        vQH2C1 ~ v20_q + v31_q + v42_q,
        vQC1 ~ -vQH2C1,
        vHresC1 ~ 4 * vQH2C1,
        vNADC1 ~ -vNADHC1,
        TNC1 ~ vNADC1 / C1_CONC,
    ]
    return System(eqs, t; name)
end
