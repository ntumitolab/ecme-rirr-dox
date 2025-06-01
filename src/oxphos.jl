"Electron transport chain (ETC)"
function get_etc_sys(;
    DOX=0μM,                    # Doxorubicin concentration
    MT_PROT=1,                  # OXPHOS protein content scale factor
    O2=6μM,                     # Oxygen concentration
    h_i=exp10(-7) * Molar,      # IMS proton concentration
    h_m=exp10(-7.6) * Molar,    # Matrix proton concentration
    name=:etcsys,
    ROTENONE_BLOCK=0,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    CYANIDE_BLOCK=0,
    nad_m=nothing,
    nadh_m=nothing,
    dpsi=nothing,
    sox_m=nothing,
    suc=nothing,
    fum=nothing,
    oaa=nothing)

    # Define missing variables
    isnothing(nad_m) && @variables nad_m(t)
    isnothing(nadh_m) && @variables nadh_m(t)
    isnothing(dpsi) && @variables dpsi(t)
    isnothing(sox_m) && @variables sox_m(t)
    isnothing(suc) && @variables suc(t)
    isnothing(fum) && @variables fum(t)
    isnothing(oaa) && @variables oaa(t)

    # Q cycle variables
    @variables begin
        Q_n(t)
        QH2_n(t)
        QH2_p(t)
        Q_p(t) ## Conserved
        SQn(t)
    end

    # Complex I using a simplified Markevich model
    # Rapid equlibrium in the flavin site
    # QSSA for the catalytic cycle in the quinone site
    @parameters begin
        ET_C1 = 1μM              ## Activity of complex I
        KI_DOX_C1 = 400μM         ## DOX IC50 on complex I
        K_RC_DOX = 1000 / 15mM    ## DOX redox cycling constant
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## 800mV (?) in Markevich, 2015
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
        kf7_C1 = 10000Hz / μM
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
        SQ_C1(t)
        QH2_C1(t)
        rKEQ_N2r_SQ(t)
        ## Reaction rates
        TNC1(t)
        vQC1(t)
        vQH2C1(t)
        vNADHC1(t)
        vNADC1(t)
        vROSIf(t)
        vROSIq(t)
        vROSC1(t)
        vHresC1(t)
    end

    c1eqs = let
        C1_CONC = ET_C1 * MT_PROT
        ## complex I inhibition by DOX and rotenone
        C1_INHIB = hil(KI_DOX_C1, DOX, 3) * (1 - ROTENONE_BLOCK)
        n2 = C1_INHIB * N2_C1
        n2r = C1_INHIB * N2r_C1
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
        ## First electron transfer
        v7 = kf7_C1 * C1_INHIB * (N2r_C1 * Q_C1 - N2_C1 * SQ_C1 * rKEQ_N2r_Q)
        ## Second electron transfer
        v13 = kf13_C1 * C1_INHIB * (N2r_C1 * SQ_C1 * fhm^2 - N2_C1 * QH2_C1 * rKEQ_N2r_SQ)
        ## Q binding and QH2 unbinding
        q = Q_n * rKD_Q_C1
        qh2 = QH2_n * rKD_QH2_C1
        v14 = kf14_C1 * (QH2_C1 * q - Q_C1 * qh2)
        ## Flavin site ROS generation
        v16 = kf16_C1 * E_LEAK_C1 * (FMNH * O2 - FMNsq * sox_m * rKEQ16_C1)
        ## Quinone site ROS generation
        v17 = kf17_C1 * (SQ_C1 * O2 - Q_C1 * sox_m * rKEQ17_C1)

        ## State transition rates in the quinone site
        ## 1 = IqQ, 2 = IqSQ, 3 = IqQH2
        ## First electron transfer
        b12a = kf7_C1 * n2r
        b21a = kf7_C1 * rKEQ_N2r_Q * n2
        v7 = b12a * Q_C1 - b21a * SQ_C1
        ## Quinone site ROS generation
        b21b = kf17_C1 * O2
        b12b = kf17_C1 * rKEQ17_C1 * sox_m
        v17 = b21b * SQ_C1 - b12b * Q_C1
        b12 = b12a + b12b
        b21 = b21a + b21b
        ## Second electron transfer
        b23 = kf13_C1 * n2r * fhm^2
        b32 = kf13_C1 * rKEQ_N2r_SQ * n2
        v13 = b23 * SQ_C1 - b32 * QH2_C1
        ## Q binding and QH2 unbinding
        q = Q_n * rKD_Q_C1
        qh2 = QH2_n * rKD_QH2_C1
        b31 = kf14_C1 * q
        b13 = kf14_C1 * qh2
        v14 = b31 * QH2_C1 - b13 * Q_C1

        w1 = b21 * (b31 + b32) + b23 * b31
        w2 = b12 * (b31 + b32) + b13 * b32
        w3 = b12 * b23 + b13 * (b21 + b23)
        qDen = w1 + w2 + w3
        qC1 = C1_CONC / qDen

        eqs = [
            rKEQ_N2r_SQ ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
            FMN ~ wFMN * fC1,
            FMN_NAD ~ wFMN_NAD * fC1,
            FMNH ~ wFMNH * fC1,
            FMNsq ~ wFMNsq * fC1,
            FMNH_NADH ~ wFMNH_NADH * fC1,
            FMN_NADH ~ wFMN_NADH * fC1,
            FMNH_NAD ~ wFMNH_NAD * fC1,
            N2_C1 ~ FMNsq / (FMNsq + FMNH * KEQ_FMNH_N2),
            N2r_C1 ~ 1 - N2_C1,
            Q_C1 ~ w1 * qC1,
            SQ_C1 ~ w2 * qC1,
            QH2_C1 ~ C1_CONC - Q_C1 - SQ_C1,
            vQC1 ~ -v14,
            vNADHC1 ~ -0.5 * (v7 + v13 + v16),
            vROSIf ~ v16,
            vROSIq ~ v17,
            vROSC1 ~ vROSIf + vROSIq,
            vQH2C1 ~ v14,
            vHresC1 ~ 4 * v13,
            vNADC1 ~ -vNADHC1,
            TNC1 ~ vNADC1 / C1_CONC,
        ]
    end

    # Reversible complex II (SDH)
    @parameters begin
        KI_DOX_C2 = 2000μM # DOX inhibition concentration (IC50) on complex II
        K_C2 = 250 / (minute * mM)  # Reaction rate constant of SDH (complex II)
        KI_OAA_C2 = 150μM           # Inhibition constant for OAA
        Em_FUM_SUC = 40mV           # midpoint potential of FUM -> SUC
        Em_Q_QH2 = 100mV            # midpoint potential of Q -> QH2
        rKEQ_C2 = exp(-2iVT * (Em_Q_QH2 - Em_FUM_SUC)) # (Reverse) equlibrium constant of SDH
    end

    @variables vSDH(t)
    kc2 = K_C2 * hil(KI_OAA_C2, oaa) * hil(KI_DOX_C2, DOX, 3)
    c2eqs = [vSDH ~ kc2 * (Q_n * suc - QH2_n * fum * rKEQ_C2)]

    # complex IV (CCO)
    @parameters begin
        rhoC4 = 325μM
        δ₅ = 0.5
        K34_C4 = 2.9445e10Hz / mM^3 # @pH7
        K43_C4 = 2.9E-6Hz / mM^3
        K35_C4 = 750Hz / mM
        K36_C4 = 4.826e11Hz / mM
        K63_C4 = 4.826Hz / mM
        K37_C4 = 2.92367e6Hz
        K73_C4 = 0.029236Hz # @pH7
        KI_DOX_C4 = 165μM   # DOX inhibition concentration (IC50) on complex IV
    end

    @variables begin
        cytc_rd(t)
        cytc_ox(t)
        vO2(t)
        vHresC4(t)
        C4_Y(t)
        C4_Yr(t)
        C4_YO(t)
        C4_YOH(t)
    end

    c4eqs = let
        C4_INHIB = hil(KI_DOX_C4, DOX, 3) * (1 - CYANIDE_BLOCK)  # complex IV inhibition by DOX
        C4_CONC = rhoC4 * MT_PROT
        aδ = exp(-iVT * δ₅ * dpsi)
        a1mδ = exp(iVT * (1 - δ₅) * dpsi)
        f_hm = h_m * inv(1E-7Molar) * aδ
        f_hi = h_i * inv(1E-7Molar) * a1mδ
        f_cr = cytc_rd
        f_co = cytc_ox * a1mδ
        a12 = K34_C4 * f_cr^3 * f_hm^4  # a12 = K34 * exp(-δ₅ * 4 * vfrt) * cytc_rd^3 * h_m^4
        a21 = K43_C4 * f_co^3 * f_hi    # K43 * exp((1 - δ₅) * 4 * vfrt) * cytc_ox^3 * h_i
        a23 = C4_INHIB * K35_C4 * O2
        a32 = 0
        a34 = K36_C4 * f_cr * f_hm^3    # K36 * exp(-δ₅ * 3 * vfrt) * cytc_rd * h_m^3
        a43 = K63_C4 * f_co * f_hi^2    # K63 * exp((1 - δ₅) * 3 * vfrt) * cytc_ox * h_i^2
        a41 = K37_C4 * f_hm  # K37 * exp(-δ₅ * vfrt) * h_m
        a14 = K73_C4 * f_hi  # K73 * exp((1 - δ₅) * F_RT * dpsi) * h_i

        # Weight of each state (from KA pattern)
        C4_e1 = a41 * a34 * (a21 + a23)
        C4_e2 = a12 * a41 * a34
        C4_e3 = a23 * a12 * (a41 + a43) + a43 * a14 * (a21 + a23)
        C4_e4 = a34 * (a14 * (a21 + a23) + a23 * a12)
        den = C4_e1 + C4_e2 + C4_e3 + C4_e4
        # Reaction rates
        # v34 = C4_CONC * (C4_Y * a12 - C4_Yr * a21)
        v35 = C4_CONC * C4_Yr * a23
        # v36 = C4_CONC * (a34 * C4_YO - a43 * C4_YOH)
        # v37 = C4_CONC * (a41 * C4_YOH - a14 * C4_Y)

        c4eqs = [
            C4_CONC ~ cytc_rd + cytc_ox,
            C4_Y ~ C4_e1 / den,
            C4_Yr ~ C4_e2 / den,
            C4_YO ~ C4_e3 / den,
            C4_YOH ~ C4_e4 / den,
            vO2 ~ v35,
            vHresC4 ~ 8vO2,
        ]
    end

    # complex III and the Q cycle
    @parameters begin
        rhoC3 = 325μM
        Q_T = 4mM        # Total CoQ pool
        KI_DOX_C3 = 185μM  # DOX inhibition concentration (IC50) on complex III
        EmQ_C3 = +73mV
        EmQn_SQn = +70mV
        EmSQn_QH2n = +170mV
        EmbL_bHo = -40mV
        EmbL_bHr = EmbL_bHo - 60mV
        EmbH_bLo = +40mV
        EmbH_bLr = EmbH_bLo - 60mV
        EmFeS = +280mV
        Emcytc1 = +240mV
        Emcytc = +260mV
        ## Split of electrical potentials
        δ₁_C3 = 0.5
        δ₂_C3 = 0.5
        δ₃_C3 = 0.5
        ## Split of the electrical distance across the IMM
        α_C3 = 0.25
        β_C3 = 0.5
        γ_C3 = 0.25
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        K04_C3 = 70Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmFeS + EmbL_bHo - 2EmQ_C3))
        KEQ4_RD_C3 = exp(iVT * (EmFeS + EmbL_bHr - 2EmQ_C3))
        ## Q_p = Q_n; QH2_p = QH2_n
        KD_Q = 22000Hz
        ## bL- + bH = bL + bH-
        K06_C3 = 166.67Hz
        KEQ6_C3 = exp(iVT * (EmbH_bLo - EmbL_bHo)) ## +70mV
        ## bH- + Q = bH + Q-
        K07_OX_C3 = 13.33Hz / mM
        K07_RD_C3 = 1.67Hz / mM
        KEQ7_OX_C3 = exp(iVT * (EmQn_SQn - EmbH_bLo)) ## +30mV
        KEQ7_RD_C3 = exp(iVT * (EmQn_SQn - EmbH_bLr)) ## +90mV
        ## bH- + Q- + 2H+ = bH + QH2
        K08_OX_C3 = 83.33Hz / mM
        K08_RD_C3 = 8.33Hz / mM
        KEQ8_OX_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLo)) ## +130mV
        KEQ8_RD_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLr)) ## +190mV
        ## FeS- + c1_3+ = FeS + c1_2+
        K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))  ## -40mV
        ## bL- + O2 + Q = bL + O2- + Q
        K010_C3 = 1400Hz / mM
        KEQ10_OX_C3 = exp(iVT * (Em_O2_SOX - EmbL_bHo)) ## -130mV
        KEQ10_RD_C3 = exp(iVT * (Em_O2_SOX - EmbL_bHr)) ## -70mV
        ## c1_2+ + c_3+ = c1_3+ + c_2+
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1)) ## +20mV
    end

    C3_CONC = rhoC3 * MT_PROT

    @variables begin
        fes_ox(t) = C3_CONC
        fes_rd(t) # Conserved
        cytc1_ox(t) = C3_CONC
        cytc1_rd(t) # Conserved
        blo_bho(t) = C3_CONC
        blr_bho(t) = 0
        blo_bhr(t) = 0
        blr_bhr(t) ## Conserved
        fracbLrd(t)
        fracbHrd(t)
        vROSC3(t)
        vHresC3(t)
        vHres(t)
        vROS(t)
    end

    c3eqs = let
        fHi = h_i * inv(1E-7Molar)
        fHm = h_m * inv(1E-7Molar)
        # complex III inhibition by DOX and antimycin
        C3_INHIB = hil(KI_DOX_C3, DOX, 3) * (1 - ANTIMYCIN_BLOCK)
        # Q reduction
        v1 = vQH2C1 + vSDH
        # QH2 diffusion
        v2 = KD_Q * (QH2_n - QH2_p)
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        qh2p = QH2_p * (1 - MYXOTHIAZOLE_BLOCK)
        qp = Q_p * (1 - MYXOTHIAZOLE_BLOCK)
        FeS = fes_ox / (fes_ox + fes_rd)
        FeSm = fes_rd / (fes_ox + fes_rd)
        el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
        er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
        k4ox = K04_C3 * KEQ4_OX_C3 * el4
        k4rd = K04_C3 * KEQ4_RD_C3 * el4
        km4 = K04_C3 * er4 * fHi^2
        v4ox = k4ox * qh2p * FeS * blo_bho - km4 * qp * FeSm * blr_bho
        v4rd = k4rd * qh2p * FeS * blo_bhr - km4 * qp * FeSm * blr_bhr
        # v5 = Q diffusion (p-side -> n-side)
        v5 = KD_Q * (Q_p - Q_n)
        ## bL- + bH = bL + bH-
        el6 = exp(-iVT * β_C3 * δ₂_C3 * dpsi)
        er6 = exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi)
        k6 = K06_C3 * KEQ6_C3 * el6
        km6 = K06_C3 * er6
        v6 = k6 * blr_bho - km6 * blo_bhr
        ## bH- + Q = bH + Q-
        ## bH- + Q- + 2Hi = bH + QH2
        Qi_avail = (C3_CONC - SQn) / C3_CONC * C3_INHIB
        el7 = exp(-iVT * γ_C3 * δ₃_C3 * dpsi)
        er7 = exp(iVT * γ_C3 * (1 - δ₃_C3) * dpsi)
        qn = Q_n * Qi_avail
        qh2n = QH2_n * Qi_avail
        k7ox = K07_OX_C3 * KEQ7_OX_C3 * el7
        k7rd = K07_RD_C3 * KEQ7_RD_C3 * el7
        km7ox = K07_OX_C3 * er7
        km7rd = K07_RD_C3 * er7
        k8ox = K08_OX_C3 * KEQ8_OX_C3 * el7 * fHm^2
        k8rd = K08_RD_C3 * KEQ8_RD_C3 * el7 * fHm^2
        km8ox = K08_OX_C3 * er7
        km8rd = K08_RD_C3 * er7
        v7ox = k7ox * blo_bhr * qn - km7ox * blo_bho * SQn
        v7rd = k7rd * blr_bhr * qn - km7rd * blr_bho * SQn
        v8ox = k8ox * blo_bhr * SQn - km8ox * blo_bho * qh2n
        v8rd = k8rd * blr_bhr * SQn - km8rd * blr_bho * qh2n
        ## FeS- + c1_3+ = FeS + c1_2+
        v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
        ## cytc1_2+ + cytc_3+ = cytc1_3+  + cytc_2+
        v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)
        ## bL- + O2 + Qp = bL + O2- + Qp
        k10ox = K010_C3 * KEQ10_OX_C3 * er4
        k10rd = K010_C3 * KEQ10_RD_C3 * er4
        km10 = K010_C3 * el4
        fq = Q_p / (Q_p + QH2_p)
        v10ox = fq * (k10ox * O2 * blr_bho - km10 * sox_m * blo_bho)
        v10rd = fq * (k10rd * O2 * blr_bhr - km10 * sox_m * blo_bhr)

        eqs = [
            C3_CONC ~ blo_bho + blr_bho + blo_bhr + blr_bhr,
            C3_CONC ~ fes_ox + fes_rd,
            C3_CONC ~ cytc1_ox + cytc1_rd,
            Q_T ~ Q_n + SQn + QH2_n + QH2_p + Q_p,
            D(Q_n) ~ v5 - v7ox - v7rd - v1,
            D(SQn) ~ v7ox + v7rd - v8ox - v8rd,
            D(QH2_n) ~ v8ox + v8rd - v2 + v1,
            D(QH2_p) ~ v2 - v4ox - v4rd,
            # D(Q_p) ~ v10 + v4_ox + v4_rd - v5,
            D(blo_bho) ~ v7ox + v8ox - v4ox,
            D(blr_bho) ~ v4ox + v7rd + v8rd - v6,
            D(blo_bhr) ~ v6 - v4rd - v7ox - v8ox,
            ## D(blr_bhr) = v4rd - v7rd - v8rd
            D(fes_ox) ~ v9 - v4ox - v4rd,
            D(cytc1_ox) ~ v33 - v9,
            D(cytc_ox) ~ 4 * vO2 - v33,
            vHresC3 ~ v4ox + v4rd,
            vHres ~ vHresC1 + vHresC3 + vHresC4,
            vROSC3 ~ v10ox + v10rd,
            vROS ~ vROSC3 + vROSC1,
        ]
    end

    return System([c1eqs; c2eqs; c4eqs; c3eqs], t; name)
end

function get_c5_sys(; dpsi, h_i, h_m, atp_i, adp_i, atp_m, adp_m, pi_m, MT_PROT=1, C5_INHIB=1, use_mg=false, mg_i=1mM, mg_m=0.4mM, name=:c5sys)
    @parameters begin
        ρF1 = 5.0mM                 # Concentration of ATP synthase
        P1_C5 = 1.346E-8
        P2_C5 = 7.739E-7
        P3_C5 = 6.65E-15
        PA_C5 = 1.656E-5Hz
        PB_C5 = 3.373E-7Hz
        PC1_C5 = 9.651E-14Hz
        PC2_C5 = 4.585E-19Hz        # Magnus model
        # Equilibrium constant of ATP synthase (ΔG=-30kJ/mol)
        # 1.71E6mM in Magnus model was inconsistent to Caplan's model (1.71E6 Molar)
        KEQ_C5 = 2E5Molar
        G_H_MITO = 2E-6mM / ms / mV # Proton leak rate constant
        VMAX_ANT = 5E-3mM / ms      # Max rate of ANT, (Wei, 2011)
        H_ANT = 0.5  # Voltage steepness
    end

    @variables begin
        AF1(t)      # Relative electrochemical activity of ATP + H2O <-> ADP + Pi
        vC5(t)      # ATP synthesis rate
        vHu(t)      # Porton flux via ATP synthase
        vANT(t)     # ANT reaction rate
        vHleak(t)   # Proton leak rate
        ΔμH(t)      # proton motive force
        E_PHOS(t)   # phosphorylation potential
    end

    # F1-Fo ATPase
    if use_mg
        adp3_m, hadp_m, _, poly_adp_m = breakdown_adp(adp_m, h_m, mg_m)
        _, _, mgatp_m, poly_atp_m = breakdown_atp(atp_m, h_m, mg_m)
        poly_pi_m = pipoly(h_m)
        ke_app = KEQ_C5 * poly_atp_m / (poly_adp_m * poly_pi_m)
        v_af1 = ke_app * mgatp_m / (pi_m * (hadp_m + adp3_m))
    else
        v_af1 = KEQ_C5 * atp_m / (pi_m * adp_m)
    end

    vb = exp(iVT * 3 * 50mV) # Boundary potential
    vh = exp(iVT * 3 * dpsi)
    fh = (h_m / h_i)^3
    common = -ρF1 * MT_PROT * C5_INHIB / ((1 + P1_C5 * AF1) * vb + (P2_C5 + P3_C5 * AF1) * vh)
    v_c5 = common * ((PA_C5 * fh + PC1_C5 * vb) * AF1 - (PA_C5 + PC2_C5 * AF1) * vh)
    v_hu = 3 * common * (PA_C5 * fh * AF1 - vh * (PA_C5 + PB_C5))

    # Adenine nucleotide translocator (ANT)
    # Free adenylates
    if use_mg
        atp4_i = atp_i / (1 + h_i / KA_ATP + mg_i / KMG_ATP)
        atp4_m = atp_m / (1 + h_m / KA_ATP + mg_m / KMG_ATP)
        adp3_i = adp_i / (1 + h_i / KA_ADP + mg_i / KMG_ADP)
        adp3_m = adp_m / (1 + h_m / KA_ADP + mg_m / KMG_ADP)
    else
        atp4_i = 0.25 * atp_i
        atp4_m = 0.025 * atp_m
        adp3_i = 0.45 * adp_i
        adp3_m = 0.17 * adp_m
    end

    f_i = atp4_i / adp3_i
    f_m = adp3_m / atp4_m
    v_ant = VMAX_ANT * (1 - f_i * f_m * exp(-iVT * dpsi)) / ((1 + f_i * exp(-iVT * H_ANT * dpsi)) * (1 + f_m))

    eqs = [
        ΔμH ~ dpsi + nernst(h_i, h_m, 1),
        E_PHOS ~ VT * NaNMath.log(AF1),
        vANT ~ v_ant,
        AF1 ~ v_af1,
        vHleak ~ G_H_MITO * ΔμH,
        vC5 ~ v_c5,
        vHu ~ v_hu,
    ]
    return System(eqs, t; name)
end
