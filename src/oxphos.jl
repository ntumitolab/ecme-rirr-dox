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

    # Define variables when needed
    isnothing(nad_m) && @variables nad_m(t)
    isnothing(nadh_m) && @variables nadh_m(t)
    isnothing(dpsi) && @variables dpsi(t)
    isnothing(sox_m) && @variables sox_m(t)
    isnothing(suc) && @variables suc(t)
    isnothing(fum) && @variables fum(t)
    isnothing(oaa) && @variables oaa(t)

    # Q cycle variables
    @variables begin
        Q_n(t) = 1805μM
        QH2_n(t) = 123μM
        QH2_p(t) = 123μM
        Q_p(t) ## Conserved
        SQn(t) = 142μM
        SQp(t) = 0μM
    end

    # Complex I using a simplified Markevich model
    # Rapid equlibrium in the flavin site
    # QSSA for the catalytic cycle in the quinone site
    @parameters begin
        K_RC_DOX = 1000 / 15mM    ## DOX redox cycling constant
        ET_C1 = 17μM              ## Activity of complex I
        KI_DOX_C1 = 400μM         ## DOX IC50 on complex I
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## ~800mV (?) in Markevich, 2015
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
        rKEQ7_C1 = exp(-iVT * (Em_N2 - Em_N3))
        kb7_C1 = kf7_C1 * rKEQ7_C1
        kf8_C1 = 10Hz / μM
        rKEQ8_C1 = 10μM
        kb8_C1 = kf8_C1 * rKEQ8_C1
        kf9_C1 = 4E5Hz / μM
        rKEQ9_C1 = exp(-iVT * (Em_Q_SQ_C1 - Em_N2))
        kb9_C1 = kf9_C1 * rKEQ9_C1
        kf13_C1 = 2.7e6Hz / μM
        kf14_C1 = 1000Hz
        rKEQ14_C1 = inv(20μM)
        kb14_C1 = kf14_C1 * rKEQ14_C1
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        rKEQ16_C1 = exp(-iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kb16_C1 = kf16_C1 * rKEQ16_C1
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        rKEQ17_C1 = exp(-iVT * (Em_O2_SOX - Em_Q_SQ_C1))
        kb17_C1 = kf17_C1 * rKEQ17_C1
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
        ## FeS cluster
        N3_C1(t)
        N3r_C1(t)
        N2_C1(t)
        N2r_C1(t) = 0
        ## Quinone site
        C1(t)
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
        ## Electron leak scaling factor from complex I
        E_LEAK_C1 = 1 + K_RC_DOX * DOX
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
        v16 = kf16_C1 * FMNH * O2 - kb16_C1 * FMNsq * sox_m

        ## N3− + N2 = N3 + N2−
        v7 = kf7_C1 * N3r_C1 * N2_C1 - kb7_C1 * N3_C1 * N2r_C1
        v12 = v7

        ## Quinone site state transition rates
        ## C1 + Q = Q_C1
        b12 = kf8_C1 * Q_n * C1_INHIB
        b21 = kb8_C1
        v8 = b12 * C1 - b21 * Q_C1
        ## Q_C1 + N2r = SQ_C1 + N2
        b23a = kf9_C1 * N2r_C1
        b32a = kb9_C1 * N2_C1
        v9 = b23a * Q_C1 - b32a * SQ_C1
        ## C1_SQ + N2r + 6Hm = C1_QH2 + N2 + 4Hi
        b34 = kf13_C1 * N2r_C1 * fhm^2
        b43 = kf13_C1 * rKEQ_N2r_SQ * N2_C1
        v13 = b34 * SQ_C1 - b43 * QH2_C1
        ## C1_QH2 = C1 + QH2
        b41 = kf14_C1
        b14 = kb14_C1 * QH2_n * C1_INHIB
        v14 = b41 * QH2_C1 - b14 * C1
        ## C1_SQ + O2 = C1_Q + sox
        b32b = kf17_C1 * O2
        b23b = kb17_C1 * sox_m
        v17 = b32b * SQ_C1 - b23b * Q_C1
        b23 = b23a + b23b
        b32 = b32a + b32b

        ## KA pattern
        wC1 = b21*b32*b41 + b21*b32*b43 + b21*b34*b41 + b23*b34*b41
        wC1_Q = b12*b32*b41 + b12*b32*b43 + b12*b34*b41 + b14*b32*b43
        wC1_SQ = b12*b23*b41 + b12*b23*b43 + b14*b21*b43 + b14*b23*b43
        wC1_QH2 = b12*b23*b34 + b14*b21*b32 + b14*b21*b34 + b14*b23*b34
        qDen = wC1 + wC1_Q + wC1_SQ + wC1_QH2
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

    C4_CONC = rhoC4 * MT_PROT

    @variables begin
        cytc_rd(t) ## Conserved
        cytc_ox(t) = C4_CONC
        vO2(t)
        vHresC4(t)
        C4_Y(t)
        C4_Yr(t)
        C4_YO(t)
        C4_YOH(t)
    end

    c4eqs = let
        C4_INHIB = hil(KI_DOX_C4, DOX, 3) * (1 - CYANIDE_BLOCK)  # complex IV inhibition by DOX
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

    ## Semireverse bc1 complex model adapted from Gauthier, 2013
    @parameters begin
        rhoC3 = 325μM    ## Complex III activity
        KI_DOX_C3 = 185μM  ## DOX inhibition concentration (IC50) on complex III
        Q_T = 4mM        ## Total CoQ pool
        EmQ_C3 = +60mV   ## Ubiquinone redox potential at complex III Qo
        EmSQp_QH2p = +290mV
        EmQp_SQp = -170mV
        EmQn_SQn = +50mV
        EmSQn_QH2n = +150mV
        EmbL_bHo = -40mV
        EmbL_bHr = EmbL_bHo - 60mV
        EmbH_bLo = +20mV
        EmbH_bLr = EmbH_bLo - 60mV
        EmFeS = +280mV
        Emcytc1 = +245mV
        EmO2 = -160mV
        Emcytc = +255mV
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmFeS + EmbL_bHo - 2EmQ_C3))
        KEQ4_RD_C3 = exp(iVT * (EmFeS + EmbL_bHr - 2EmQ_C3))
        ## QH2 and Q flip rate in the IMM
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
        ## bL- + Q = bL + Q-
        K010_C3 = 28.33Hz / mM
        KEQ10_OX_C3 = exp(iVT * (EmQp_SQp - EmbL_bHo)) ## -130mV
        KEQ10_RD_C3 = exp(iVT * (EmQp_SQp - EmbL_bHr)) ## -70mV
        ## Q- + O2 = Q + O2-
        K011_C3 = 100Hz / mM
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
        ## Split of electrical potentials in a reaction (Eyring equation)
        δ₁_C3 = 0.5
        δ₂_C3 = 0.5
        δ₃_C3 = 0.5
        ## Split of the electrical distance across the IMM
        α_C3 = 0.25
        β_C3 = 0.5
        γ_C3 = 0.25
        ## pH factors
        fHi = h_i * inv(1E-7Molar)
        fHm = h_m * inv(1E-7Molar)
        # complex III inhibition by DOX and antimycin
        C3_INHIB = hil(KI_DOX_C3, DOX, 3) * (1 - ANTIMYCIN_BLOCK)
        # Q reduction
        v1 = vQH2C1 + vSDH
        # QH2 diffusion
        v2 = KD_Q * (QH2_n - QH2_p)
        ## Lumped v3 and v4
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        FeS = fes_ox / C3_CONC * (1 - MYXOTHIAZOLE_BLOCK)
        FeSm = fes_rd / C3_CONC * (1 - MYXOTHIAZOLE_BLOCK)
        el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
        er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
        k4ox = K04_C3 * KEQ4_OX_C3 * el4
        k4rd = K04_C3 * KEQ4_RD_C3 * el4
        km4 = K04_C3 * er4 * fHi^2
        v4ox = k4ox * QH2_p * FeS * blo_bho - km4 * Q_p * FeSm * blr_bho
        v4rd = k4rd * QH2_p * FeS * blo_bhr - km4 * Q_p * FeSm * blr_bhr
        ## v5 = Q diffusion (p-side -> n-side)
        v5 = KD_Q * (Q_p - Q_n)
        ## bL- + bH = bL + bH-
        el6 = exp(-iVT * β_C3 * δ₂_C3 * dpsi)
        er6 = exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi)
        k6 = K06_C3 * KEQ6_C3 * el6
        km6 = K06_C3 * er6
        v6 = k6 * blr_bho - km6 * blo_bhr
        ## v7 = bH to Qn; v8: bH to SQn
        ## bH- + Q = bH + Q-
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

        ## cytc1_2+  + cytc_3+ = cytc1_3+  + cytc_2+
        v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

        ## bL- + Qp = bL + Qp-
        k10ox = K010_C3 * KEQ10_OX_C3 * er4
        k10rd = K010_C3 * KEQ10_RD_C3 * er4
        km10 = K010_C3 * el4
        v10ox = k10ox * Q_p * blr_bho - km10 * SQp * blo_bho
        v10rd = k10rd * Q_p * blr_bhr - km10 * SQp * blo_bhr

        ## Qp- + O2 = Qp + O2-
        v11 = K011_C3 * SQp * O2

        eqs = [
            C3_CONC ~ blo_bho + blr_bho + blo_bhr + blr_bhr,
            C3_CONC ~ fes_ox + fes_rd,
            C3_CONC ~ cytc1_ox + cytc1_rd,
            Q_T ~ Q_n + SQn + QH2_n + QH2_p + Q_p + SQp,
            fracbLrd ~ (blr_bho + blr_bhr) / C3_CONC,
            fracbHrd ~ (blo_bhr + blr_bhr) / C3_CONC,
            D(SQp) ~  v10ox + v10rd - v11,
            D(SQn) ~ v7ox + v7rd - v8ox - v8rd,
            D(Q_n) ~ v5 - v7ox - v7rd - v1,
            D(QH2_n) ~ v8ox + v8rd - v2 + v1,
            D(QH2_p) ~ v2 - v4ox - v4rd,
            # D(Q_p) ~ v10 + v4ox + v4rd - v5,
            D(blo_bho) ~ v7ox + v8ox - v4ox + v10ox,
            D(blr_bho) ~ v4ox + v7rd + v8rd - v6 - v10ox,
            D(blo_bhr) ~ v6 - v4rd - v7ox - v8ox + v10rd,
            ## D(blr_bhr) = v4rd - v7rd - v8rd - v10rd
            D(fes_ox) ~ v9 - v4ox - v4rd,
            D(cytc1_ox) ~ v33 - v9,
            D(cytc_ox) ~ 4 * vO2 - v33,
            ## Counting charge movement across the IMM
            vHresC3 ~ α_C3 * (v4ox + v4rd) + β_C3 * v6 + γ_C3 * (v7ox + v7rd + v8ox + v8rd),
            vHres ~ vHresC1 + vHresC3 + vHresC4,
            vROSC3 ~ v11,
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
