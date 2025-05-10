#=
Detailed Oxidative phosphorylation model by Gauthier et al. (2013)
With ROS generation
Default parameter values from Kembro et al. and Gauthier et al.
Some are adjusted by An-Chi Wei to prevent negative concentrations
=#

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
        SQp(t)
        SQn(t)
    end

    # Complex I using a simplified Markevich model
    # Rapid equlibrium in the flavin site
    # QSSA for the catalytic cycle in the quinone site
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

    c1eqs = let
        C1_CONC = ET_C1 * MT_PROT
        ## complex I inhibition by DOX and rotenone
        C1_INHIB = hil(KI_DOX_C1, DOX, 3) * (1 - ROTENONE_BLOCK)
        ## Electron leak scaling factor from complex I
        E_LEAK_C1 = 1 + K_RC_DOX * DOX
        fhm = h_m * inv(1E-7Molar)
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
        ## FMN + NADH = FMNH- + NAD+
        rFMNH_FMN = (nadh_m / nad_m) * KEQ2_C1
        ## FMNHsq + N3 = FMN + N3− + Hi+
        rFMNHsq_FMN = rN3_C1 * fhm * rKEQ11_C1
        ## Weights in the flavin site
        denf = 1 + rFMNH_FMN + rFMNHsq_FMN
        fFMN = KI_NAD_C1 / (nad_m + KI_NAD_C1)
        fFMNH = KI_NADH_C1 / (nadh_m + KI_NADH_C1)

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
        qC1 = ET_C1 / qDen

        eqs = [
            D(N2r_C1) ~ v7 + v12 - v9 - v13,
            rKEQ13_C1 ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
            wIq ~ b21 * b32 * b41 + b21 * b32 * b43 + b21 * b34 * b41 + b23 * b34 * b41,
            wIqQ ~ b12 * b32 * b41 + b12 * b32 * b43 + b12 * b34 * b41 + b14 * b32 * b43,
            wIqSQ ~ b12 * b23 * b41 + b12 * b23 * b43 + b14 * b21 * b43 + b14 * b23 * b43,
            wIqQH2 ~ b12 * b23 * b34 + b14 * b21 * b32 + b14 * b21 * b34 + b14 * b23 * b34,
            Iq_C1 ~ wIq * qC1,
            Q_C1 ~ wIqQ * qC1,
            SQ_C1 ~ wIqSQ * qC1,
            QH2_C1 ~ wIqQH2 * qC1,
            ET_C1 ~ N2_C1 + N2r_C1,
            ET_C1 ~ N3_C1 + N3r_C1,
            N3_C1 ~ ET_C1 / (1 + rN3_C1),
            rN3_C1 ~ KEQ_NADH_N3 * NaNMath.sqrt(nadh_m / (nad_m * fhm)),
            FMN ~ ET_C1 / denf * fFMN,
            FMN_NAD ~ ET_C1 / denf * (1 - fFMN),
            FMNsq ~ ET_C1 / denf * rFMNHsq_FMN,
            FMNH ~ ET_C1 / denf * rFMNH_FMN * fFMNH,
            FMNH_NADH ~ ET_C1 / denf * rFMNH_FMN * (1 - fFMNH),
            vQC1 ~ -v8,
            vNADHC1 ~ -0.5 * (v7 + v12 + v16),
            vROSIf ~ v16,
            vROSIq ~ v17,
            vROSC1 ~ vROSIf + vROSIq,
            vQH2C1 ~ v14,
            vHresC1 ~ 4v13,
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
        vCytcOx(t)
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
            vCytcOx ~ 4vO2,
            vHresC4 ~ vCytcOx,
        ]
    end

    # complex III and the Q cycle
    @parameters begin
        rhoC3 = 325μM
        rhoQo = rhoC3    # Qo seat
        rhoQi = rhoC3    # Qi seat
        Q_T = 4mM        # Total CoQ pool
        KI_DOX_C3 = 185μM  # DOX inhibition concentration (IC50) on complex III
        EmSQp_QH2p = +290mV
        EmQp_SQp = -160mV
        EmQn_SQn = +70mV
        EmSQn_QH2n = +170mV
        EmbL_bHo = -40mV
        EmbL_bHr = EmbL_bHo - 60mV
        EmbH_bLo = +40mV
        EmbH_bLr = EmbH_bLo - 60mV
        EmFeS = +280mV
        Emcytc1 = +245mV
        Emcytc = +255mV
        K03_C3 = 1666.63Hz / mM
        KEQ3_C3 = exp(iVT * (EmFeS - EmSQp_QH2p)) # -10mV
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = 129.9853 # +130mV
        KEQ4_RD_C3 = 13.7484  # +70mV
        # Split of electrical potentials
        δ₁_C3 = 0.5
        δ₂_C3 = 0.5
        δ₃_C3 = 0.5
        # Split of the electrical distance across the IMM
        α_C3 = 0.25
        β_C3 = 0.5
        γ_C3 = 0.25
        KD_Q = 22000Hz
        K06_C3 = 166.67Hz
        KEQ6_C3 = 9.4546 # +60mV
        K07_OX_C3 = 13.33Hz / mM
        K07_RD_C3 = 1.67Hz / mM
        KEQ7_OX_C3 = 3.0748 # +30mV
        KEQ7_RD_C3 = 29.0714 # +90mV
        K08_OX_C3 = 83.33Hz / mM
        K08_RD_C3 = 8.33Hz / mM
        KEQ8_OX_C3 = 129.9853 # +130mV
        KEQ8_RD_C3 = 9.4546   # +60mV
        K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))  # -35mV
        K010_C3 = 28.33Hz / mM
        KEQ10_C3 = 1.4541 # +10mV
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = 2.1145 # +20mV
    end

    @variables begin
        fes_ox(t)
        fes_rd(t) # Conserved
        cytc1_ox(t)
        cytc1_rd(t) # Conserved
        cytb_1(t)
        cytb_2(t)
        cytb_3(t)
        cytb_4(t) # Conserved
        vROSC3(t)
        vHres(t)
        vHresC3(t)
        vROS(t)
    end

    c3eqs = let
        fHi = h_i * inv(1E-7Molar)
        fHm = h_m * inv(1E-7Molar)
        # complex III inhibition by DOX and antimycin
        C3_INHIB = hil(KI_DOX_C3, DOX, 3) * (1 - ANTIMYCIN_BLOCK)
        C3_CONC = rhoC3 * MT_PROT
        # Q reduction
        v1 = vQH2C1 + vSDH
        # QH2 diffusion
        v2 = KD_Q * (QH2_n - QH2_p)
        ## QH2 + FeS = SQp + FeS- + 2H+
        Qo_avail = (rhoQo - SQp) / rhoQo * (1 - MYXOTHIAZOLE_BLOCK)
        v3 = K03_C3 * (KEQ3_C3 * Qo_avail * fes_ox * QH2_p - fes_rd * SQp * fHi^2)
        # v4: SQp + bL = Qp + bL-
        el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
        er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
        v4_ox = K04_C3 * (KEQ4_OX_C3 * SQp * el4 * cytb_1 - Q_p * er4 * cytb_2)
        v4_rd = K04_C3 * (KEQ4_RD_C3 * SQp * el4 * cytb_3 - Q_p * er4 * cytb_4)
        # v5 = Q diffusion (p-side -> n-side)
        v5 = KD_Q * (Q_p - Q_n)
        # v6 = bL to bH
        v6 = K06_C3 * (KEQ6_C3 * cytb_2 * exp(-iVT * β_C3 * δ₂_C3 * dpsi) - cytb_3 * exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi))
        # v7 = bH to Qn; v8: bH to SQn
        Qi_avail = (rhoQi - SQn) / rhoQi * C3_INHIB
        el7 = exp(-iVT * γ_C3 * δ₃_C3 * dpsi)
        er7 = exp(iVT * γ_C3 * (1 - δ₃_C3) * dpsi)
        v7_ox = K07_OX_C3 *  (KEQ7_OX_C3 * cytb_3 * Q_n * Qi_avail * el7 - cytb_1 * SQn * er7)
        v7_rd = K07_RD_C3 * (KEQ7_RD_C3 * cytb_4 * Q_n * Qi_avail * el7 - cytb_2 * SQn * er7)
        v8_ox = K08_OX_C3 * (KEQ8_OX_C3 * cytb_3 * SQn * fHm^2 * el7 - cytb_1 * QH2_n * Qi_avail * er7)
        v8_rd = K08_RD_C3 * (KEQ8_RD_C3 * cytb_4 * SQn * fHm^2 * el7 - cytb_2 * QH2_n * Qi_avail * er7)
        # v9 = fes -> cytc1
        v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
        # v10: SQp + O2 -> O2- + Q
        v10 = K010_C3 * (KEQ10_C3 * O2 * SQp - sox_m * Q_p)
        # cytc1_2+  + cytc_3+ = cytc1_3+  + cytc_2+
        v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)
        [
            C3_CONC ~ cytb_1 + cytb_2 + cytb_3 + cytb_4,
            C3_CONC ~ fes_ox + fes_rd,
            C3_CONC ~ cytc1_ox + cytc1_rd,
            Q_T ~ Q_n + SQn + QH2_n + QH2_p + SQp + Q_p,
            D(Q_n) ~ v5 - v7_ox - v7_rd - v1,
            D(SQn) ~ v7_ox + v7_rd - v8_ox - v8_rd,
            D(QH2_n) ~ v8_ox + v8_rd - v2 + v1,
            D(QH2_p) ~ v2 - v3,
            D(SQp) ~ v3 - v10 - v4_ox - v4_rd,
            # D(Q_p) ~ v10 + v4_ox + v4_rd - v5,
            D(cytb_1) ~ v7_ox + v8_ox - v4_ox,
            D(cytb_2) ~ v4_ox + v7_rd + v8_rd - v6,
            D(cytb_3) ~ v6 - v4_rd - v7_ox - v8_ox,
            # D(cytb_4) = v4_rd - v7_rd - v8_rd
            D(fes_ox) ~ v9 - v3,
            D(cytc1_ox) ~ v33 - v9,
            D(cytc_ox) ~ vCytcOx - v33,
            vHresC3 ~ v3,
            vHres ~ vHresC1 + vHresC3 + vHresC4,
            vROSC3 ~ v10,
            vROS ~ vROSC3 + vROSC1,
        ]
    end

    return ODESystem([c1eqs; c2eqs; c4eqs; c3eqs], t; name)
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
    return ODESystem(eqs, t; name)
end
