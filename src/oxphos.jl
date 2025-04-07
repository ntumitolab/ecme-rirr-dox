#=
Detailed Oxidative phosphorylation model by Gauthier et al. (2013)
With ROS generation
Default parameter values from Kembro et al. and Gauthier et al.
Some are adjusted by An-Chi Wei to prevent negative concentrations
=#
"Electron transport chain (ETC)"
function get_etc_sys(;
    DOX=0,       # Doxorubcin
    MT_PROT=1,   # OXPHOS protein content scale
    O2=6μM,      # Oxygen
    h_i=exp10(-7) * Molar,      # IMS proton conc
    h_m=exp10(-7.6) * Molar,    # Matrix proton conc
    name=:etcsys,
    ROTENONE_BLOCK = 0,
    ANTIMYCIN_BLOCK = 0,
    MYXOTHIAZOLE_BLOCK = 0,
    CYANIDE_BLOCK = 0,
    OLIGOMYCIN_BLOCK = 0
)
    @parameters begin
        KI_DOX_C3 = 185μM  # DOX inhibition concentration (IC50) on complex III
        KI_DOX_C4 = 165μM  # DOX inhibition concentration (IC50) on complex IV
        E_O2_SOX = -160mV  # O2/Superoxide redox potential
        E_FMN = -375mV     # FMN/FMN- redox potential
    end

    @variables begin
        nad_m(t)
        nadh_m(t)
        dpsi(t)
        sox_m(t)
        suc(t)
        fum(t)
        oaa(t)
        Q_n(t)
        QH2_n(t)
        QH2_p(t)
        Q_p(t) # Conserved
        Qdot_p(t)
        Qdot_n(t)
    end

    # Complex I
    @parameters begin
        KI_DOX_C1 = 400μM  # DOX inhibition concentration (IC50) on complex I
        K_RC_DOX = 1000 / 15mM  # DOX redox cycling constant
        ρC1 = 5mM # Adjusted # 8.85mM, Concentration of complex I, from Gauthier et al. (2013)
        dpsi_B_C1 = 50mV   # Phase boundary potential
        # Transition rates
        K12_C1 = 6.3396e11Hz/mM^2
        K21_C1 = 5Hz
        K56_C1 = 100Hz
        K65_C1 = 2.5119e13Hz/mM^2
        K61_C1 = 1e7Hz
        K16_C1 = 130Hz
        K23_C1 = 3886.7Hz / sqrt(mM)
        K32_C1 = 9.1295e6Hz
        K34_C1 = 639.1364Hz
        K43_C1 = 3.2882Hz / sqrt(mM)
        K47_C1 = 1.5962E7Hz/mM
        K74_C1 = 65.2227Hz
        K75_C1 = 24.615E3Hz
        K57_C1 = 1.1667E3Hz / sqrt(mM)
        K42_C1 = 6.0318Hz / mM
    end

    @variables begin
        vQC1(t)
        vHresC1(t)
        vROSC1(t)
        vNADHC1(t)
        C1_1(t)
        C1_2(t)
        C1_3(t)
        C1_4(t)
        C1_5(t)
        C1_6(t)
        C1_7(t)
        C1_a12(t)
        C1_a21(t)
        C1_a56(t)
        C1_a65(t)
        C1_a61(t)
        C1_a16(t)
        C1_a23(t)
        C1_a32(t)
        C1_a34(t)
        C1_a43(t)
        C1_a47(t)
        C1_a74(t)
        C1_a57(t)
        C1_a75(t)
        C1_a42(t)
        C1_a24(t)
        C1_e1(t)
        C1_e2(t)
        C1_e3(t)
        C1_e4(t)
        C1_e5(t)
        C1_e6(t)
        C1_e7(t)
    end

    c1eqs = let
        C1_CONC = ρC1 * MT_PROT
        # complex I inhibition by DOX and rotenone
        C1_INHIB = hil(KI_DOX_C1, DOX, 3) * (1 - ROTENONE_BLOCK)
        # Electron leak scaling factor from complex I
        E_LEAK_C1 = 1 + K_RC_DOX * DOX
        fv = exp(iVT * (dpsi - dpsi_B_C1))

        # State transition rates
        a12 = C1_a12
        a21 = C1_a21
        a65 = C1_a65
        a56 = C1_a56
        a61 = C1_a61
        a16 = C1_a16
        a23 = C1_a23
        a32 = C1_a32
        a34 = C1_a34
        a43 = C1_a43
        a47 = C1_a47
        a74 = C1_a74
        a57 = C1_a57
        a75 = C1_a75
        a42 = C1_a42
        a24 = C1_a24

        # Fraction of each state in Complex I derived from KA pattern
        denom = C1_e1 + C1_e2 + C1_e3 + C1_e4 + C1_e5 + C1_e6 + C1_e7

        c1eqs = [
            C1_e1 ~ a21 * a32 * a42 * a56 * a61 * a74 + a21 * a32 * a42 * a56 * a61 * a75 + a21 * a32 * a42 * a57 * a61 * a74 + a21 * a32 * a42 * a57 * a65 * a74 + a21 * a32 * a43 * a56 * a61 * a74 + a21 * a32 * a43 * a56 * a61 * a75 + a21 * a32 * a43 * a57 * a61 * a74 + a21 * a32 * a43 * a57 * a65 * a74 + a21 * a32 * a47 * a56 * a61 * a75 + a21 * a34 * a42 * a56 * a61 * a74 + a21 * a34 * a42 * a56 * a61 * a75 + a21 * a34 * a42 * a57 * a61 * a74 + a21 * a34 * a42 * a57 * a65 * a74 + a21 * a34 * a47 * a56 * a61 * a75 + a23 * a34 * a47 * a56 * a61 * a75 + a24 * a32 * a47 * a56 * a61 * a75 + a24 * a34 * a47 * a56 * a61 * a75,
            C1_e2 ~ a12 * a32 * a42 * a56 * a61 * a74 + a12 * a32 * a42 * a56 * a61 * a75 + a12 * a32 * a42 * a57 * a61 * a74 + a12 * a32 * a42 * a57 * a65 * a74 + a12 * a32 * a43 * a56 * a61 * a74 + a12 * a32 * a43 * a56 * a61 * a75 + a12 * a32 * a43 * a57 * a61 * a74 + a12 * a32 * a43 * a57 * a65 * a74 + a12 * a32 * a47 * a56 * a61 * a75 + a12 * a34 * a42 * a56 * a61 * a74 + a12 * a34 * a42 * a56 * a61 * a75 + a12 * a34 * a42 * a57 * a61 * a74 + a12 * a34 * a42 * a57 * a65 * a74 + a12 * a34 * a47 * a56 * a61 * a75 + a16 * a32 * a42 * a57 * a65 * a74 + a16 * a32 * a43 * a57 * a65 * a74 + a16 * a34 * a42 * a57 * a65 * a74,
            C1_e3 ~ a12 * a23 * a42 * a56 * a61 * a74 + a12 * a23 * a42 * a56 * a61 * a75 + a12 * a23 * a42 * a57 * a61 * a74 + a12 * a23 * a42 * a57 * a65 * a74 + a12 * a23 * a43 * a56 * a61 * a74 + a12 * a23 * a43 * a56 * a61 * a75 + a12 * a23 * a43 * a57 * a61 * a74 + a12 * a23 * a43 * a57 * a65 * a74 + a12 * a23 * a47 * a56 * a61 * a75 + a12 * a24 * a43 * a56 * a61 * a74 + a12 * a24 * a43 * a56 * a61 * a75 + a12 * a24 * a43 * a57 * a61 * a74 + a12 * a24 * a43 * a57 * a65 * a74 + a16 * a21 * a43 * a57 * a65 * a74 + a16 * a23 * a42 * a57 * a65 * a74 + a16 * a23 * a43 * a57 * a65 * a74 + a16 * a24 * a43 * a57 * a65 * a74,
            C1_e4 ~ a12 * a23 * a34 * a56 * a61 * a74 + a12 * a23 * a34 * a56 * a61 * a75 + a12 * a23 * a34 * a57 * a61 * a74 + a12 * a23 * a34 * a57 * a65 * a74 + a12 * a24 * a32 * a56 * a61 * a74 + a12 * a24 * a32 * a56 * a61 * a75 + a12 * a24 * a32 * a57 * a61 * a74 + a12 * a24 * a32 * a57 * a65 * a74 + a12 * a24 * a34 * a56 * a61 * a74 + a12 * a24 * a34 * a56 * a61 * a75 + a12 * a24 * a34 * a57 * a61 * a74 + a12 * a24 * a34 * a57 * a65 * a74 + a16 * a21 * a32 * a57 * a65 * a74 + a16 * a21 * a34 * a57 * a65 * a74 + a16 * a23 * a34 * a57 * a65 * a74 + a16 * a24 * a32 * a57 * a65 * a74 + a16 * a24 * a34 * a57 * a65 * a74,
            C1_e5 ~ a12 * a23 * a34 * a47 * a61 * a75 + a12 * a23 * a34 * a47 * a65 * a75 + a12 * a24 * a32 * a47 * a61 * a75 + a12 * a24 * a32 * a47 * a65 * a75 + a12 * a24 * a34 * a47 * a61 * a75 + a12 * a24 * a34 * a47 * a65 * a75 + a16 * a21 * a32 * a42 * a65 * a74 + a16 * a21 * a32 * a42 * a65 * a75 + a16 * a21 * a32 * a43 * a65 * a74 + a16 * a21 * a32 * a43 * a65 * a75 + a16 * a21 * a32 * a47 * a65 * a75 + a16 * a21 * a34 * a42 * a65 * a74 + a16 * a21 * a34 * a42 * a65 * a75 + a16 * a21 * a34 * a47 * a65 * a75 + a16 * a23 * a34 * a47 * a65 * a75 + a16 * a24 * a32 * a47 * a65 * a75 + a16 * a24 * a34 * a47 * a65 * a75,
            C1_e6 ~ a12 * a23 * a34 * a47 * a56 * a75 + a12 * a24 * a32 * a47 * a56 * a75 + a12 * a24 * a34 * a47 * a56 * a75 + a16 * a21 * a32 * a42 * a56 * a74 + a16 * a21 * a32 * a42 * a56 * a75 + a16 * a21 * a32 * a42 * a57 * a74 + a16 * a21 * a32 * a43 * a56 * a74 + a16 * a21 * a32 * a43 * a56 * a75 + a16 * a21 * a32 * a43 * a57 * a74 + a16 * a21 * a32 * a47 * a56 * a75 + a16 * a21 * a34 * a42 * a56 * a74 + a16 * a21 * a34 * a42 * a56 * a75 + a16 * a21 * a34 * a42 * a57 * a74 + a16 * a21 * a34 * a47 * a56 * a75 + a16 * a23 * a34 * a47 * a56 * a75 + a16 * a24 * a32 * a47 * a56 * a75 + a16 * a24 * a34 * a47 * a56 * a75,
            C1_e7 ~ a12 * a23 * a34 * a47 * a56 * a61 + a12 * a23 * a34 * a47 * a57 * a61 + a12 * a23 * a34 * a47 * a57 * a65 + a12 * a24 * a32 * a47 * a56 * a61 + a12 * a24 * a32 * a47 * a57 * a61 + a12 * a24 * a32 * a47 * a57 * a65 + a12 * a24 * a34 * a47 * a56 * a61 + a12 * a24 * a34 * a47 * a57 * a61 + a12 * a24 * a34 * a47 * a57 * a65 + a16 * a21 * a32 * a42 * a57 * a65 + a16 * a21 * a32 * a43 * a57 * a65 + a16 * a21 * a32 * a47 * a57 * a65 + a16 * a21 * a34 * a42 * a57 * a65 + a16 * a21 * a34 * a47 * a57 * a65 + a16 * a23 * a34 * a47 * a57 * a65 + a16 * a24 * a32 * a47 * a57 * a65 + a16 * a24 * a34 * a47 * a57 * a65,
            C1_a12 ~ K12_C1 * h_m^2,
            C1_a21 ~ K21_C1,
            C1_a65 ~ K65_C1 * h_i^2,
            C1_a56 ~ K56_C1,
            C1_a61 ~ K61_C1 / fv,
            C1_a16 ~ K16_C1 * fv,
            C1_a23 ~ K23_C1 * NaNMath.sqrt(nadh_m),
            C1_a32 ~ K32_C1,
            C1_a34 ~ K34_C1,
            C1_a43 ~ K43_C1 * NaNMath.sqrt(nad_m),
            C1_a47 ~ C1_INHIB * K47_C1 * NaNMath.sqrt(Q_n * h_m),
            C1_a74 ~ K74_C1,
            C1_a57 ~ C1_INHIB * K57_C1 * NaNMath.sqrt(QH2_n),
            C1_a75 ~ K75_C1,
            C1_a42 ~ K42_C1 * E_LEAK_C1 * O2,
            C1_a24 ~ K42_C1 * exp(iVT * (E_FMN - E_O2_SOX)) * sox_m,
            vQC1 ~ 0.5 * C1_CONC * (C1_4 * a47 - C1_7 * a74),
            vHresC1 ~ 4 * vQC1,
            vROSC1 ~ C1_CONC * (C1_4 * a42 - C1_2 * a24),
            vNADHC1 ~ 0.5 * C1_CONC * (C1_2 * a23 - C1_3 * a32),
            C1_1 ~ C1_e1 / denom,
            C1_2 ~ C1_e2 / denom,
            C1_3 ~ C1_e3 / denom,
            C1_4 ~ C1_e4 / denom,
            C1_5 ~ C1_e5 / denom,
            C1_6 ~ C1_e6 / denom,
            C1_7 ~ C1_e7 / denom,
        ]
    end

    # Reversible complex II (SDH)
    @parameters begin
        KI_DOX_C2 = 2000μM # DOX inhibition concentration (IC50) on complex II
        K_C2 = 250 / (minute * mM)   # Reaction rate constant of SDH (complex II)
        KI_OAA_C2 = 150μM       # Inhibition constant for OAA
        EmFUM = 40mV            # midpoint potential of FUM -> SUC
        EmQ = 100mV             # midpoint potential of Q -> QH2
        KEQ_C2 = exp(-2iVT * (EmQ - EmFUM)) # (Reverse) equlibrium constant of SDH
    end

    @variables vSDH(t)
    c2eqs = [vSDH ~ K_C2 * (Q_n * suc - QH2_n * fum / KEQ_C2) * hil(KI_OAA_C2, oaa) * hil(KI_DOX_C2, DOX, 3)]

    # complex IV (CCO)
    @parameters begin
        ρC4 = 325μM
        δ₅ = 0.5
        K34_C4 = 1.7667E28 / minute / mM^7
        K43_C4 = 1.7402 / minute / mM^4
        K35_C4 = 45000 / minute / mM
        K36_C4 = 2.8955E25 / minute / mM^4
        K63_C4 = 2.8955E10 / minute / mM^3
        K37_C4 = 1.7542E12 / minute / mM
        K73_C4 = 1.7542E4 / minute / mM
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
        C4_e1(t)
        C4_e2(t)
        C4_e3(t)
        C4_e4(t)
    end

    c4eqs = let
        C4_INHIB = hil(KI_DOX_C4, DOX, 3) * (1 - CYANIDE_BLOCK)  # complex IV inhibition by DOX
        C4_CONC = ρC4 * MT_PROT
        aδ = exp(-iVT * δ₅ * dpsi)
        a1mδ = exp(iVT * (1 - δ₅) * dpsi)
        f_hm = h_m * aδ
        f_hi = h_i * a1mδ
        f_cr = cytc_rd
        f_co = cytc_ox * a1mδ
        a12 = K34_C4 * f_cr^3 * f_hm^4  # a12 = K34 * exp(-δ₅ * 4 * vfrt) * cytc_rd^3 * h_m^4
        a14 = K73_C4 * f_hi  # K73 * exp((1 - δ₅) * F_RT * dpsi) * h_i
        a21 = K43_C4 * f_co^3 * f_hi # K43 * exp((1 - δ₅) * 4 * vfrt) * cytc_ox^3 * h_i
        a23 = C4_INHIB * K35_C4 * O2
        a34 = K36_C4 * f_cr * f_hm^3 # K36 * exp(-δ₅ * 3 * vfrt) * cytc_rd * h_m^3
        a41 = K37_C4 * f_hm  # K37 * exp(-δ₅ * vfrt) * h_m
        a43 = K63_C4 * f_co * f_hi^2  # K63 * exp((1 - δ₅) * 3 * vfrt) * cytc_ox * h_i^2

        # Weight of each state (from KA pattern)
        den = C4_e1 + C4_e2 + C4_e3 + C4_e4
        # Reaction rates
        v34 = C4_CONC * (C4_Y * a12 - C4_Yr * a21)
        v35 = C4_CONC * C4_Yr * a23
        v36 = C4_CONC * (a34 * C4_YO - a43 * C4_YOH)
        v37 = C4_CONC * (a41 * C4_YOH - a14 * C4_Y)

        c4eqs = [
            C4_e1 ~ a21 * a41 * a34 + a41 * a34 * a23,
            C4_e2 ~ a12 * a41 * a34,
            C4_e3 ~ a23 * a12 * a41 + a43 * a14 * a21 + a23 * a43 * a12 + a23 * a43 * a14,
            C4_e4 ~ a14 * a34 * a21 + a34 * a23 * a12 + a34 * a23 * a14,
            C4_Y ~ C4_e1 / den,
            C4_Yr ~ C4_e2 / den,
            C4_YO ~ C4_e3 / den,
            C4_YOH ~ C4_e4 / den,
            vO2 ~ v35,
            vHresC4 ~ v34 + 2 * v36 + v37,
            vCytcOx ~ 3 * v34 + v35,
            C4_CONC ~ cytc_rd + cytc_ox
        ]
    end

    # complex III and the Q cycle
    @parameters begin
        ρC3 = 325μM
        Q_T = 4.0mM  # Total CoQ pool
        E_QH_QH2p = +290mV
        E_Qp_Qdotp = -160mV
        E_Qn_Qdotn = +70mV
        E_Qdotn_QH2n = +170mV
        E_bL = -40mV
        E_bH = +40mV
        E_FeS = +280mV
        E_cytc1 = +245mV
        E_cytc = +255mV
        K03_C3 = 99998 / minute / mM
        KEQ3_C3 = 0.6877
        K04_C3 = 3640.2 / minute / mM
        KEQ4_OX_C3 = 129.9853
        KEQ4_RD_C3 = 13.7484
        δ₁_C3 = 0.5
        δ₂_C3 = 0.5
        δ₃_C3 = 0.5
        β_C3 = 0.5006
        α_C3 = 0.5 * (1 - β_C3)
        γ_C3 = 1 - α_C3 - β_C3
        KD_Q = 1.32E6 / minute
        K06_C3 = 10000 / minute
        KEQ6_C3 = 9.4546
        K07_OX_C3 = 800 / minute / mM
        K07_RD_C3 = 100 / minute / mM
        KEQ7_OX_C3 = 3.0748
        KEQ7_RD_C3 = 29.0714
        K08_OX_C3 = 5000 / minute / mM
        K08_RD_C3 = 500 / minute / mM
        KEQ8_OX_C3 = 129.9853
        KEQ8_RD_C3 = 9.4546
        K09_C3 = 49949 / minute / mM
        KEQ9_C3 = 0.2697
        K010_C3 = 1700 / minute / mM
        KEQ10_C3 = 1.4541
        K33_C3 = 148148 / minute / mM
        KEQ33_C3 = 2.1145
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
        fHi = h_i / 1E-7Molar
        C3_INHIB = hil(KI_DOX_C3, DOX, 3) * (1 - ANTIMYCIN_BLOCK)  # complex III inhibition by DOX and antimycin
        C3_CONC = ρC3 * MT_PROT
        v1 = vQC1 + vSDH # Q reduction
        v2 = KD_Q * (QH2_n - QH2_p) # QH2 diffusion
        # v3 = QH2 to FeS
        Qo_avail = (C3_CONC - Qdot_p) / C3_CONC
        v3 = K03_C3 * (1 - MYXOTHIAZOLE_BLOCK) * (KEQ3_C3 * Qo_avail * fes_ox * QH2_p - fes_rd * Qdot_p * fHi^2)
        # v4 = Qdot_p and bH
        el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
        er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
        v4_ox = K04_C3 * (KEQ4_OX_C3 * Qdot_p * el4 * cytb_1 - Q_p * er4 * cytb_2)
        v4_rd = K04_C3 * (KEQ4_RD_C3 * Qdot_p * el4 * cytb_3 - Q_p * er4 * cytb_4)
        # v5 = Q diffusion (p-side -> n-side)
        v5 = KD_Q * (Q_p - Q_n)
        # v6 = bH to bL
        v6 = K06_C3 * (KEQ6_C3 * cytb_2 * exp(-iVT * β_C3 * δ₂_C3 * dpsi) - cytb_3 * exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi))
        # v7 = bL to Qn; v8: bL to Qdot_n
        el7 = exp(-iVT * γ_C3 * δ₃_C3 * dpsi)
        er7 = exp(iVT * γ_C3 * (1 - δ₃_C3) * dpsi)
        v7_ox = K07_OX_C3 * C3_INHIB * (KEQ7_OX_C3 * cytb_3 * Q_n * el7 - cytb_1 * Qdot_n * er7)
        v7_rd = K07_RD_C3 * C3_INHIB * (KEQ7_RD_C3 * cytb_4 * Q_n * el7 - cytb_2 * Qdot_n * er7)
        FAC_PH = h_m / 1E-7Molar
        v8_ox = K08_OX_C3 * C3_INHIB * (KEQ8_OX_C3 * cytb_3 * Qdot_n * FAC_PH^2 * el7 - cytb_1 * QH2_n * er7)
        v8_rd = K08_RD_C3 * C3_INHIB * (KEQ8_RD_C3 * cytb_4 * Qdot_n * FAC_PH^2 * el7 - cytb_2 * QH2_n * er7)
        # v9 = fes -> cytc1
        v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
        # v10 = Qdot + O2 -> O2- + Q  (ROS produced by complex III)
        v10 = K010_C3 * (KEQ10_C3 * O2 * Qdot_p - sox_m * Q_p)
        # Redox reaction between cytc1 and cytc
        v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)
        [
            D(Q_n) ~ v5 - v7_ox - v7_rd - v1,
            D(Qdot_n) ~ v7_ox + v7_rd - v8_ox - v8_rd,
            D(QH2_n) ~ v8_ox + v8_rd - v2 + v1,
            D(QH2_p) ~ v2 - v3,
            D(Qdot_p) ~ v3 - v10 - v4_ox - v4_rd,
            # D(Q_p) ~ v10 + v4_ox + v4_rd - v5,
            Q_T ~ Q_n + Qdot_n + QH2_n + QH2_p + Qdot_p + Q_p,
            D(cytb_1) ~ v7_ox + v8_ox - v4_ox,
            D(cytb_2) ~ v4_ox + v7_rd + v8_rd - v6,
            D(cytb_3) ~ v6 - v4_rd - v7_ox - v8_ox,
            # D(cytb_4) = v4_rd - v7_rd - v8_rd
            C3_CONC ~ cytb_1 + cytb_2 + cytb_3 + cytb_4,
            D(fes_ox) ~ v9 - v3,
            C3_CONC ~ fes_ox + fes_rd,
            D(cytc1_ox) ~ v33 - v9,
            C3_CONC ~ cytc1_ox + cytc1_rd,
            D(cytc_ox) ~ vCytcOx - v33,
            vHresC3 ~ 2 * v3,
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
        KEQ_C5 = 2.2E5Molar
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
