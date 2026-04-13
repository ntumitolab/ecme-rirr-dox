_inhibit_c3(DOX=0μM, KI_DOX_C3=185μM) = hil(KI_DOX_C3, DOX, 3)

## Semireverse bc1 complex model adapted from Gauthier, 2013
function get_eqs_c3(;
    vQH2C1, vSDH,
    cytc_ox, cytc_rd,
    dpsi, sox_m,
    C3_INHIB_DOX=_inhibit_c3(),     ## Doxorubicin concentration
    MT_PROT=1,                      ## OXPHOS protein content scale factor
    O2=6μM,                         ## Oxygen concentration
    h_i=exp10(-7) * Molar,          ## IMS proton concentration
    h_m=exp10(-7.6) * Molar,        ## Matrix proton concentration
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOL_BLOCK=0,
    STIGMATELLIN_BLOCK=0,
)

    @parameters begin
        rhoC3 = 325μM    ## Complex III activity
        EmQ_C3 = +60mV   ## Ubiquinone redox potential at complex III Qo
        EmSQp_QH2p = +390mV
        EmQp_SQp = -270mV
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
        ## Split of electrical potentials
        δ₁_C3 = 0.5
        δ₂_C3 = 0.5
        δ₃_C3 = 0.5
        ## Split of the electrical distance across the IMM
        α_C3 = 0.25
        β_C3 = 0.5
        γ_C3 = 0.25
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmFeS + EmbL_bHo - 2EmQ_C3))
        KEQ4_RD_C3 = exp(iVT * (EmFeS + EmbL_bHr - 2EmQ_C3))
        ## QH2 and Q flip rate in the IMM
        KD_Q = 22000Hz
        ## bL- + bH = bL + bH-
        K06_C3 = 10000Hz ## 166.67Hz
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
        KEQ11_C3 = exp(iVT * (EmO2 - EmQp_SQp))
        ## c1_2+ + c_3+ = c1_3+ + c_2+
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1)) ## +20mV
    end

    C3_INHIB = _inhibit_c3(C3_INHIB_DOX) * (1 - ANTIMYCIN_BLOCK)
    C3_CONC = rhoC3 * MT_PROT

    sts = @variables begin
        SQn(t) = 142μM
        SQp(t) = 0
        fes_rd(t) = 0
        cytc1_rd(t) = 0
        blr_bhr(t) = 0
        blr_bho(t) = 0
        blo_bhr(t) = 0
    end

    @variables begin
        fes_ox(t) ## Conserved
        cytc1_ox(t) ## Conserved
        bl0_bh0(t) ## Conserved
        fracbLrd(t)
        fracbHrd(t)
        vQpC3(t)
        vQH2pC3(t)
        vQnC3(t)
        vQH2nC3(t)
        vROSC3(t)
        vHresC3(t)
        vCytcC3(t)
    end

    ## pH factors
    fHi = h_i * inv(1E-7Molar)
    fHm = h_m * inv(1E-7Molar)

    ## v1: Q(n) = QH2(n)
    v1 = vQH2C1 + vSDH
    ## v2: QH2(n) = QH2(p)
    v2 = KD_Q * (QH2_n - QH2_p)
    ## v5: Q(p) = Q(n)
    v5 = KD_Q * (Q_p - Q_n)
    ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
    ## Lumped v3 and v4
    FeS = fes_ox / C3_CONC * (1 - MYXOTHIAZOL_BLOCK)
    FeSm = fes_rd / C3_CONC * (1 - MYXOTHIAZOL_BLOCK)
    el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
    er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
    k4ox = K04_C3 * KEQ4_OX_C3 * el4
    k4rd = K04_C3 * KEQ4_RD_C3 * el4
    km4 = K04_C3 * er4 * fHi^2
    v4ox = k4ox * QH2_p * FeS * blo_bho - km4 * Q_p * FeSm * blr_bho
    v4rd = k4rd * QH2_p * FeS * blo_bhr - km4 * Q_p * FeSm * blr_bhr
    ## v6: bL- + bH = bL + bH-
    el6 = exp(-iVT * β_C3 * δ₂_C3 * dpsi)
    er6 = exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi)
    v6 = K06_C3 * (KEQ6_C3 * el6 * blr_bho - blo_bhr * er6)
    ## v7: bH- + Q = bH + Q-
    ## v8: bH- + Q- = bH + QH2
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

    ## v9: FeS- + c1_3+ = FeS + c1_2+
    v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)

    ## v33: cytc1_2+  + cytc_3+ = cytc1_3+  + cytc_2+
    v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

    ## v10: bL- + Qp = bL + Qp-
    k10ox = K010_C3 * KEQ10_OX_C3 * er4
    k10rd = K010_C3 * KEQ10_RD_C3 * er4
    km10 = K010_C3 * el4
    v10ox = k10ox * Q_p * blr_bho - km10 * SQp * blo_bho
    v10rd = k10rd * Q_p * blr_bhr - km10 * SQp * blo_bhr

    ## v11: Qp- + O2 = Qp + O2-
    v11 = K011_C3 * (KEQ11_C3 * SQp * O2 - Q_p * sox_m)

    ## ODEs
    dSQn = v7ox + v7rd - v8ox - v8rd
    dSQp = v10ox + v10rd - v11
    dQn = -v1 + v5 - v7ox - v7rd
    dQp = v4ox + v4rd - v5 - v10ox - v10rd + v11
    dQH2n = v1 - v2 + v8ox + v8rd
    dQH2p = v2 - v4ox - v4rd

    eqs_c3 = [
        C3_CONC ~ cytb_1 + cytb_2 + cytb_3 + cytb_4,
        C3_CONC ~ fes_ox + fes_rd,
        C3_CONC ~ cytc1_ox + cytc1_rd,
        fracbLrd ~ (blr_bho + blr_bhr) / C3_CONC,
        fracbHrd ~ (blo_bhr + blr_bhr) / C3_CONC,

        # TODO: QSSA for SQp (SQp proportion is very small)
        # vROS = Qp * ((k11 * O2 * k10 * bL-) - (km11 * sox * km10 * bL)) / (km10 * bL + k11 * O2)
        # D(Q_p) ~ dQp
        D(Q_n) ~ dQn,
        D(QH2_n) ~ dQH2n,
        D(QH2_p) ~ dQH2p,
        D(SQn) ~ dSQn,
        D(SQp) ~ dSQp,
        D(blo_bho) ~ v7ox + v8ox - v4ox + v10ox,
        D(blr_bho) ~ v4ox + v7rd + v8rd - v6 - v10ox,
        D(blo_bhr) ~ v6 - v4rd - v7ox - v8ox + v10rd,
        ## D(blr_bhr) = v4rd - v7rd - v8rd - v10rd
        D(fes_ox) ~ v9 - v4ox - v4rd,
        D(cytc1_ox) ~ v33 - v9,
        vHresC3 ~ v6, ## Charge movement across the IMM
        vROSC3 ~ v11,
        vCytcC3 ~ v33,
        vQpC3 ~ dQp,
        vQH2pC3 ~ dQH2p,
        vQnC3 ~ dQn,
        vQH2nC3 ~ dQH2n
    ]

    return (; eqs_c3, vROSC3, vHresC3, vCytcC3, SQn, SQp)
end
