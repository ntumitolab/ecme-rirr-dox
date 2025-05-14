# # Complex III model
# Comparing Gauthier and my complex III models
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using Plots
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms, minute

## complex III and the Q cycle from Gauthier, 2013
## Adapted from Demin, 2001
function c3_gauthier(suc, fum, dpsi;
    MT_PROT=1,
    O2=6μM,
    sox_m=0.001μM,
    h_i=exp10(-7) * Molar,
    h_m=exp10(-7.6) * Molar,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    vQH2C1 = 0,
    cytc_ox = 208μM,
    cytc_rd = 325μM - cytc_ox,
    name = :c3_gauthier
    )

    @parameters begin
        ## Reaction rate constant of SDH (complex II)
        K_C2 = 250 / (minute * mM)
        ## midpoint potential of FUM -> SUC
        Em_FUM_SUC = 30mV
        ## midpoint potential of Q -> QH2
        Em_Q_QH2 = 100mV
        ## (Reverse) equlibrium constant of SDH (-140meV)
        rKEQ_C2 = exp(-2iVT * (Em_Q_QH2 - Em_FUM_SUC))
        rhoC3 = 325μM    # Complex III activity
        rhoQo = rhoC3    # Qo seat
        rhoQi = rhoC3    # Qi seat
        Q_T = 4mM        # Total CoQ pool
        KI_DOX_C3 = 185μM  # DOX inhibition concentration (IC50) on complex III
        EmSQp_QH2p = +290mV
        EmQp_SQp = -170mV
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
        KEQ4_OX_C3 = exp(iVT * (EmbL_bHo - EmQp_SQp)) # +130mV
        KEQ4_RD_C3 = exp(iVT * (EmbL_bHr - EmQp_SQp)) # +70mV
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
        KEQ8_RD_C3 = 9.4546   # +60mV??? should be +190mV?
        K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))  # -35mV
        K010_C3 = 28.33Hz / mM
        KEQ10_C3 = 1.4541 # +10mV
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = 2.1145 # +20mV
    end

    ## complex III inhibition by DOX and antimycin
    C3_CONC = rhoC3 * MT_PROT

    @variables begin
        UQ(t) = Q_T
        UQH2(t) ## Conserved
        Q_n(t)
        QH2_n(t)
        QH2_p(t)
        Q_p(t)
        SQp(t) = 0
        SQn(t) = 0
        fes_ox(t) = C3_CONC
        fes_rd(t) ## Conserved
        cytc1_ox(t) = C3_CONC
        cytc1_rd(t) ## Conserved
        cytb_1(t) = C3_CONC
        cytb_2(t) = 0
        cytb_3(t) = 0
        cytb_4(t) ## Conserved
        vSDH(t)
        vROSC3(t)
        vHresC3(t)
    end

    ## Split of electrical potentials
    δ₁_C3 = 0.5
    δ₂_C3 = 0.5
    δ₃_C3 = 0.5
    ## Split of the electrical distance across the IMM
    α_C3 = 0.25
    β_C3 = 0.5
    γ_C3 = 0.25
    fHi = h_i * inv(1E-7Molar)
    fHm = h_m * inv(1E-7Molar)
    ## Q reduction
    v1 = vQH2C1 + vSDH
    ## QH2 diffusion
    v2 = KD_Q * (QH2_n - QH2_p)
    ## QH2 + FeS = SQp + FeS- + 2H+
    Qo_avail = (rhoQo - SQp) / rhoQo * (1 - MYXOTHIAZOLE_BLOCK)
    v3 = K03_C3 * (KEQ3_C3 * Qo_avail * fes_ox * QH2_p - fes_rd * SQp * fHi^2)
    ## v4: SQp + bL = Qp + bL- (rds?)
    el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
    er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
    v4_ox = K04_C3 * (KEQ4_OX_C3 * SQp * el4 * cytb_1 - Q_p * er4 * cytb_2)
    v4_rd = K04_C3 * (KEQ4_RD_C3 * SQp * el4 * cytb_3 - Q_p * er4 * cytb_4)
    ## v5 = Q diffusion (p-side -> n-side)
    v5 = KD_Q * (Q_p - Q_n)
    ## v6 = bL to bH
    v6 = K06_C3 * (KEQ6_C3 * cytb_2 * exp(-iVT * β_C3 * δ₂_C3 * dpsi) - cytb_3 * exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi))
    ## v7 = bH to Qn; v8: bH to SQn
    Qi_avail = (rhoQi - SQn) / rhoQi * (1 - ANTIMYCIN_BLOCK)
    el7 = exp(-iVT * γ_C3 * δ₃_C3 * dpsi)
    er7 = exp(iVT * γ_C3 * (1 - δ₃_C3) * dpsi)
    qn = Q_n * Qi_avail
    qh2n = QH2_n * Qi_avail
    v7_ox = K07_OX_C3 *  (KEQ7_OX_C3 * cytb_3 * qn * el7 - cytb_1 * SQn * er7)
    v7_rd = K07_RD_C3 * (KEQ7_RD_C3 * cytb_4 * qn * el7 - cytb_2 * SQn * er7)
    v8_ox = K08_OX_C3 * (KEQ8_OX_C3 * cytb_3 * SQn * fHm^2 * el7 - cytb_1 * qh2n * er7)
    v8_rd = K08_RD_C3 * (KEQ8_RD_C3 * cytb_4 * SQn * fHm^2 * el7 - cytb_2 * qh2n * er7)
    ## v9 = fes -> cytc1
    v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
    ## v10: SQp + O2 -> O2- + Q
    v10 = K010_C3 * (KEQ10_C3 * O2 * SQp - sox_m * Q_p)
    ## cytc1_2+  + cytc_3+ = cytc1_3+  + cytc_2+
    v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

    dQn = v5 - v7_ox - v7_rd - v1
    dQp = v10 + v4_ox + v4_rd - v5
    dQH2n = v8_ox + v8_rd - v2 + v1
    dQH2p = v2 - v3
    eqs = [
        vSDH ~ K_C2 * (Q_n * suc - QH2_n * fum * rKEQ_C2),
        C3_CONC ~ cytb_1 + cytb_2 + cytb_3 + cytb_4,
        C3_CONC ~ fes_ox + fes_rd,
        C3_CONC ~ cytc1_ox + cytc1_rd,
        Q_T ~ UQ + SQn + UQH2 + SQp,
        Q_n ~ 0.5 * UQ,
        Q_p ~ 0.5 * UQ,
        QH2_n ~ 0.5 * UQH2,
        QH2_p ~ 0.5 * UQH2,
        D(UQ) ~ dQn + dQp,
        # D(UQH2) ~ dQH2n + dQH2p,
        D(SQn) ~ v7_ox + v7_rd - v8_ox - v8_rd,
        D(SQp) ~ v3 - v10 - v4_ox - v4_rd,
        D(cytb_1) ~ v7_ox + v8_ox - v4_ox,
        D(cytb_2) ~ v4_ox + v7_rd + v8_rd - v6,
        D(cytb_3) ~ v6 - v4_rd - v7_ox - v8_ox,
        # D(cytb_4) = v4_rd - v7_rd - v8_rd
        D(fes_ox) ~ v9 - v3,
        D(cytc1_ox) ~ v33 - v9,
        vHresC3 ~ v3,
        vROSC3 ~ v10,
    ]
    return ODESystem(eqs, t; name)
end

## Semireverse bc1 complex model adapted from Gauthier, 2013
function c3_semireverse(suc, fum, dpsi;
    MT_PROT=1,
    O2=6μM,
    sox_m=0.001μM,
    h_i=exp10(-7) * Molar,
    h_m=exp10(-7.6) * Molar,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    vQH2C1 = 0,
    cytc_ox = 208μM,
    cytc_rd = 325μM - cytc_ox,
    name = :c3_semireverse)

    @parameters begin
        ## Reaction rate constant of SDH (complex II)
        K_C2 = 250 / (minute * mM)
        ## midpoint potential of FUM -> SUC
        Em_FUM_SUC = 30mV
        ## midpoint potential of Q -> QH2
        Em_Q_QH2 = 100mV
        ## (Reverse) equlibrium constant of SDH (-140meV)
        rKEQ_C2 = exp(-2iVT * (Em_Q_QH2 - Em_FUM_SUC))
        rhoC3 = 325μM    # Complex III activity
        Q_T = 4mM        # Total CoQ pool
        KI_DOX_C3 = 185μM  # DOX inhibition concentration (IC50) on complex III
        EmQ_C3 = +100mV  # Ubiquinone redox potential
        EmSQp_QH2p = +290mV
        EmQp_SQp = -170mV
        EmQn_SQn = +70mV
        EmSQn_QH2n = +170mV
        EmbL_bHo = -30mV
        EmbL_bHr = EmbL_bHo - 60mV
        EmbH_bLo = +40mV
        EmbH_bLr = EmbH_bLo - 60mV
        EmFeS = +280mV
        Emcytc1 = +240mV
        EmO2 = -160mV
        Emcytc = +260mV
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmFeS + EmbL_bHo - 2EmQ_C3))
        KEQ4_RD_C3 = exp(iVT * (EmFeS + EmbL_bHr - 2EmQ_C3))
        ## bL- + bH = bL + bH-
        K06_C3 = 166.67Hz
        KEQ6_C3 = exp(iVT * (EmbH_bLo - EmbL_bHo)) # +80mV
        ## bH- + Q = bH + Q-
        K07_OX_C3 = 13.33Hz / mM
        K07_RD_C3 = 1.67Hz / mM
        KEQ7_OX_C3 = exp(iVT * (EmQn_SQn - EmbH_bLo)) # +30mV
        KEQ7_RD_C3 = exp(iVT * (EmQn_SQn - EmbH_bLr)) # +90mV
        ## bH- + Q- + 2H+ = bH + QH2
        K08_OX_C3 = 83.33Hz / mM
        K08_RD_C3 = 8.33Hz / mM
        KEQ8_OX_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLo)) # +130mV
        KEQ8_RD_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLr)) # +60mV??? should be +190mV?
        ## FeS- + c1_3+ = FeS + c1_2+
        K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))  # -40mV
        ## bL- + O2 + Q = bL + O2- + Q
        K010_C3 = 28.33Hz / mM
        KEQ10_OX_C3 = exp(iVT * (EmO2 - EmbL_bHo)) # -130mV
        KEQ10_RD_C3 = exp(iVT * (EmO2 - EmbL_bHr)) # -70mV
        ## c1_2+ + c_3+ = c1_3+ + c_2+
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1)) # +20mV
    end

    ## complex III inhibition by DOX and antimycin
    C3_INHIB = 1 - ANTIMYCIN_BLOCK
    C3_CONC = rhoC3 * MT_PROT

    @variables begin
        UQ(t) = Q_T
        UQH2(t) ## Conserved
        Q_n(t)
        QH2_n(t)
        QH2_p(t)
        Q_p(t)
        SQn(t) = 0
        fes_ox(t) = C3_CONC
        fes_rd(t) ## Conserved
        cytc1_ox(t) = C3_CONC
        cytc1_rd(t) ## Conserved
        blo_bho(t) = C3_CONC
        blr_bho(t) = 0
        blo_bhr(t) = 0
        blr_bhr(t) ## Conserved
        vSDH(t)
        vROSC3(t)
        vHresC3(t)
    end

    ## Split of electrical potentials
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

    ## Q reduction
    v1 = vQH2C1 + vSDH

    ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
    qh2p = QH2_p * (1 - MYXOTHIAZOLE_BLOCK)
    qp = Q_p * (1 - MYXOTHIAZOLE_BLOCK)
    FeS = fes_ox / (fes_ox + fes_rd)
    FeSm = fes_rd / (fes_ox + fes_rd)
    el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
    er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
    k4ox = K04_C3 * KEQ4_OX_C3 * el4
    k4rd = K04_C3 * KEQ4_RD_C3 * el4
    km4 = K04_C3 * er4
    v4ox = k4ox * qh2p * FeS * blo_bho - km4 * qp * FeSm * blr_bho * fHi^2
    v4rd = k4rd * qh2p * FeS * blo_bhr - km4 * qp * FeSm * blr_bhr * fHi^2

    ## bL- + bH = bL + bH-
    el6 = exp(-iVT * β_C3 * δ₂_C3 * dpsi)
    er6 = exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi)
    k6 = K06_C3 * KEQ6_C3 * el6
    km6 = K06_C3 * er6
    v6 = k6 * blr_bho - km6 * blo_bhr

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
    k8ox = K08_OX_C3 * KEQ8_OX_C3 * el7
    k8rd = K08_RD_C3 * KEQ8_RD_C3 * el7
    km8ox = K08_OX_C3 * er7
    km8rd = K08_RD_C3 * er7
    v7ox = k7ox * blo_bhr * qn - km7ox * blo_bho * SQn
    v7rd = k7rd * blr_bhr * qn - km7rd * blr_bho * SQn
    v8ox = k8ox * blo_bhr * SQn * fHm^2 - km8ox * blo_bho * qh2n
    v8rd = k8rd * blr_bhr * SQn * fHm^2 - km8rd * blr_bho * qh2n

    ## FeS- + c1_3+ = FeS + c1_2+
    v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)

    ## cytc1_2+  + cytc_3+ = cytc1_3+  + cytc_2+
    v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

    ## bL- + O2 + Qp = bL + O2- + Qp
    k10ox = K010_C3 * KEQ10_OX_C3 * er4
    k10rd = K010_C3 * KEQ10_RD_C3 * er4
    km10 = K010_C3 * el4
    fq = 2Q_p / Q_T
    v10ox = fq * (k10ox * O2 * blr_bho - km10 * sox_m * blo_bho)
    v10rd = fq * (k10rd * O2 * blr_bhr - km10 * sox_m * blo_bhr)

    dQn = -v7ox - v7rd - v1
    dQp = v10ox + v10rd + v4ox + v4rd
    dQH2n = v8ox + v8rd + v1
    dQH2p = -v4ox - v4rd

    eqs = [
        vSDH ~ K_C2 * (Q_n * suc - QH2_n * fum * rKEQ_C2),
        C3_CONC ~ blo_bho + blr_bho + blo_bhr + blr_bhr,
        C3_CONC ~ fes_ox + fes_rd,
        C3_CONC ~ cytc1_ox + cytc1_rd,
        Q_T ~ UQ + SQn + UQH2,
        Q_n ~ 0.5 * UQ,
        Q_p ~ 0.5 * UQ,
        QH2_n ~ 0.5 * UQH2,
        QH2_p ~ 0.5 * UQH2,
        D(UQ) ~ dQn + dQp,
        # D(UQH2) ~ dQH2n + dQH2p,
        D(SQn) ~ v7ox + v7rd - v8ox - v8rd,
        D(blo_bho) ~ v7ox + v8ox - v4ox,
        D(blr_bho) ~ v4ox + v7rd + v8rd - v6,
        D(blo_bhr) ~ v6 - v4rd - v7ox - v8ox,
        # D(blr_bhr) = v4rd - v7rd - v8rd
        D(fes_ox) ~ v9 - v4ox - v4rd,
        D(cytc1_ox) ~ v33 - v9,
        vHresC3 ~ v4ox + v4rd,
        vROSC3 ~ v10ox + v10rd,
    ]
    return ODESystem(eqs, t; name)
end

#---
@parameters begin
    suc = 20μM
    fum = 10μM
    dpsi = 150mV
    cytc_ox = 208μM
    cytc_rd = 325μM - cytc_ox
end

#---
gsys = c3_gauthier(suc, fum, dpsi; cytc_ox, cytc_rd) |> structural_simplify
rsys = c3_semireverse(suc, fum, dpsi; cytc_ox, cytc_rd) |> structural_simplify
#---
prob_g = SteadyStateProblem(gsys, [])
prob_r = SteadyStateProblem(rsys, [rsys.K010_C3 => 340Hz / mM])

alg = DynamicSS(Rodas5P())
ealg = EnsembleThreads()
extract(sim, k) = map(s -> s[k], sim)

# ## Varying MMP
dpsirange = 100mV:1mV:200mV
alter_dpsi = (prob, i, repeat) -> begin
    prob.ps[dpsi] = dpsirange[i]
    prob
end

eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi)
eprob_r = EnsembleProblem(prob_r; prob_func=alter_dpsi)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)
@time sim_r = solve(eprob_r, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

xs = dpsirange
ys = [extract(sim_g, gsys.vHresC3) extract(sim_r, rsys.vHresC3)]
plot(xs, ys, xlabel="MMP (mV)", ylabel="Resp. Rate (mM/s)", label=["G" "R"])

#---
ys = [extract(sim_g, gsys.vROSC3) extract(sim_r, rsys.vROSC3)]
plot(xs, ys, xlabel="MMP (mV)", ylabel="ROS Rate (mM/s)", label=["G" "R"])

# ## Varying succinate
sucrange = 10μM:10μM:1000μM
alter_suc = (prob, i, repeat) -> begin
    prob.ps[suc] = sucrange[i]
    prob
end

eprob_g = EnsembleProblem(prob_g; prob_func=alter_suc)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(sucrange), abstol=1e-8, reltol=1e-8)

xs = sucrange
ys = [extract(sim_g, gsys.vHresC3) extract(sim_g, gsys.vROSC3) extract(sim_g, gsys.vSDH)]
plot(xs, ys[:, 2], xlabel="Succinate (μM)", ylabel="Rate (mM/s)", label=["G (ROS)"])

plot(xs, extract(sim_g, gsys.QH2_n * 2), xlabel="Succinate (μM)", ylabel="QH2 (μM)", label=["G"])

plot(xs, ys, xlabel="Succinate (μM)", ylabel="Rate (mM/s)", label=["G (Resp.)" "G (ROS)" "G (SDH)"])
