# # Complex I model
# Comparing Gauthier, Markevich, and simplified complex I models
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using NaNMath
using Plots
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms

# From Markevich, 2015
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4426091/
function c1_markevich_full(; name=:c1markevich_full,
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
        Em_N1a = −370mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## 800mV (?) in Markevich, 2015
        ET_C1 = 3μM               ## Activity of complex I
        ## DOX IC50 on complex I
        KI_DOX_C1 = 400μM
        K1_C1 = 83Hz/ μM
        KEQ1_C1 = 0.01 / μM
        K2_C1 = 1.44e12Hz
        K3_C1 = 1e6Hz
        KEQ3_C1 = 25μM
        KEQ2_C1 = exp(iVT * (Em_FMN_FMNH - Em_NAD)) / KEQ1_C1 / KEQ3_C1
        K4_C1 = 1Hz/ μM
        KEQ4_C1 = 0.001 / μM
        K5_C1 = 2Hz/ μM
        KEQ5_C1 = 0.02 / μM
        K6_C1 = 5e8Hz/ μM
        KEQ6_C1 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
        K7_C1 = 1E4Hz/ μM
        KEQ7_C1 = exp(iVT * (Em_N2 - Em_N3))
        K8_C1 = 10Hz / μM
        KEQ8_C1 = 0.1 / μM         ## Association constant for Q
        K9_C1 = 4E5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        K10_C1 = 2e6Hz / μM
        KEQ10_C1 = exp(iVT * (Em_N1a - Em_FMN_FMNsq))
        K11_C1 = 1e9Hz / μM
        KEQ11_C1 = exp(iVT * (Em_N3 - Em_FMN_FMNsq))
        K13_C1 = 2.7e6Hz / μM
        K14_C1 = 1000Hz
        KEQ14_C1 = 20μM            ## Dissociation constant for QH2
        K16_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        K17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
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
    end

    fhm = h_m / 1E-7Molar
    ## NADH + FMN = FMN.NADH
    v1 = K1_C1 * (nadh * FMN - FMN_NADH / KEQ1_C1)
    ## FMN.NADH = FMNH−.NAD+
    v2 = K2_C1 * (FMN_NADH - FMNH_NAD / KEQ2_C1)
    ## FMNH−.NAD+ = FMNH− + NAD+
    v3 = K3_C1 * (FMNH_NAD - FMNH * nad / KEQ3_C1)
    ## FMN + NAD+ = FMN.NAD+
    v4 = K4_C1 * (nad * FMN - FMN_NAD / KEQ4_C1)
    ## FMNH− + NADH = FMNH−.NADH
    v5 = K5_C1 * (FMNH * nadh - FMNH_NADH / KEQ5_C1)
    ## FMNH− + N3 = FMNHsq + N3−
    v6 = K6_C1 * (FMNH * N3_C1 - FMNsq * N3r_C1 / KEQ6_C1)
    ## N3− + N2 = N3 + N2−
    v7 = K7_C1 * (N3r_C1 * N2_C1 - N3_C1 * N2r_C1 / KEQ7_C1)
    ## Q association
    v8 = K8_C1 * (Iq_C1 * Q_n - Q_C1 / KEQ8_C1)
    ## CI.Q + N2− = CIQsq + N2
    v9 = K9_C1 * (Q_C1 * N2r_C1 - SQ_C1 * N2_C1 / KEQ9_C1)
    ## FMNHsq + N1a = FMN + N1a− + Hi+
    v10 = K10_C1 * (FMNsq * N1a_C1 - FMN * N1ar_C1 * fhm / KEQ10_C1)
    ## FMNHsq + N3 = FMN + N3− + Hi+
    v11 = K11_C1 * (FMNsq * N3_C1 - FMN * N3r_C1 * fhm / KEQ11_C1)
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = K13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ13_C1)
    ## QH2 dissociation
    v14 = K14_C1 * (QH2_C1 - Iq_C1 * QH2_n / KEQ14_C1)
    ## Flavin site ROS generation
    v16 = K16_C1 * (FMNH * O2 - FMNsq * sox_m / KEQ16_C1)
    ## Quinone site ROS generation
    v17 = K17_C1 * (SQ_C1 * O2 - Q_C1 * sox_m / KEQ17_C1)

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
        D(FMNH) ~ v3 - v5 - v16 - v6,
        D(FMNsq) ~ v6 + v16 - v10 - v11,
        D(N1ar_C1) ~ v10,
        D(N3r_C1) ~ v6 + v11 - v7 - v12,
        D(N2r_C1) ~ v7 + v12 -v9 - v13,
        D(Q_C1) ~ v8 - v9 + v17,
        D(SQ_C1) ~ v9 - v17 - v13,
        D(QH2_C1) ~ v13 - v14,

        ## Positive: production; negative: consumption
        vNADH_C1 ~ -v1,
        vNAD_C1 ~ v3,
        vQ_C1 ~ -v8,
        vQH2_C1 ~ v14,
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROS_C1 ~ vROSIf + vROSIq,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

# simplified from Markevich, 2015
# Assuming immediate eletron transfer from NADH to iron-sulfur chain
function c1_markevich(; name=:c1sys,
    Q_n=3.0mM, QH2_n=0.3mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)
    @parameters begin
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMNsq_FMNH = -375mV    ## FMN semiquinone/FMNH- redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_FMN_FMNH = -340mV      ##FMN/FMNH- avg redox potential
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## 800mV in Markevich, 2015
        ET_C1 = 3μM               ## Activity of complex I
        KI_NAD_C1 = 1mM
        KI_NADH_C1 = 50μM
        KD_NAD_C1 = 25μM
        KD_NADH_C1 = 100μM
        ## DOX IC50 on complex I
        KI_DOX_C1 = 400μM
        K1_C1 = 10Hz / μM
        KEQ1_C1 = 0.1 / μM        ## Association constant for Q
        K2_C1 = 4E5Hz / μM
        KEQ2_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        K3_C1 = 2.7E6Hz / μM
        K4_C1 = 1000Hz
        KEQ4_C1 = 20μM            ## Dissociation constant for QH2
        K5_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ5_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        K6_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        KEQ6_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        I_C1(t) # Conserved
        Q_C1(t) = 0
        SQ_C1(t) = 0
        QH2_C1(t) = 0
        fFMNH_C1(t)  ## fraction of fully-reduced FMN
        fFMNsq_C1(t) ## fraction of FMN semiquinone
        N2_C1(t)
        N2m_C1(t)
        KEQ3_C1(t)
        vQ_C1(t)
        vQH2_C1(t)
        vNADH_C1(t)
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
        TN_C1(t) ## NADH turnover number
    end

    ## Q association
    v1 = K1_C1 * (I_C1 * Q_n - Q_C1 / KEQ1_C1)
    ## First electron transfer
    v2 = K2_C1 * (Q_C1 * N2m_C1 - SQ_C1 * N2_C1 / KEQ2_C1)
    ## Second electron transfer
    fhm = h_m / 1E-7Molar
    v3 = K3_C1 * (SQ_C1 * N2m_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ3_C1)
    ## QH2 dissociation
    v4 = K4_C1 * (QH2_C1 - I_C1 * QH2_n / KEQ4_C1)
    ## Flavin site ROS generation
    v5 = K5_C1 * ET_C1 * (fFMNH_C1 * O2 - fFMNsq_C1 * sox_m / KEQ5_C1)
    ## Quinone site ROS generation
    v6 = K6_C1 * (SQ_C1 * O2 - Q_C1 * sox_m / KEQ6_C1)
    ## FeS N2 reduction ratio
    rN2 = exp(iVT * (Em_N2 - Em_NAD)) * (nadh / nad)

    eqs = [
        fFMNH_C1 ~ inv(1 + nad / KD_NAD_C1 + nadh / KI_NADH_C1 + exp(2iVT * (Em_FMN_FMNH - Em_NAD) * (nad / nadh) * (1 + nad / KI_NAD_C1 + nadh / KD_NADH_C1))),
        fFMNsq_C1 ~ exp(iVT * (Em_FMNsq_FMNH - Em_FMN_FMNH)) * fFMNH_C1,
        KEQ3_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_m / h_i)^4,
        N2m_C1 ~ ET_C1 * rN2 / (1 + rN2),
        ET_C1 ~ N2m_C1 + N2_C1,
        ## Positive: production; negative: consumption
        vQ_C1 ~ -v1,
        vQH2_C1 ~ v4,
        vROSIf ~ v5,
        vROSIq ~ v6,
        vROS_C1 ~ vROSIf + vROSIq,
        vNADH_C1 ~ -0.5 * (v2 + v3 + v5),
        TN_C1 ~ -vNADH_C1 / ET_C1,
        ET_C1 ~ I_C1 + Q_C1 + SQ_C1 + QH2_C1,
        D(Q_C1) ~ v1 - v2 + v6,
        D(SQ_C1) ~ v2 - v3 - v6,
        D(QH2_C1) ~ v3 - v4,
    ]
    return ODESystem(eqs, t; name)
end

# From Gauthier 2012
function c1_gauthier(; name=:c1gauthier,
    Q_n=3.0mM, QH2_n=0.3mM, nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    C1_INHIB=1)
    @parameters begin
        ET_C1 = 8.85mM      ## Activity of complex I
        dpsi_B_C1 = 50mV    ## Phase boundary potential
        ## Transition rates
        K12_C1 = 6339.6Hz ## pH = 7
        K21_C1 = 5Hz
        K56_C1 = 100Hz
        K65_C1 = 251190Hz ## pH = 7
        K61_C1 = 1e7Hz
        K16_C1 = 130Hz
        K23_C1 = 3886.7Hz / sqrt(mM)
        K32_C1 = 9.1295e6Hz
        K34_C1 = 639.1364Hz
        K43_C1 = 3.2882Hz / sqrt(mM)
        K47_C1 = 1.5962E7Hz / mM
        K74_C1 = 65.2227Hz
        K75_C1 = 24.615E3Hz
        K57_C1 = 1.1667E3Hz / sqrt(mM)
        K42_C1 = 6.0318Hz / mM
        Em_O2_SOX = -160mV         ## O2/Superoxide redox potential
        Em_FMNH2_FMNH = -375mV     ## FMNH/FMNH2 redox potential
    end

    @variables begin
        C1_1(t) = 0
        C1_2(t) ## Conserved
        C1_3(t) = 0
        C1_4(t) = 0
        C1_5(t) = 0
        C1_6(t) = 0
        C1_7(t) = 0
        vQ_C1(t)
        vNADH_C1(t)
        vROS_C1(t)
        TN_C1(t) ## Turnover number
    end

    fhi = h_i / 1E-7Molar
    fhm = h_m / 1E-7Molar
    fv = exp(iVT * (dpsi - dpsi_B_C1))
    ## State transition rates
    a12 = K12_C1 * fhm^2
    a21 = K21_C1
    a65 = K65_C1 * fhi^2
    a56 = K56_C1
    a61 = K61_C1 / fv
    a16 = K16_C1 * fv
    a23 = K23_C1 * NaNMath.sqrt(nadh)
    a32 = K32_C1
    a34 = K34_C1
    a43 = K43_C1 * NaNMath.sqrt(nad)
    a47 = C1_INHIB * K47_C1 * NaNMath.sqrt(Q_n * h_m)
    a74 = K74_C1
    a57 = C1_INHIB * K57_C1 * NaNMath.sqrt(QH2_n)
    a75 = K75_C1
    a42 = K42_C1 * O2
    a24 = K42_C1 * exp(iVT * (Em_FMNH2_FMNH - Em_O2_SOX)) * sox_m

    v12 = a12 * C1_1 - a21 * C1_2
    v23 = a23 * C1_2 - a32 * C1_3
    v34 = a34 * C1_3 - a43 * C1_4
    v47 = a47 * C1_4 - a74 * C1_7
    v75 = a75 * C1_7 - a57 * C1_5
    v56 = a56 * C1_5 - a65 * C1_6
    v61 = a61 * C1_6 - a16 * C1_1
    v42 = a42 * C1_4 - a24 * C1_2

    eqs = [
        ET_C1 ~ C1_1 + C1_2 + C1_3 + C1_4 + C1_5 + C1_6 + C1_7,
        D(C1_1) ~ v61 - v12,
        D(C1_3) ~ v23 - v34,
        D(C1_4) ~ v34 - v47 - v42,
        D(C1_7) ~ v47 - v75,
        D(C1_5) ~ v75 - v56,
        D(C1_6) ~ v56 - v61,
        vQ_C1 ~ -0.5v47,
        vROS_C1 ~ v42,
        vNADH_C1 ~ -v23,
        TN_C1 ~ -vNADH_C1 / ET_C1,
    ]
    return ODESystem(eqs, t; name)
end

@parameters begin
    Q_n = 1.8mM
    QH2_n = 0.2mM
    nad = 500μM
    nadh = 500μM
    dpsi = 150mV
end

#---
markevich = c1_markevich(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
m_full = c1_markevich_full(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
gauthier = c1_gauthier(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify

prob_m = SteadyStateProblem(markevich, [markevich.ET_C1 => 3μM, markevich.K5_C1 => 0.004Hz / μM, markevich.K6_C1 => 0.0001Hz / μM, markevich.KEQ1_C1 => 0.02 / μM, markevich.KEQ4_C1 => 50μM])
prob_mf = SteadyStateProblem(m_full, [m_full.ET_C1 => 17μM, m_full.K16_C1 => 0.001Hz / μM, m_full.K17_C1 => 0.001Hz / μM / 50])
prob_g = SteadyStateProblem(gauthier, [])
alg = DynamicSS(Rodas5P())
ealg = EnsembleThreads()

# ## Varying MMP
dpsirange = 100mV:5mV:200mV
alter_dpsi = (prob, i, repeat) -> remake(prob, p=[dpsi => dpsirange[i]])

eprob_mf = EnsembleProblem(prob_mf; prob_func=alter_dpsi, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_dpsi, safetycopy=false)
eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi, safetycopy=false)
@time sim_mf = solve(eprob_mf, alg, ealg; trajectories=length(dpsirange))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(dpsirange))
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(dpsirange))

# Helper function
extract(sim, k) = map(s -> s[k], sim)
# MMP vs NADH turnover
# markevich model has a steeper dependence
xs = dpsirange
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_mf = extract(sim_mf, m_full.vNADH_C1)

plot(xs, [ys_g ys_m ys_mf], xlabel="MMP (mV)", ylabel="NADH rate (μM/ms)", label=["Gauthier" "Markevich" "Markevich (full)"])

# MMP vs ROS production
xs = dpsirange
ys_g = extract(sim_g, gauthier.vROS_C1) .* 1000
ys_m = extract(sim_m, markevich.vROS_C1) .* 1000
ys_mf = extract(sim_mf, m_full.vROS_C1) .* 1000
plot(xs, [ys_g ys_m ys_mf], xlabel="MMP (mV)", ylabel="ROS production (μM/s)", label=["Gauthier" "Markevich" "Markevich (full)"])

#---
ys_if = extract(sim_mf, m_full.vROSIf)
ys_iq = extract(sim_mf, m_full.vROSIq)
plot(xs, [ys_if ys_iq], title="M model (full)", xlabel="MMP (mV)", ylabel="ROS production", label=["IF" "IQ"])

# Inside the ODE system
ys = stack(extract.(Ref(sim_mf), [m_full.Q_C1, m_full.SQ_C1, m_full.QH2_C1, m_full.Iq_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Fraction", label=["Q_C1" "SQ_C1" "QH2_C1" "Iq_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_mf), [m_full.N1ar_C1, m_full.N3r_C1, m_full.N2r_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Fraction", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
@unpack FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq= m_full
ys = stack(extract.(Ref(sim_mf), [FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Conc (μM)", label=["FMN" "FMN_NADH" "FMNH_NAD" "FMN_NAD" "FMNH_NADH" "FMNH" "FMNsq"], legend=:right, lw=1.5)

# ## Varying NADH
nadhrange = 10μM:10μM:990μM
nadrange = 1000μM .- nadhrange
alter_nadh = (prob, i, repeat) -> remake(prob, p=[nadh => nadhrange[i], nad => nadrange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_nadh, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_nadh, safetycopy=false)
eprob_mf = EnsembleProblem(prob_mf; prob_func=alter_nadh, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(nadhrange))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(nadhrange))
@time sim_mf = solve(eprob_mf, alg, ealg; trajectories=length(nadhrange))

# NADH vs turnover
xs = nadhrange
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_mf = extract(sim_mf, m_full.vNADH_C1)

plot(xs, [ys_g ys_m ys_mf], xlabel="NADH (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich" "Markevich (full)"])

# NADH vs ROS production
xs = nadhrange
ys_g = extract(sim_g, gauthier.vROS_C1)
ys_m = extract(sim_m, markevich.vROS_C1)
ys_mf = extract(sim_mf, m_full.vROS_C1)

plot(xs, [ys_g ys_m ys_mf], xlabel="NADH (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "Markevich (full)"])

#---
ys_if = extract(sim_mf, m_full.vROSIf)
ys_iq = extract(sim_mf, m_full.vROSIq)
plot(xs, [ys_if ys_iq], title="M model (full)", xlabel="NADH (μM)", ylabel="ROS production", label=["IF" "IQ"])

# Inside the ODE system
ys = stack(extract.(Ref(sim_mf), [m_full.Q_C1, m_full.SQ_C1, m_full.QH2_C1, m_full.Iq_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Fraction", label=["Q_C1" "SQ_C1" "QH2_C1" "Iq_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_mf), [m_full.N1ar_C1, m_full.N3r_C1, m_full.N2r_C1]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Fraction", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
@unpack FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq= m_full
ys = stack(extract.(Ref(sim_mf), [FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq]), dims=2)
plot(xs, ys, xlabel="NADH (μM)", ylabel="Conc (μM)", label=["FMN" "FMN_NADH" "FMNH_NAD" "FMN_NAD" "FMNH_NADH" "FMNH" "FMNsq"], legend=:right, lw=1.5)

# ## Varying Q
qh2range = 10μM:10μM:1990μM
qrange = 2000μM .- qh2range
alter_qh2 = (prob, i, repeat) -> remake(prob, p=[QH2_n => qh2range[i], Q_n => qrange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_qh2, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_qh2, safetycopy=false)
eprob_mf = EnsembleProblem(prob_mf; prob_func=alter_qh2, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(qh2range))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(qh2range))
@time sim_mf = solve(eprob_mf, alg, ealg; trajectories=length(qh2range))
# QH2 vs NADH turnover
xs = qh2range
ys_g = extract(sim_g, gauthier.vNADH_C1)
ys_m = extract(sim_m, markevich.vNADH_C1)
ys_mf = extract(sim_mf, m_full.vNADH_C1)

plot(xs, [ys_g ys_m ys_mf], xlabel="QH2 (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich" "Markevich (full)"])

# QH2 vs ROS production
# Gauthier model produces a lot of SOX on high QH2
xs = qh2range
ys_g = extract(sim_g, gauthier.vROS_C1)
ys_m = extract(sim_m, markevich.vROS_C1)
ys_mf = extract(sim_mf, m_full.vROS_C1)
plot(xs, [ys_g ys_m ys_mf], xlabel="QH2 (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "Markevich (full)"])

#---
@unpack C1_1, C1_2, C1_3, C1_4, C1_5, C1_6, C1_7 = gauthier
ys = stack(extract.(Ref(sim_g), [C1_3, C1_4, C1_6]), dims=2)
plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc (μM)", label=["C1_3" "C1_4" "C1_6"], legend=:right, lw=1.5)

#---
ys_if = extract(sim_mf, m_full.vROSIf)
ys_iq = extract(sim_mf, m_full.vROSIq)
plot(xs, [ys_if ys_iq], title="M model (full)", xlabel="MMP (mV)", ylabel="ROS production", label=["IF" "IQ"])

# Inside the ODE system
ys = stack(extract.(Ref(sim_mf), [m_full.Q_C1, m_full.SQ_C1, m_full.QH2_C1, m_full.Iq_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Fraction", label=["Q_C1" "SQ_C1" "QH2_C1" "Iq_C1"], legend=:right)

#---
ys = stack(extract.(Ref(sim_mf), [m_full.N1ar_C1, m_full.N3r_C1, m_full.N2r_C1]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Fraction", label=["N1ar" "N3r" "N2r"], legend=:right)

#---
@unpack FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq= m_full
ys = stack(extract.(Ref(sim_mf), [FMN, FMN_NADH, FMNH_NAD, FMN_NAD, FMNH_NADH, FMNH, FMNsq]), dims=2)
plot(xs, ys, xlabel="MMP (mV)", ylabel="Conc (μM)", label=["FMN" "FMN_NADH" "FMNH_NAD" "FMN_NAD" "FMNH_NADH" "FMNH" "FMNsq"], legend=:right, lw=1.5)
