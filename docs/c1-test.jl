# # Complex I model
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SteadyStateDiffEq
using OrdinaryDiffEq
using NaNMath
using Plots
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms

function c1_2state(; name=:c1sys,
    Q_n=3.0mM, QH2_n=0.3mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,)

    @parameters begin
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMNsq_FMNH = -375mV    ## FMN semiquinone/FMNH- redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -350mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +550mV     ## 800mV in Markevich, 2015
        ET_C1 = 3μM               ## Activity of complex I
        KI_NAD_C1 = 1000μM
        KI_NADH_C1 = 50μM
        KD_NAD_C1 = 25μM
        KD_NADH_C1 = 100μM
        ## First electron transfer
        K1_C1 = 10Hz / μM
        KdQ_C1 = 1000μM         ## Association constant for Q
        ## Second electron transfer
        K2_C1 = 1000Hz
        KdQH2_C1 = 1000μM          ## Dissociation constant for QH2
        ## If site SOX production
        K5_C1 = 2Hz / μM
        KEQ5_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        ## Iq site SOX production
        K6_C1 = 0.04Hz / μM
        KEQ6_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        I_C1(t) ## Conserved
        SQ_C1(t) = 0
        fFMNH_C1(t)  ## fraction of fully-reduced FMN
        fFMNsq_C1(t) ## fraction of FMN semiquinone
        KEQ1_C1(t)
        KEQ2_C1(t)
        vQ_C1(t)
        vQH2_C1(t)
        vNADH_C1(t)
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
        TN_C1(t) ## NADH turnover number
    end

    k1 = K1_C1
    km1 = K1_C1 / KEQ1_C1 * KdQ_C1
    v1 = k1 * Q_n * I_C1 - km1 * SQ_C1

    fhm = h_m / 1E-7Molar
    k2 = K2_C1
    km2 = K2_C1 / KEQ2_C1 / KdQH2_C1
    v2 = k2 * SQ_C1 * fhm^2 - km2 * I_C1 * QH2_n

    ## Flavin site ROS generation
    v5 = K5_C1 * ET_C1 * (fFMNH_C1 * O2 - fFMNsq_C1 * sox_m / KEQ5_C1)
    ## Quinone site ROS generation
    v6 = K6_C1 * (SQ_C1 * O2 - I_C1 * Q_n * sox_m / KEQ6_C1 / KdQ_C1)

    eqs = [
        fFMNH_C1 ~ inv(1 + nad / KD_NAD_C1 + nadh / KI_NADH_C1 + exp(2iVT * (Em_FMN_FMNH - Em_NAD) * (nad / nadh) * (1 + nad / KI_NAD_C1 + nadh / KD_NADH_C1))),
        fFMNsq_C1 ~ exp(iVT * (Em_FMNsq_FMNH - Em_FMN_FMNH)) * fFMNH_C1,
        KEQ1_C1 ~ exp(iVT * (Em_Q_SQ_C1 - Em_NAD)) * (nadh / nad),
        KEQ2_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_NAD - 4dpsi)) * (h_m / h_i)^4 * (nadh / nad),
        vQ_C1 ~ -v1 + v6,
        vQH2_C1 ~ v2,
        vROSIf ~ v5,
        vROSIq ~ v6,
        vROS_C1 ~ vROSIf + vROSIq,
        vNADH_C1 ~ -0.5 * (v1 + v2 + v5),
        TN_C1 ~ -vNADH_C1 / ET_C1,
        ET_C1 ~ I_C1 + SQ_C1,
        D(SQ_C1) ~ v1 - v2 - v6,
    ]
    return ODESystem(eqs, t; name)
end

# Adapted from Markevich, 2015
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

function c1_gauthier(; name=:c1sys,
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
    Q_n = 3.0mM
    QH2_n = 0.3mM
    nad = 500μM
    nadh = 500μM
    dpsi = 150mV
end


#---
markevich = c1_markevich(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify
gauthier = c1_gauthier(; Q_n, QH2_n, nad, nadh, dpsi) |> structural_simplify

prob_m = SteadyStateProblem(markevich, [markevich.ET_C1 => 3μM, markevich.K5_C1 => 0.02Hz / μM, markevich.K6_C1 => 0.0004Hz / μM, markevich.KEQ1_C1 => 0.02 / μM, markevich.KEQ4_C1 => 50μM])
prob_g = SteadyStateProblem(gauthier, [])
alg = DynamicSS(Rodas5P())
ealg = EnsembleThreads()

# ## Varying MMP
dpsirange = 100mV:5mV:200mV
alter_dpsi = (prob, i, repeat) -> remake(prob, p=[dpsi => dpsirange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_dpsi, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(dpsirange))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(dpsirange))

# MMP vs NADH turnover
# markevich model has a steeper dependence
xs = dpsirange
ys_g = map(sim_g) do sol
    sol[gauthier.vNADH_C1]
end
ys_m = map(sim_m) do sol
    sol[markevich.vNADH_C1]
end
plot(xs, [ys_g ys_m], xlabel="MMP (mV)", ylabel="NADH rate (μM/ms)", label=["Gauthier" "Markevich"])

# MMP vs ROS production
xs = dpsirange
ys_g = map(sim_g) do sol
    sol[gauthier.vROS_C1] * 1000
end
ys_m = map(sim_m) do sol
    sol[markevich.vROS_C1] * 1000
end
plot(xs, [ys_g ys_m], xlabel="MMP (mV)", ylabel="ROS production (μM/s)", label=["Gauthier" "Markevich"])

#---
ys_if = map(sim_m) do sol
    sol[markevich.vROSIf]
end
ys_iq = map(sim_m) do sol
    sol[markevich.vROSIq]
end
plot(xs, [ys_if ys_iq], xlabel="MMP (mV)", ylabel="ROS production", label=["IF" "IQ"])

#---
ys_qc1 = map(sim_m) do sol
    sol[markevich.Q_C1]
end
ys_sqc1 = map(sim_m) do sol
    sol[markevich.SQ_C1]
end
ys_qh2c1 = map(sim_m) do sol
    sol[markevich.QH2_C1]
end
ys_c1 = map(sim_m) do sol
    sol[markevich.I_C1]
end

plot(xs, [ys_c1 ys_qc1 ys_sqc1 ys_qh2c1], xlabel="MMP (mV)", ylabel="Fraction", label=["I_C1" "Q_C1" "SQ_C1" "QH2_C1"])

# ## Varying NADH
nadhrange = 10μM:10μM:990μM
alter_nadh = (prob, i, repeat) -> remake(prob, p=[nadh => nadhrange[i], nad=>1000μM - nadhrange[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_nadh, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_nadh, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(nadhrange))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(nadhrange))

# NADH vs turnover
xs = nadhrange
ys_g = map(sim_g) do sol
    sol[gauthier.vNADH_C1]
end
ys_m = map(sim_m) do sol
    sol[markevich.vNADH_C1]
end

plot(xs, [ys_g ys_m], xlabel="NADH (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich"])

# NADH vs ROS production
xs = nadhrange
ys_g = map(sim_g) do sol
    sol[gauthier.vROS_C1]
end
ys_m = map(sim_m) do sol
    sol[markevich.vROS_C1]
end

plot(xs, [ys_g ys_m], xlabel="NADH (μM)", ylabel="ROS production", label=["Gauthier" "Markevich"])

#---
ys_qc1 = map(sim_m) do sol
    sol[markevich.Q_C1]
end
ys_sqc1 = map(sim_m) do sol
    sol[markevich.SQ_C1]
end
ys_qh2c1 = map(sim_m) do sol
    sol[markevich.QH2_C1]
end
ys_c1 = map(sim_m) do sol
    sol[markevich.I_C1]
end

plot(xs, [ys_c1 ys_qc1 ys_sqc1 ys_qh2c1], xlabel="MMP (mV)", ylabel="Fraction", label=["I_C1" "Q_C1" "SQ_C1" "QH2_C1"])

# ## Varying Q
qh2range = 10μM:10μM:1990μM
alter_qh2 = (prob, i, repeat) -> remake(prob, p=[QH2_n => qh2range[i], Q_n=>2000μM - qh2range[i]])

eprob_g = EnsembleProblem(prob_g; prob_func=alter_qh2, safetycopy=false)
eprob_m = EnsembleProblem(prob_m; prob_func=alter_qh2, safetycopy=false)
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(qh2range))
@time sim_m = solve(eprob_m, alg, ealg; trajectories=length(qh2range))

# QH2 vs turnover
xs = qh2range
ys_g = map(sim_g) do sol
    sol[gauthier.vNADH_C1]
end
ys_m = map(sim_m) do sol
    sol[markevich.vNADH_C1]
end

plot(xs, [ys_g ys_m], xlabel="QH2 (μM)", ylabel="NADH consumption (μM/ms)", label=["Gauthier" "Markevich"])

# QH2 vs ROS production
xs = qh2range
ys_g = map(sim_g) do sol
    sol[gauthier.vROS_C1]
end
ys_m = map(sim_m) do sol
    sol[markevich.vROS_C1]
end

plot(xs, [ys_g ys_m], xlabel="QH2 (μM)", ylabel="ROS production", label=["Gauthier" "Markevich"])

#---
ys_if = map(sim_m) do sol
    sol[markevich.vROSIf]
end
ys_iq = map(sim_m) do sol
    sol[markevich.vROSIq]
end
plot(xs, [ys_if ys_iq], xlabel="QH2 (μM)", ylabel="ROS production", label=["IF" "IQ"])

#---
ys_qc1 = map(sim_m) do sol
    sol[markevich.Q_C1]
end
ys_sqc1 = map(sim_m) do sol
    sol[markevich.SQ_C1]
end
ys_qh2c1 = map(sim_m) do sol
    sol[markevich.QH2_C1]
end
ys_c1 = map(sim_m) do sol
    sol[markevich.I_C1]
end

plot(xs, [ys_c1 ys_qc1 ys_sqc1 ys_qh2c1], xlabel="QH2 (μM)", ylabel="Fraction", label=["I_C1" "Q_C1" "SQ_C1" "QH2_C1"])
