using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ECMEDox
using ECMEDox: mM, μM, iVT, mV, Molar, Hz, ms
using Plots

# Adapted from Markevich, 2015
function get_c1_sys(; name=:c1sys,
    Qn=3.0mM, QH2_n=0.3mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.001μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)
    @parameters begin
        Em_O2_SOX = -160mV        # O2/Superoxide redox potential
        Em_FMNsq_FMNH = -375mV    # FMN semiquinone/FMNH- redox potential
        Em_NAD = -320mV           # NAD/NADH redox potential
        Em_FMN_FMNH = -340mV      # FMN/FMNH- redox potential
        Em_N2 = -80mV
        Em_Q_SQ_C1 = -400mV       # -213mV
        Em_SQ_QH2_C1 = +600mV
        ET_C1 = 0.4mM  # Activity of complex I
        KI_NAD_C1 = 1mM
        KI_NADH_C1 = 50μM
        KD_NAD_C1 = 25μM
        KD_NADH_C1 = 100μM
        KI_DOX_C1 = 400μM  # DOX inhibition concentration (IC50) on complex I
        K1_C1 = 10Hz / μM
        KEQ1_C1 = 0.1 / μM     # Association constant for Q
        K2_C1 = 4E5Hz / μM
        KEQ2_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        K3_C1 = 2.7E6Hz / μM
        K4_C1 = 1000Hz
        KEQ4_C1 = 20μM
        K5_C1 = 2Hz / μM
        KEQ5_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        K6_C1 = 0.04Hz / μM
        KEQ6_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        I_C1(t) # Conserved
        Q_C1(t) = 0
        SQ_C1(t) = 0
        QH2_C1(t) = 0
        fFMNH_C1(t)  # fraction of fully-reduced FMN
        fFMNsq_C1(t) # fraction of FMN semiquinone
        N2_C1(t)
        N2m_C1(t)
        KEQ3_C1(t)
        vQ_C1(t)
        vQH2_C1(t)
        vNADH_C1(t)
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
    end

    # Q association
    v1 = K1_C1 * (I_C1 * Qn - Q_C1 / KEQ1_C1)
    # First electron transfer
    v2 = K2_C1 * (Q_C1 * N2m_C1 - SQ_C1 * N2_C1 / KEQ2_C1)
    # Second electron transfer
    fhm = h_m / 1E-7Molar
    v3 = K3_C1 * (SQ_C1 * N2m_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ3_C1)
    # QH2 dissociation
    v4 = K4_C1 * (QH2_C1 - I_C1 * QH2_n / KEQ4_C1)
    # Flavin site ROS generation
    v5 = K5_C1 * ET_C1 * (fFMNH_C1 * O2 - fFMNsq_C1 * sox_m / KEQ5_C1)
    # Quinone site ROS generation
    v6 = K6_C1 * (SQ_C1 * O2 - Q_C1 * sox_m / KEQ6_C1)

    rN2 = exp(iVT * (Em_N2 - Em_NAD)) * (nadh / nad)

    eqs = [
        fFMNH_C1 ~ inv(1 + nad / KD_NAD_C1 + nadh / KI_NADH_C1 + exp(2iVT * (Em_FMN_FMNH - Em_NAD) * (nad / nadh) * (1 + nad / KI_NAD_C1 + nadh / KD_NADH_C1))),
        fFMNsq_C1 ~ exp(iVT * (Em_FMNsq_FMNH - Em_FMN_FMNH)) * fFMNH_C1,
        KEQ3_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_m / h_i)^4,
        N2m_C1 ~ ET_C1 * rN2 / (1 + rN2),
        ET_C1 ~ N2m_C1 + N2_C1,
        # Positive: production; negative: consumption
        vQ_C1 ~ -v1,
        vQH2_C1 ~ v4,
        vROSIf ~ v5,
        vROSIq ~ v6,
        vROS_C1 ~ vROSIf + vROSIq,
        vNADH_C1 ~ -0.5 * (v2 + v3 + v5),
        ET_C1 ~ I_C1 + Q_C1 + SQ_C1 + QH2_C1,
        D(Q_C1) ~ v1 - v2 + v6,
        D(SQ_C1) ~ v2 - v3 - v6,
        D(QH2_C1) ~ v3 - v4,
    ]
    return ODESystem(eqs, t; name)
end

@parameters begin
    Qn = 3.0mM
    QH2_n = 0.3mM
    nad = 500μM
    nadh = 500μM
    dpsi = 150mV
end

sys = get_c1_sys(; Qn, QH2_n, nad, nadh, dpsi) |> structural_simplify

tend = 100.0ms
prob = ODEProblem(sys, [], tend, [sys.nad => 10μM, sys.nadh => 990μM])

sol = solve(prob)

sol[sys.fFMNH_C1][end]
sol[sys.fFMNsq_C1][end]
sol[sys.KEQ3_C1][end]
sol[sys.SQ_C1][end]
sol[sys.vROS_C1][end]
sol[sys.vNADH_C1][end]

@unpack vQ_C1, vQH2_C1, vROS_C1, vNADH_C1 = sys
plot(sol, idxs=[vQ_C1, vQH2_C1, vROS_C1, vNADH_C1], tspan=(20, 40))
