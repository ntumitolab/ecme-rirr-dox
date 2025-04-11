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
        Em_Q_SQ_C1 = -400mV
        Em_SQ_QH2_C1 = +600mV
        ET_C1 = 0.4mM  # Activity of complex I
        KI_NAD_C1 = 1mM
        KI_NADH_C1 = 50μM
        KD_NAD_C1 = 25μM
        KD_NADH_C1 = 100μM
        KI_DOX_C1 = 400μM  # DOX inhibition concentration (IC50) on complex I
        KQ_C1 = 0.1 / μM     # Association constant for Q
        KQH2_C1 = 0.05 / μM   # Association constant for QH2
        K1_C1 = 100Hz
        K2_C1 = 1000Hz
        K3_C1 = 0.01Hz / μM
        KEQ3_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        K4_C1 = 0.005Hz / μM
        KEQ4_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
    end

    @variables begin
        C1(t) # Conserved
        SQ_C1(t) = 0
        fFMNH_C1(t)  # fraction of fully-reduced FMN
        fFMNsq_C1(t) # fraction of FMN semiquinone
        KEQ1_C1(t)
        KEQ2_C1(t)
        vQ_C1(t)
        vQH2_C1(t)
        vNADH_C1(t)
        vROSIf(t)
        vROSIq(t)
        vROS_C1(t)
        v1_C1(t)
        v2_C1(t)
        v3_C1(t)
        v4_C1(t)
    end

    # First electron transfer
    v1 = K1_C1 * (KEQ1_C1 * KQ_C1 * (ET_C1 - SQ_C1) * Qn - SQ_C1)
    fhm = h_m / 1E-7Molar
    # Second electron transfer
    v2 = K2_C1 * (KEQ2_C1 * SQ_C1 * fhm^2 - KQH2_C1 * (ET_C1 - SQ_C1) * QH2_n)
    # If site ROS generation
    v3 = K3_C1 * ET_C1 * (KEQ3_C1 * fFMNH_C1 * O2 - fFMNsq_C1 * sox_m)
    # IQ site ROS generation
    v4 = K4_C1 * (KEQ4_C1 * SQ_C1 * O2 - KQ_C1 * (ET_C1 - SQ_C1) * Qn * sox_m)

    eqs = [
        fFMNH_C1 ~ inv(1 + nad / KD_NAD_C1 + nadh / KI_NADH_C1 + exp(2iVT * (Em_FMN_FMNH - Em_NAD) * (nad / nadh) * (1 + nad / KI_NAD_C1 + nadh / KD_NADH_C1))),
        fFMNsq_C1 ~ exp(iVT * (Em_FMNsq_FMNH - Em_FMN_FMNH)) * fFMNH_C1,
        KEQ1_C1 ~ exp(iVT * (Em_Q_SQ_C1 - Em_NAD)) * (nadh / nad),
        KEQ2_C1 ~ exp(iVT * (Em_SQ_QH2_C1 - Em_NAD - 4dpsi)) * (nadh / nad) * (h_m/h_i)^4 ,
        v1_C1 ~ v1,
        v2_C1 ~ v2,
        v3_C1 ~ v3,
        v4_C1 ~ v4,
        # Positive: production; negative: consumption
        vQ_C1 ~ -v1_C1 + v4_C1,
        vQH2_C1 ~ v2_C1,
        vROSIf ~ v3_C1,
        vROSIq ~ v4_C1,
        vROS_C1 ~ vROSIf + vROSIq,
        vNADH_C1 ~ -0.5 * (v1_C1 + v2_C1 + v3_C1),
        D(SQ_C1) ~ v1_C1 - v2_C1 - v4_C1,
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
prob = ODEProblem(sys, [], tend)

sol = solve(prob)

sol[sys.fFMNH_C1][end]
sol[sys.fFMNsq_C1][end]
sol[sys.KEQ1_C1][end]
sol[sys.KEQ2_C1][end]
sol[sys.SQ_C1][end]
sol[sys.vROS_C1][end]

@unpack v1_C1, v2_C1, v3_C1, v4_C1 = sys
@unpack vQ_C1, vQH2_C1, vROS_C1, vNADH_C1 = sys
plot(sol, idxs=[vQ_C1, vQH2_C1, vROS_C1, vNADH_C1], tspan=(20, 40))
