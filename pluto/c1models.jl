### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 6e998cd6-a007-45aa-b05b-4439da7c371d
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	using ModelingToolkit
	using PlutoUI
	using OrdinaryDiffEq
	using SteadyStateDiffEq
	using Plots
	using NaNMath
	using DisplayAs: PNG
	md"Import packages"
end

# ╔═╡ edece65e-8635-11f0-39b5-c5c17d0ad6f0
md"""
# Complex I model

Comparing Gauthier, Markevich, and simplified complex I models
"""

# ╔═╡ a43f0927-bc4d-45d4-bdf8-626de501c8cf
PlutoUI.TableOfContents()

# ╔═╡ 9e16194f-f20c-41f5-a271-b8aee964b0df
begin
	const mM = 1
	const Molar = 1000mM
	const μM = mM / 1000
	const second = 1
	const Hz = inv(second)
	const ms = second / 1000
	const minute = 60second
	const Volt = 1
	const mV = Volt / 1000
	const RGAS = 8.314 		# Joule / Kelvin / mol
	const Temp = 310 		# Kelvin
	const Faraday = 96485  	# Columb / mol
	const iVT = Faraday / RGAS / Temp
	md"Constants"
end

# ╔═╡ f01351ba-e5d6-42e1-90b3-aed341352f3e
# Gauthier 2012 7-state QSSA model
function c1_gauthier(; name=:c1gauthier,
    Q_n=1800μM, QH2_n=200μM, nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.01μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    C1_INHIB=1)

	@independent_variables t
	D = Differential(t)
	
    @parameters begin
        ET_C1 = 8.85mM      ## Activity of complex I
        dpsi_B_C1 = 50mV    ## Phase boundary potential
        ## Transition rates
        K12_C1 = 6.3396E11Hz / mM^2
        K21_C1 = 5Hz
        K56_C1 = 100Hz
        K65_C1 = 2.5119E13Hz / mM^2
        K61_C1 = 1e7Hz
        K16_C1 = 130Hz
        K23_C1 = 3886.7Hz / sqrt(mM)
        K32_C1 = 9.1295e6Hz
        K34_C1 = 639.1364Hz
        K43_C1 = 3.2882Hz / sqrt(mM)
        K47_C1 = 1.5962E7Hz / mM
        K74_C1 = 65.2227Hz
        K75_C1 = 24615Hz
        K57_C1 = 1166.7Hz / sqrt(mM)
        K42_C1 = 6.0318Hz / mM
        Em_O2_SOX = -160mV         ## O2/Superoxide redox potential
        Em_FMNH2_FMNH = -375mV     ## FMNH/FMNH2 redox potential
        rKEQ_ROS_C1 = exp(iVT * (Em_FMNH2_FMNH - Em_O2_SOX))
    end

    @variables begin
        C1_1(t)
        C1_2(t)
        C1_3(t)
        C1_4(t)
        C1_5(t)
        C1_6(t)
        C1_7(t)
        wC1_1(t)
        wC1_2(t)
        wC1_3(t)
        wC1_4(t)
        wC1_5(t)
        wC1_6(t)
        wC1_7(t)
        vQC1(t)
        vNADHC1(t)
        vROSC1(t)
        vHresC1(t)
        TNC1(t) ## Turnover number
    end

    fv = exp(iVT * (dpsi - dpsi_B_C1))
    ## State transition rates
    a12 = K12_C1 * h_m^2
    a21 = K21_C1
    a65 = K65_C1 * h_i^2
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
    a24 = K42_C1 * rKEQ_ROS_C1 * sox_m

    ## KA pattern
    w1 = a21 * a32 * a42 * a56 * a61 * a74 + a21 * a32 * a42 * a56 * a61 * a75 + a21 * a32 * a42 * a57 * a61 * a74 + a21 * a32 * a42 * a57 * a65 * a74 + a21 * a32 * a43 * a56 * a61 * a74 + a21 * a32 * a43 * a56 * a61 * a75 + a21 * a32 * a43 * a57 * a61 * a74 + a21 * a32 * a43 * a57 * a65 * a74 + a21 * a32 * a47 * a56 * a61 * a75 + a21 * a34 * a42 * a56 * a61 * a74 + a21 * a34 * a42 * a56 * a61 * a75 + a21 * a34 * a42 * a57 * a61 * a74 + a21 * a34 * a42 * a57 * a65 * a74 + a21 * a34 * a47 * a56 * a61 * a75 + a23 * a34 * a47 * a56 * a61 * a75 + a24 * a32 * a47 * a56 * a61 * a75 + a24 * a34 * a47 * a56 * a61 * a75
    w2 = a12 * a32 * a42 * a56 * a61 * a74 + a12 * a32 * a42 * a56 * a61 * a75 + a12 * a32 * a42 * a57 * a61 * a74 + a12 * a32 * a42 * a57 * a65 * a74 + a12 * a32 * a43 * a56 * a61 * a74 + a12 * a32 * a43 * a56 * a61 * a75 + a12 * a32 * a43 * a57 * a61 * a74 + a12 * a32 * a43 * a57 * a65 * a74 + a12 * a32 * a47 * a56 * a61 * a75 + a12 * a34 * a42 * a56 * a61 * a74 + a12 * a34 * a42 * a56 * a61 * a75 + a12 * a34 * a42 * a57 * a61 * a74 + a12 * a34 * a42 * a57 * a65 * a74 + a12 * a34 * a47 * a56 * a61 * a75 + a16 * a32 * a42 * a57 * a65 * a74 + a16 * a32 * a43 * a57 * a65 * a74 + a16 * a34 * a42 * a57 * a65 * a74
    w3 = a12 * a23 * a42 * a56 * a61 * a74 + a12 * a23 * a42 * a56 * a61 * a75 + a12 * a23 * a42 * a57 * a61 * a74 + a12 * a23 * a42 * a57 * a65 * a74 + a12 * a23 * a43 * a56 * a61 * a74 + a12 * a23 * a43 * a56 * a61 * a75 + a12 * a23 * a43 * a57 * a61 * a74 + a12 * a23 * a43 * a57 * a65 * a74 + a12 * a23 * a47 * a56 * a61 * a75 + a12 * a24 * a43 * a56 * a61 * a74 + a12 * a24 * a43 * a56 * a61 * a75 + a12 * a24 * a43 * a57 * a61 * a74 + a12 * a24 * a43 * a57 * a65 * a74 + a16 * a21 * a43 * a57 * a65 * a74 + a16 * a23 * a42 * a57 * a65 * a74 + a16 * a23 * a43 * a57 * a65 * a74 + a16 * a24 * a43 * a57 * a65 * a74
    w4 = a12 * a23 * a34 * a56 * a61 * a74 + a12 * a23 * a34 * a56 * a61 * a75 + a12 * a23 * a34 * a57 * a61 * a74 + a12 * a23 * a34 * a57 * a65 * a74 + a12 * a24 * a32 * a56 * a61 * a74 + a12 * a24 * a32 * a56 * a61 * a75 + a12 * a24 * a32 * a57 * a61 * a74 + a12 * a24 * a32 * a57 * a65 * a74 + a12 * a24 * a34 * a56 * a61 * a74 + a12 * a24 * a34 * a56 * a61 * a75 + a12 * a24 * a34 * a57 * a61 * a74 + a12 * a24 * a34 * a57 * a65 * a74 + a16 * a21 * a32 * a57 * a65 * a74 + a16 * a21 * a34 * a57 * a65 * a74 + a16 * a23 * a34 * a57 * a65 * a74 + a16 * a24 * a32 * a57 * a65 * a74 + a16 * a24 * a34 * a57 * a65 * a74
    w5 = a12 * a23 * a34 * a47 * a61 * a75 + a12 * a23 * a34 * a47 * a65 * a75 + a12 * a24 * a32 * a47 * a61 * a75 + a12 * a24 * a32 * a47 * a65 * a75 + a12 * a24 * a34 * a47 * a61 * a75 + a12 * a24 * a34 * a47 * a65 * a75 + a16 * a21 * a32 * a42 * a65 * a74 + a16 * a21 * a32 * a42 * a65 * a75 + a16 * a21 * a32 * a43 * a65 * a74 + a16 * a21 * a32 * a43 * a65 * a75 + a16 * a21 * a32 * a47 * a65 * a75 + a16 * a21 * a34 * a42 * a65 * a74 + a16 * a21 * a34 * a42 * a65 * a75 + a16 * a21 * a34 * a47 * a65 * a75 + a16 * a23 * a34 * a47 * a65 * a75 + a16 * a24 * a32 * a47 * a65 * a75 + a16 * a24 * a34 * a47 * a65 * a75
    w6 = a12 * a23 * a34 * a47 * a56 * a75 + a12 * a24 * a32 * a47 * a56 * a75 + a12 * a24 * a34 * a47 * a56 * a75 + a16 * a21 * a32 * a42 * a56 * a74 + a16 * a21 * a32 * a42 * a56 * a75 + a16 * a21 * a32 * a42 * a57 * a74 + a16 * a21 * a32 * a43 * a56 * a74 + a16 * a21 * a32 * a43 * a56 * a75 + a16 * a21 * a32 * a43 * a57 * a74 + a16 * a21 * a32 * a47 * a56 * a75 + a16 * a21 * a34 * a42 * a56 * a74 + a16 * a21 * a34 * a42 * a56 * a75 + a16 * a21 * a34 * a42 * a57 * a74 + a16 * a21 * a34 * a47 * a56 * a75 + a16 * a23 * a34 * a47 * a56 * a75 + a16 * a24 * a32 * a47 * a56 * a75 + a16 * a24 * a34 * a47 * a56 * a75
    w7 = a12 * a23 * a34 * a47 * a56 * a61 + a12 * a23 * a34 * a47 * a57 * a61 + a12 * a23 * a34 * a47 * a57 * a65 + a12 * a24 * a32 * a47 * a56 * a61 + a12 * a24 * a32 * a47 * a57 * a61 + a12 * a24 * a32 * a47 * a57 * a65 + a12 * a24 * a34 * a47 * a56 * a61 + a12 * a24 * a34 * a47 * a57 * a61 + a12 * a24 * a34 * a47 * a57 * a65 + a16 * a21 * a32 * a42 * a57 * a65 + a16 * a21 * a32 * a43 * a57 * a65 + a16 * a21 * a32 * a47 * a57 * a65 + a16 * a21 * a34 * a42 * a57 * a65 + a16 * a21 * a34 * a47 * a57 * a65 + a16 * a23 * a34 * a47 * a57 * a65 + a16 * a24 * a32 * a47 * a57 * a65 + a16 * a24 * a34 * a47 * a57 * a65

    den = wC1_1 + wC1_2 + wC1_3 + wC1_4 + wC1_5 + wC1_6 + wC1_7

    v47 = a47 * C1_4 - a74 * C1_7
    v42 = a42 * C1_4 - a24 * C1_2
    v23 = a23 * C1_2 - a32 * C1_3
    v61 = a61 * C1_1 - a16 * C1_6
    eqs = [
        wC1_1 ~ w1,
        wC1_2 ~ w2,
        wC1_3 ~ w3,
        wC1_4 ~ w4,
        wC1_5 ~ w5,
        wC1_6 ~ w6,
        wC1_7 ~ w7,
        C1_1 ~ wC1_1 / den * ET_C1,
        C1_2 ~ wC1_2 / den * ET_C1,
        C1_3 ~ wC1_3 / den * ET_C1,
        C1_4 ~ wC1_4 / den * ET_C1,
        C1_7 ~ wC1_5 / den * ET_C1,
        C1_5 ~ wC1_6 / den * ET_C1,
        C1_6 ~ wC1_7 / den * ET_C1,
        vQC1 ~ -0.5 * v47,
        vROSC1 ~ v42,
        vNADHC1 ~ -0.5 * v23,
        vHresC1 ~ 2 * v61,
        TNC1 ~ -vNADHC1 / ET_C1,
    ]
    return System(eqs, t; name)
end

# ╔═╡ 91e8801e-0aaa-454f-8a52-7fbf7f70072c
# Markevich 2015 mass action law model
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4426091/
function c1_markevich_full(; name=:c1markevich_full,
    Q_n=1.8mM, QH2_n=0.2mM,
    nad=500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.01μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0)

	@independent_variables t
	D = Differential(t)
	
    @parameters begin
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -150mV
        Em_N1a = -370mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich, 2015
        Em_SQ_QH2_C1 = +500mV     ## ~800mV (?) in Markevich, 2015
        ET_C1 = 17μM              ## Activity of complex I
        ## DOX IC50 on complex I
        KI_DOX_C1 = 400μM
        kf1_C1 = 83Hz / μM
        KEQ1_C1 = 0.01 / μM
        kr1_C1 = kf1_C1 / KEQ1_C1
        kf3_C1 = 1e6Hz
        KEQ3_C1 = 25μM
        kr3_C1 = kf3_C1 / KEQ3_C1
        kf2_C1 = 1.44e12Hz
        KEQ2_C1 = exp(2iVT * (Em_FMN_FMNH - Em_NAD)) / KEQ1_C1 / KEQ3_C1
        kr2_C1 = kf2_C1 / KEQ2_C1
        kf4_C1 = 1Hz / μM
        KEQ4_C1 = 0.001 / μM
        kr4_C1 = kf4_C1 / KEQ4_C1
        kf5_C1 = 2Hz / μM
        KEQ5_C1 = 0.02 / μM
        kr5_C1 = kf5_C1 / KEQ5_C1
        kf6_C1 = 5e8Hz / μM
        KEQ6_C1 = exp(iVT * (Em_N3 - Em_FMNsq_FMNH))
        kr6_C1 = kf6_C1 / KEQ6_C1
        kf7_C1 = 10000Hz / μM
        KEQ7_C1 = exp(iVT * (Em_N2 - Em_N3))
        kr7_C1 = kf7_C1 / KEQ7_C1
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = 0.1 / μM         ## Association constant for Q
        kr8_C1 = kf8_C1 / KEQ8_C1
        kf9_C1 = 4e5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kr9_C1 = kf9_C1 / KEQ9_C1
        kf10_C1 = 2e6Hz / μM
        KEQ10_C1 = exp(iVT * (Em_N1a - Em_FMN_FMNsq))
        KEQ10B_C1 = exp(iVT * (Em_FMNsq_FMNH - Em_N1a))
        kr10_C1 = kf10_C1 / KEQ10_C1
        kr10b_C1 = kf10_C1 / KEQ10B_C1
        kf11_C1 = 1e9Hz / μM
        KEQ11_C1 = exp(iVT * (Em_N3 - Em_FMN_FMNsq))
        kr11_C1 = kf11_C1 / KEQ11_C1
        kf13_C1 = 2.7e6Hz / μM
        kf14_C1 = 1000Hz
        KEQ14_C1 = 20μM             ## Dissociation constant for QH2
        kr14_C1 = kf14_C1 / KEQ14_C1
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kr16_C1 = kf16_C1 / KEQ16_C1
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        KEQ17_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
        kr17_C1 = kf17_C1 / KEQ17_C1
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
        vQC1(t)
        vQH2C1(t)
        vNADHC1(t)
        vNADC1(t)
        vROSIf(t)
        vROSIq(t)
        vROSC1(t)
        TNC1(t) ## NADH turnover number
        vHresC1(t)
    end

    fhm = h_m / 1E-7Molar
    C1_INHIB = (1 - ROTENONE_BLOCK)
    ## NADH + FMN = FMN.NADH
    v1 = kf1_C1 * nadh * FMN - kr1_C1 * FMN_NADH
    ## FMN.NADH = FMNH−.NAD+
    v2 = kf2_C1 * FMN_NADH - kr2_C1 * FMNH_NAD
    ## FMNH−.NAD+ = FMNH− + NAD+
    v3 = kf3_C1 * FMNH_NAD - kr3_C1 * FMNH * nad
    ## FMN + NAD+ = FMN.NAD+
    v4 = kf4_C1 * nad * FMN - kr4_C1 * FMN_NAD
    ## FMNH− + NADH = FMNH−.NADH
    v5 = kf5_C1 * FMNH * nadh - kr5_C1 * FMNH_NADH
    ## FMNH− + N3 = FMNHsq + N3−
    v6 = kf6_C1 * FMNH * N3_C1 - kr6_C1 * FMNsq * N3r_C1
    ## N3− + N2 = N3 + N2−
    v7 = kf7_C1 * N3r_C1 * N2_C1 - kr7_C1 * N3_C1 * N2r_C1
    ## Q association
    q = Q_n * C1_INHIB
    v8 = kf8_C1 * Iq_C1 * q - kr8_C1 * Q_C1
    ## CI.Q + N2− = CIQsq + N2
    v9 = kf9_C1 * Q_C1 * N2r_C1 - kr9_C1 * SQ_C1 * N2_C1
    ## FMNHsq + N1a = FMN + N1a− + Hi+
    v10 = kf10_C1 * FMNsq * N1a_C1 - kr10_C1 * FMN * N1ar_C1 * fhm
    ## FMNsq + N1a− = FMNH- + N1a
    v10b = kf10_C1 * FMNsq * N1ar_C1 - kr10b_C1 * FMNH * N1a_C1
    ## FMNHsq + N3 = FMN + N3− + Hi+
    v11 = kf11_C1 * FMNsq * N3_C1 - kr11_C1 * FMN * N3r_C1 * fhm
    ## N2 + N3− = N2− + N3
    v12 = v7
    ## Second electron transfer
    v13 = kf13_C1 * (SQ_C1 * N2r_C1 * fhm^2 - QH2_C1 * N2_C1 / KEQ13_C1)
    ## QH2 dissociation
    qh2 = QH2_n * C1_INHIB
    v14 = kf14_C1 * QH2_C1 - kr14_C1 * Iq_C1 * qh2
    ## Flavin site ROS generation
    v16 = kf16_C1 * FMNH * O2 - kr16_C1 * FMNsq * sox_m
    ## Quinone site ROS generation
    v17 = kf17_C1 * SQ_C1 * O2 - kr17_C1 * Q_C1 * sox_m

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
        D(FMNH) ~ v3 - v5 - v16 - v6 + v10b,
        D(FMNsq) ~ v6 + v16 - v10 - v11 - v10b,
        D(N1ar_C1) ~ v10 - v10b,
        D(N3r_C1) ~ v6 + v11 - v7 - v12,
        D(N2r_C1) ~ v7 + v12 - v9 - v13,
        D(Q_C1) ~ v8 - v9 + v17,
        D(SQ_C1) ~ v9 - v17 - v13,
        D(QH2_C1) ~ v13 - v14,

        ## Positive: production; negative: consumption
        vNADHC1 ~ -(v1 + v5),
        vNADC1 ~ v3 - v4,
        vQC1 ~ -v8,
        vQH2C1 ~ v14,
        vROSIf ~ v16,
        vROSIq ~ v17,
        vHresC1 ~ 4 * v13,
        vROSC1 ~ vROSIf + vROSIq,
        TNC1 ~ -vNADHC1 / ET_C1,
    ]
    return System(eqs, t; name)
end

# ╔═╡ 68cf8251-063d-4130-abef-a931172d6b1f
# Simplified Markevich complex I model
# Rapid equlibrium in the flavin site
# QSSA for the catalytic cycle in the quinone site
function c1_markevich_s(; name=:c1s,
    Q_n=1800μM, QH2_n=200μM,
    nad=2500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.01μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0, MT_PROT=1)

	@independent_variables t
	D = Differential(t)

    @parameters begin
        ET_C1 = 17μM                ## Activity of complex I
        KI_DOX_C1 = 400μM           ## DOX IC50 on complex I
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -150mV
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
        vQC1(t)
        vQH2C1(t)
        vROSC1(t)
        vROSIf(t)
        vROSIq(t)
        vNADHC1(t)
        vNADC1(t)
        TNC1(t)
        vHresC1(t)
    end

    C1_CONC = ET_C1 * MT_PROT
    C1_INHIB = (1 - ROTENONE_BLOCK) / (1 + (DOX / KI_DOX_C1)^3)
    ## Mitochondrial pH factor
    fhm = h_m * inv(1E-7Molar)

    ## Flavin site in rapid equilibrium
    ## Weights in the flavin site
    wFMN = 1
    wFMN_NAD = wFMN * nad / KI_NAD_C1
    wFMN_NADH = wFMN * nadh / KD_NADH_C1
    wFMNH = wFMN * (nadh / nad) * KEQ_NADH_FMN
    wFMNH_NAD = wFMNH * nad / KD_NAD_C1
    wFMNH_NADH = wFMNH * nadh / KI_NADH_C1
    wFMNsq = NaNMath.sqrt(wFMN * wFMNH * rKEQ_FMNsq_Dis * fhm)
    fDen = wFMN + wFMN_NAD + wFMNH + wFMNH_NADH + wFMNsq + wFMN_NADH + wFMNH_NAD
    fC1 = C1_CONC / fDen

    ## Flavin site ROS generation
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
    wC1 = b21 * b32 * b41 + b21 * b32 * b43 + b21 * b34 * b41 + b23 * b34 * b41
    wC1_Q = b12 * b32 * b41 + b12 * b32 * b43 + b12 * b34 * b41 + b14 * b32 * b43
    wC1_SQ = b12 * b23 * b41 + b12 * b23 * b43 + b14 * b21 * b43 + b14 * b23 * b43
    wC1_QH2 = b12 * b23 * b34 + b14 * b21 * b32 + b14 * b21 * b34 + b14 * b23 * b34
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
    return System(eqs, t; name)
end

# ╔═╡ a793d5fe-64f8-4fbe-9297-df7ee8ba96e1
# Quinone site complex I model, assuming N2 at equlibrium with the flavin site
# Adapted from the simplified Markevich complex I model
function c1_q(; name=:c1q,
    Q_n=1800μM, QH2_n=200μM,
    nad=2500μM, nadh=500μM,
    dpsi=150mV, O2=6μM, sox_m=0.01μM,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    DOX=0μM, ROTENONE_BLOCK=0, MT_PROT=1)

	@independent_variables t
	D = Differential(t)

    @parameters begin
        ET_C1 = 17μM              ## Activity of complex I
        KI_DOX_C1 = 400μM         ## DOX IC50 on complex I
        Em_O2_SOX = -160mV        ## O2/Superoxide redox potential
        Em_FMN_FMNsq = -387mV     ## FMN/FMNH- avg redox potential
        Em_FMNsq_FMNH = -293mV    ## FMN semiquinone/FMNH- redox potential
        Em_FMN_FMNH = -340mV      ## FMN/FMNH- avg redox potential
        Em_NAD = -320mV           ## NAD/NADH avg redox potential
        Em_N3 = -250mV
        Em_N2 = -150mV
        Em_Q_SQ_C1 = -300mV       ## -213mV in Markevich's model
        Em_SQ_QH2_C1 = +500mV
        KI_NADH_C1 = 50μM
        KD_NADH_C1 = 100μM
        KI_NAD_C1 = 1000μM
        KD_NAD_C1 = 25μM
        ## NADH + FMN = NAD+ + FMNH-
        KEQ_NADH_FMN = exp(2iVT * (Em_FMN_FMNH - Em_NAD))
        ## 2FMNsq = (ISC) = FMN + FMNH- + H+
        rKEQ_FMNsq_Dis = exp(-iVT * (Em_FMNsq_FMNH - Em_FMN_FMNsq))
        ## FMNH- + N2 = FMNsq + N2-
        KEQ_FMNH_N2 = exp(iVT * (Em_N2 - Em_FMNsq_FMNH))
        ## I + Q = IQ
        kf8_C1 = 10Hz / μM
        KEQ8_C1 = 0.1 / μM
        kr8_C1 = kf8_C1 / KEQ8_C1
        ## Q + e- = Q-
        kf9_C1 = 1e4Hz / μM ## 4e5Hz / μM
        KEQ9_C1 = exp(iVT * (Em_Q_SQ_C1 - Em_N2))
        kr9_C1 = kf9_C1 / KEQ9_C1
        ## Q- + e- + 6Hm = QH2 + 4Hi
        kf13_C1 = 6e4Hz / μM ## 2.7e6Hz / μM
        ## C1_QH2 = C1 + QH2
        kf14_C1 = 1000Hz
        KEQ14_C1 = 20μM
        kr14_C1 = kf14_C1 / KEQ14_C1 ## 50Hz / μM
        kf16_C1 = 2Hz / μM          ## SOX production rate from If site
        KEQ16_C1 = exp(iVT * (Em_O2_SOX - Em_FMNsq_FMNH))
        kr16_C1 = kf16_C1 / KEQ16_C1
        kf17_C1 = 0.04Hz / μM       ## SOX production rate from Iq site
        KEQ17_C1 = exp(iVT * (Em_O2_SOX - Em_Q_SQ_C1))
        kr17_C1 = kf17_C1 / KEQ17_C1
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
        I_C1(t) ## Conserved
        Q_C1(t) = 0
        SQ_C1(t) = 0
        QH2_C1(t) = 0
        rKEQ_N2r_SQ(t)
        ## Reaction rates
        vQC1(t)
        vQH2C1(t)
        vROSC1(t)
        vROSIf(t)
        vROSIq(t)
        vNADHC1(t)
        vNADC1(t)
        TNC1(t)
        vHresC1(t)
    end

    C1_CONC = ET_C1 * MT_PROT
    C1_INHIB = (1 - ROTENONE_BLOCK) / (1 + (DOX / KI_DOX_C1)^3)
    fhm = h_m * inv(1E-7Molar) ## Mitochondrial pH factor
    ## Flavin site in rapid equilibrium
    ## Weights in the flavin site
    wFMN = 1
    wFMN_NAD = wFMN * nad / KI_NAD_C1
    wFMN_NADH = wFMN * nadh / KD_NADH_C1
    wFMNH = wFMN * (nadh / nad) * KEQ_NADH_FMN
    wFMNH_NAD = wFMNH * nad / KD_NAD_C1
    wFMNH_NADH = wFMNH * nadh / KI_NADH_C1
    wFMNsq = NaNMath.sqrt(wFMN * wFMNH * rKEQ_FMNsq_Dis * fhm)
    fDen = wFMN + wFMN_NAD + wFMNH + wFMNH_NADH + wFMNsq + wFMN_NADH + wFMNH_NAD
    fC1 = C1_CONC / fDen

    ## Quinone site state transitions
    ## Q + C1 = Q_C1
    q = Q_n * C1_INHIB
    v8 = kf8_C1 * I_C1 * q - kr8_C1 * Q_C1
    ## Q_C1 + N2r = SQ_C1 + N2
    v9 = kf9_C1 * Q_C1 * N2r_C1 - kr9_C1 * SQ_C1 * N2_C1
    ## SQ_C1 + N2r + 6Hm = QH2_C1 + N2 + 4Hi
    v13 = kf13_C1 * (fhm^2 * N2r_C1 * SQ_C1 - rKEQ_N2r_SQ * N2_C1 * QH2_C1)
    ## QH2_C1 = QH2 + C1
    qh2 = QH2_n * C1_INHIB
    v14 = kf14_C1 * QH2_C1 - kr14_C1 * I_C1 * qh2
    ## Flavin site ROS production
    v16 = kf16_C1 * FMNH * O2 - kr16_C1 * FMNsq * sox_m
    ## Quinone site ROS production
    v17 = kf17_C1 * SQ_C1 * O2 - kr17_C1 * Q_C1 * sox_m

    eqs = [
        rKEQ_N2r_SQ ~ exp(-iVT * (Em_SQ_QH2_C1 - Em_N2 - 4dpsi)) * (h_i / h_m)^4,
        FMN ~ wFMN * fC1,
        FMN_NAD ~ wFMN_NAD * fC1,
        FMNH ~ wFMNH * fC1,
        FMNsq ~ wFMNsq * fC1,
        FMNH_NADH ~ wFMNH_NADH * fC1,
        FMN_NADH ~ wFMN_NADH * fC1,
        FMNH_NAD ~ wFMNH_NAD * fC1,
        N2_C1 ~ C1_CONC * FMNsq / (FMNsq + FMNH * KEQ_FMNH_N2),
        C1_CONC ~ N2r_C1 + N2_C1,
        C1_CONC ~ I_C1 + Q_C1 + SQ_C1 + QH2_C1,
        D(Q_C1) ~ v8 - v9 + v17,
        D(SQ_C1) ~ v9 - v13 - v17,
        D(QH2_C1) ~ v13 - v14,
        vNADHC1 ~ -0.5 * (v9 + v13 + v16),
        vROSIf ~ v16,
        vROSIq ~ v17,
        vROSC1 ~ vROSIf + vROSIq,
        vQH2C1 ~ v14,
        vQC1 ~ -v8,
        vHresC1 ~ 4 * v13,
        vNADC1 ~ -vNADHC1,
        TNC1 ~ vNADC1 / C1_CONC,
    ]
    return System(eqs, t; name)
end

# ╔═╡ a235d169-0f16-442b-b8d3-0c6990e702d0
@parameters begin
    Q_n = 1.8mM
    QH2_n = 0.2mM
    nad = 2500μM
    nadh = 500μM
    dpsi = 150mV
    sox_m = 0.01μM
end

# ╔═╡ a5cdf926-3462-4445-b647-cb808486f242
# Helper function
extract(sim, k) = map(s -> s[k], sim)

# ╔═╡ 848ab4e3-d35a-4696-97e1-c78a24f1b400
sys_q = c1_q(; Q_n, QH2_n, nad, nadh, dpsi, sox_m) |> mtkcompile

# ╔═╡ d20a3185-b1d8-4eca-ac3a-b5cad7ce929a
sys_m = c1_markevich_full(; Q_n, QH2_n, nad, nadh, dpsi, sox_m) |> mtkcompile

# ╔═╡ a588058f-cbbc-4e71-ac78-f9f78648586e
sys_g = c1_gauthier(; Q_n, QH2_n, nad, nadh, dpsi, sox_m) |> mtkcompile

# ╔═╡ 144490e3-7e81-4ea4-946c-1198acfabe47
# The parameter for ROS generation rate is adjusted (x10000) to be comparable to ROS generation from complex III
prob_g = SteadyStateProblem(sys_g, [
	sys_g.K42_C1 => 6.0318Hz / mM * 10000
])

# ╔═╡ 14fc9ce7-085d-4679-a868-9f2d2cb8124a
prob_q = SteadyStateProblem(sys_q, [
    sys_q.ET_C1 => 17μM,
    sys_q.kf8_C1 => 5Hz / μM,
    sys_q.kf9_C1 => 1000Hz / μM,
    sys_q.kf13_C1 => 10000Hz / μM,
    sys_q.kf16_C1 => 20Hz / μM,
    sys_q.kf17_C1 => 0.4Hz / μM,
])

# ╔═╡ 10c4c6a4-adba-43f9-8b19-0a5af2bd6ceb
# The redox potential of N2 is adjusted to -150mV in B. taurus mitochondria
# instead of -80mV in E. coli(?)
prob_m = SteadyStateProblem(sys_m, [
    sys_m.ET_C1 => 17μM,
    sys_m.kf16_C1 => 20Hz / μM,
    sys_m.kf17_C1 => 0.4Hz / μM,
])

# ╔═╡ 25f1b5ca-fa1b-4d3b-ac5b-53a51d025eae
alg = DynamicSS(TRBDF2())

# ╔═╡ 0f7aa8de-c7fc-4b00-b76f-e598e7ea1beb
ealg = EnsembleThreads()

# ╔═╡ 80b2e8ac-857b-4f5d-883a-89587d06ff13
md"""
## Varying MMP
"""

# ╔═╡ 1fca3b81-1811-4596-9640-df0bac1907c2
dpsirange = 100:1:200

# ╔═╡ 9d52ce00-9edb-42f1-8036-ea8262623d74
alter_dpsi = (prob, i, repeat) -> begin
    prob.ps[dpsi] = dpsirange[i] * mV
    prob
end

# ╔═╡ 1d7b6c29-2f8c-4115-b76c-c71c4fbe63ae
@time mmp_q = solve(EnsembleProblem(prob_q; prob_func=alter_dpsi), alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 7ab2c8ac-4376-4fa4-8f80-b1f91547d5f1
@time mmp_m = solve(EnsembleProblem(prob_m; prob_func=alter_dpsi), alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 2d841c84-a71c-4cee-8557-03c0e8e6faff
@time mmp_g = solve(EnsembleProblem(prob_g; prob_func=alter_dpsi), alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 95765fe6-a0ee-4938-8a1e-c968839a8500
md"""
### MMP vs NADH turnover
"""

# ╔═╡ 0a187444-73b1-48fb-8a6a-5a59d1f2ef8a
let
	xs = dpsirange
	ys = hcat(extract(mmp_g, sys_g.vNADHC1), extract(mmp_m, sys_m.vNADHC1), extract(mmp_q, sys_q.vNADHC1))
	plot(xs, ys, xlabel="MMP (mV)", ylabel="NADH rate (mM/s)", label=["Gauthier" "Markevich" "Q site"]) 
end |> PNG

# ╔═╡ 681583a3-f37a-4eb2-aa68-3a8686e3de3d
md"""
### MMP vs Q site status
"""

# ╔═╡ 0d2f1e9f-e740-47f6-998c-f9da661dff10
let
	xs = dpsirange
	ys = stack(extract.(Ref(mmp_q), [sys_q.Q_C1, sys_q.SQ_C1, sys_q.QH2_C1, sys_q.I_C1]), dims=2)
	pl1 = plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "I_C1"], title="Q site model", legend=:left, ylims=(0, 17μM))
	ys2 = stack(extract.(Ref(mmp_m), [sys_m.Q_C1, sys_m.SQ_C1, sys_m.QH2_C1, sys_m.Iq_C1]), dims=2)
	pl2 = plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "I_C1"], title="M model", ylims=(0, 17μM))
	plot(pl1, pl2)
end |> PNG

# ╔═╡ 7c7a793a-befb-443b-9af9-051a217b5646
md"""
### MMP vs F site status
"""

# ╔═╡ a3b060b3-b2b8-439e-9659-77a85e5ef87c
let
	xs = dpsirange
	ys = stack(extract.(Ref(mmp_q), [sys_q.FMN, sys_q.FMNsq, sys_q.FMNH, sys_q.FMN_NAD, sys_q.FMNH_NADH, sys_q.FMN_NADH, sys_q.FMNH_NAD]), dims=2)
	pl1 = plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH" "FMN_NADH" "FMNH_NAD"], title="Q site model", legend=:right)
	
	ys = stack(extract.(Ref(mmp_m), [sys_m.FMN, sys_m.FMNsq, sys_m.FMNH, sys_m.FMN_NAD, sys_m.FMNH_NADH, sys_m.FMN_NADH, sys_m.FMNH_NAD]), dims=2)
	pl2 = plot(xs, ys, xlabel="MMP (mV)", ylabel="Concentration", label=["FMN" "FMNsq" "FMNH" "FMN_NAD" "FMNH_NADH" "FMN_NADH" "FMNH_NAD"], title="M model", legend=:right)
	
	plot(pl1, pl2) 
end |> PNG

# ╔═╡ 077bc541-a191-4171-8a30-a7881ece3f7e
md"""
### MMP vs ROS production
"""

# ╔═╡ 44ff5991-1b2f-4a15-b767-2a1ab95a6296
let
	xs = dpsirange
	ys_g = extract(mmp_g, sys_g.vROSC1)
	ys_m = extract(mmp_m, sys_m.vROSC1)
	ys_q = extract(mmp_q, sys_q.vROSC1)
	plot(xs, [ys_g ys_m ys_q], xlabel="MMP (mV)", ylabel="ROS production (mM/s)", label=["Gauthier" "Markevich" "Q site"])
end |> PNG

# ╔═╡ 646b8cbb-24f3-4319-a537-2e3a645d97bb
md"""
## Varying NADH
"""

# ╔═╡ 6e7724bf-2a7a-4857-98a3-ff4b573a046f
nadhrange = 10:10:2990

# ╔═╡ 4ecd3d6f-08a3-47fe-8b0c-245d6b1ba0f2
alter_nadh = (prob, i, repeat) -> begin
    prob.ps[nadh] = nadhrange[i] * μM
    prob.ps[nad] = 3000μM - prob.ps[nadh]
    prob
end

# ╔═╡ 0a0ed00d-e9cc-45a7-9ac8-d9ac342ea0d4
@time nad_q = solve(EnsembleProblem(prob_q; prob_func=alter_nadh), alg, ealg; trajectories=length(nadhrange), abstol=1e-8, reltol=1e-8)

# ╔═╡ b1bc0a77-4176-4865-a35a-4054b0c00050
@time nad_g = solve(EnsembleProblem(prob_g; prob_func=alter_nadh), alg, ealg; trajectories=length(nadhrange), abstol=1e-8, reltol=1e-8)

# ╔═╡ f149b683-50b1-435e-944f-f7155c0b26ce
@time nad_m = solve(EnsembleProblem(prob_m; prob_func=alter_nadh), alg, ealg; trajectories=length(nadhrange), abstol=1e-8, reltol=1e-8)

# ╔═╡ a0756738-fed4-40e9-bd80-1c35b7a3340f
md"""
### NADH vs turnover
"""

# ╔═╡ ab4eb5e8-e97d-4221-8400-48935532a12c
let
	xs = nadhrange
	ys_g = extract(nad_g, sys_g.vNADHC1)
	ys_m = extract(nad_m, sys_m.vNADHC1)
	ys_q = extract(nad_q, sys_q.vNADHC1)
	
	plot(xs, [ys_g ys_m ys_q], xlabel="NADH (μM)", ylabel="NADH rate (mM/s)", label=["Gauthier" "Markevich" "Q site"]) 
end |> PNG

# ╔═╡ cdd88f99-d491-479d-a44d-c7d164af0792
let
	xs = nadhrange
	ys = [extract(nad_m, sys_m.N2r_C1) extract(nad_q, sys_q.N2r_C1)] ./ μM
	plot(xs, ys, xlabel="NADH (μM)", ylabel="Reduced N2 (μM)", label=["Markevich" "Q site"])
end |> PNG

# ╔═╡ b59c9225-b585-4572-a277-0fef35ceb4f5
md"""
### NADH vs ROS production
"""

# ╔═╡ 92583a37-230a-4384-9d98-7f2deb8c6738
let
	xs = nadhrange
	ys = [extract(nad_g, sys_g.vROSC1) extract(nad_m, sys_m.vROSC1) extract(nad_q, sys_q.vROSC1)]
	plot(xs, ys, xlabel="NADH (μM)", ylabel="ROS production (mM/s)", label=["Gauthier" "Markevich" "Q site"])
end |> PNG

# ╔═╡ 95b38367-025f-46d7-a695-0ce3cb4c4228
let
	xs = nadhrange
	ys = stack(extract.(Ref(nad_q), [sys_q.Q_C1, sys_q.SQ_C1, sys_q.QH2_C1, sys_q.I_C1]), dims=2) .* 1000
	pl1 = plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration (μM)", label=["Q_C1" "SQ_C1" "QH2_C1" "I_C1"], title="Q site", legend=:top, ylims=(0, 17))
	
	ys = stack(extract.(Ref(nad_m), [sys_m.Q_C1, sys_m.SQ_C1, sys_m.QH2_C1, sys_m.Iq_C1]), dims=2) .* 1000
	pl2 = plot(xs, ys, xlabel="NADH (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "I_C1"], title="M model", ylims=(0, 17))
	
	plot(pl1, pl2)
end |> PNG

# ╔═╡ 2e506373-1fc3-4494-b155-17f4b8228693
md"""
## Varying Q/QH2
"""

# ╔═╡ d994343f-d56f-4777-9582-69a3459564f6
qh2range = 10:10:1990

# ╔═╡ 53ef8337-eeb8-4560-ab64-8ed88cb9aa01
alter_qh2 = (prob, i, repeat) -> begin
    prob.ps[QH2_n] = qh2range[i] * μM
    prob.ps[Q_n] = 2000μM - prob.ps[QH2_n]
    prob
end

# ╔═╡ 05f56cd1-862d-4f83-b617-d368b94ecc37
@time qh2_q = solve(EnsembleProblem(prob_q; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 72135403-f863-4b1a-8e8b-25cac207712a
@time qh2_g = solve(EnsembleProblem(prob_g; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 531937c3-9334-4860-8ea0-b4ccb1d10b58
@time qh2_m = solve(EnsembleProblem(prob_m; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 44554235-3a11-4f42-94b5-1c65a3ae500a
md"""
### QH2 vs turnover
"""

# ╔═╡ f99bdd9c-bbab-4630-9c04-dd58730e4811
let
	xs = qh2range
	ys = [extract(qh2_g, sys_g.vNADHC1) extract(qh2_m, sys_m.vNADHC1) extract(qh2_q, sys_q.vNADHC1)]
	plot(xs, ys, xlabel="QH2 (μM)", ylabel="NADH rate (mM/s)", label=["Gauthier" "Markevich" "Q site"])
end |> PNG

# ╔═╡ ba831dcf-3986-4c2d-8687-f3f31c00dd8f
md"""
### QH2 vs ROS production
"""

# ╔═╡ f2d2545b-1660-481a-abe4-9adb0ca919e3
let
	xs = qh2range
	ys = [extract(qh2_g, sys_g.vROSC1) extract(qh2_m, sys_m.vROSC1) extract(qh2_q, sys_q.vROSC1)]
	plot(xs, ys, xlabel="QH2 (μM)", ylabel="ROS production", label=["Gauthier" "Markevich" "Q site"])
end |> PNG

# ╔═╡ fd238756-307d-42fa-87fd-e7d5d41f670b
let
	xs = qh2range
	ys = stack(extract.(Ref(qh2_q), [sys_q.Q_C1, sys_q.SQ_C1, sys_q.QH2_C1, sys_q.I_C1]), dims=2) .* 1000
	pl1 = plot(xs, ys, xlabel="QH2 (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "I_C1"], title="Q site", legend=:top, ylims=(0, 17))
	
	ys = stack(extract.(Ref(qh2_m), [sys_m.Q_C1, sys_m.SQ_C1, sys_m.QH2_C1, sys_m.Iq_C1]), dims=2) .* 1000
	pl2 = plot(xs, ys, xlabel="QH2 (μM)", ylabel="Concentration", label=["Q_C1" "SQ_C1" "QH2_C1" "I_C1"], title="Markevich", legend=:top, ylims=(0, 17))
	
	plot(pl1, pl2)
end |> PNG

# ╔═╡ 2772584f-c7d3-4051-bc26-700c75b97070
md"""
Gauthier model is sensitive to high QH2 because of slow NAD binding.
"""

# ╔═╡ d38e6ec8-1ca5-4c5a-8de5-d51149c0ce72
let
	@unpack C1_1, C1_3, C1_4, C1_5, C1_6, C1_7 = sys_g
	xs = qh2range
	ys = stack(extract.(Ref(qh2_g), [C1_3, C1_4]), dims=2) .* 1000
	plot(xs, ys, xlabel="QH2 (μM)", ylabel="Conc (μM)", label=["C1_3" "C1_4"], lw=1.5)
end |> PNG

# ╔═╡ Cell order:
# ╠═edece65e-8635-11f0-39b5-c5c17d0ad6f0
# ╠═6e998cd6-a007-45aa-b05b-4439da7c371d
# ╠═a43f0927-bc4d-45d4-bdf8-626de501c8cf
# ╠═9e16194f-f20c-41f5-a271-b8aee964b0df
# ╟─f01351ba-e5d6-42e1-90b3-aed341352f3e
# ╠═91e8801e-0aaa-454f-8a52-7fbf7f70072c
# ╠═68cf8251-063d-4130-abef-a931172d6b1f
# ╠═a793d5fe-64f8-4fbe-9297-df7ee8ba96e1
# ╠═a235d169-0f16-442b-b8d3-0c6990e702d0
# ╠═a5cdf926-3462-4445-b647-cb808486f242
# ╠═848ab4e3-d35a-4696-97e1-c78a24f1b400
# ╠═d20a3185-b1d8-4eca-ac3a-b5cad7ce929a
# ╠═a588058f-cbbc-4e71-ac78-f9f78648586e
# ╠═144490e3-7e81-4ea4-946c-1198acfabe47
# ╠═14fc9ce7-085d-4679-a868-9f2d2cb8124a
# ╠═10c4c6a4-adba-43f9-8b19-0a5af2bd6ceb
# ╠═25f1b5ca-fa1b-4d3b-ac5b-53a51d025eae
# ╠═0f7aa8de-c7fc-4b00-b76f-e598e7ea1beb
# ╠═80b2e8ac-857b-4f5d-883a-89587d06ff13
# ╠═1fca3b81-1811-4596-9640-df0bac1907c2
# ╠═9d52ce00-9edb-42f1-8036-ea8262623d74
# ╠═1d7b6c29-2f8c-4115-b76c-c71c4fbe63ae
# ╠═7ab2c8ac-4376-4fa4-8f80-b1f91547d5f1
# ╠═2d841c84-a71c-4cee-8557-03c0e8e6faff
# ╠═95765fe6-a0ee-4938-8a1e-c968839a8500
# ╠═0a187444-73b1-48fb-8a6a-5a59d1f2ef8a
# ╠═681583a3-f37a-4eb2-aa68-3a8686e3de3d
# ╠═0d2f1e9f-e740-47f6-998c-f9da661dff10
# ╠═7c7a793a-befb-443b-9af9-051a217b5646
# ╠═a3b060b3-b2b8-439e-9659-77a85e5ef87c
# ╠═077bc541-a191-4171-8a30-a7881ece3f7e
# ╠═44ff5991-1b2f-4a15-b767-2a1ab95a6296
# ╠═646b8cbb-24f3-4319-a537-2e3a645d97bb
# ╠═6e7724bf-2a7a-4857-98a3-ff4b573a046f
# ╠═4ecd3d6f-08a3-47fe-8b0c-245d6b1ba0f2
# ╠═0a0ed00d-e9cc-45a7-9ac8-d9ac342ea0d4
# ╠═b1bc0a77-4176-4865-a35a-4054b0c00050
# ╠═f149b683-50b1-435e-944f-f7155c0b26ce
# ╠═a0756738-fed4-40e9-bd80-1c35b7a3340f
# ╠═ab4eb5e8-e97d-4221-8400-48935532a12c
# ╠═cdd88f99-d491-479d-a44d-c7d164af0792
# ╠═b59c9225-b585-4572-a277-0fef35ceb4f5
# ╠═92583a37-230a-4384-9d98-7f2deb8c6738
# ╠═95b38367-025f-46d7-a695-0ce3cb4c4228
# ╠═2e506373-1fc3-4494-b155-17f4b8228693
# ╠═d994343f-d56f-4777-9582-69a3459564f6
# ╠═53ef8337-eeb8-4560-ab64-8ed88cb9aa01
# ╠═05f56cd1-862d-4f83-b617-d368b94ecc37
# ╠═72135403-f863-4b1a-8e8b-25cac207712a
# ╠═531937c3-9334-4860-8ea0-b4ccb1d10b58
# ╠═44554235-3a11-4f42-94b5-1c65a3ae500a
# ╠═f99bdd9c-bbab-4630-9c04-dd58730e4811
# ╠═ba831dcf-3986-4c2d-8687-f3f31c00dd8f
# ╠═f2d2545b-1660-481a-abe4-9adb0ca919e3
# ╠═fd238756-307d-42fa-87fd-e7d5d41f670b
# ╠═2772584f-c7d3-4051-bc26-700c75b97070
# ╠═d38e6ec8-1ca5-4c5a-8de5-d51149c0ce72
