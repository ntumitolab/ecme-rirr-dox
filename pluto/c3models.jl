### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ acd4eb16-255d-46ee-b498-7e9787837492
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

# ╔═╡ 13270f16-84b6-11f0-033c-a9652a6e1210
md"""
# Complex III model

Comparing Gauthier's and my complex III models

Reference redox potetials:

| Species    	| ``E_m`` (mV) |
| -------- 		| ------- 	|
| Q/QH2(p)  	| +60~65  	|
| Q/SQ(p) 		| <-170   	|
| SQ/QH2(p)    	| >290    	|
| Q/SQ(n)    	| +50   	|
| SQ/QH2(n)    	| +150    	|
| bL    		| -40    	|
| bH    		| +20~40    |
| ISP    		| +280    	|
| CytC1    		| +245    	|
| CytC    		| +265    	|
| O2/SOX    	| -160    	|

"""

# ╔═╡ 9061f18a-31ad-46ac-b0fa-aa632d63147c
PlutoUI.TableOfContents()

# ╔═╡ a6e0a7c4-db17-4041-94d3-ec0ae9c31c3c
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

# ╔═╡ a49251d3-c329-4491-91c2-717f54e85bd8
md"""
## Models
### Gauthier's semiforward model
Complex III and the Q cycle from Gauthier, 2013 (adapted from Demin, 2001)
"""

# ╔═╡ 9a284368-7063-4774-a2f1-f3ff892ad6c2
function c3_gauthier(;
    dpsi=150mV,
    MT_PROT=1,
    O2=6μM,
    sox_m=0.001μM,
    h_i=exp10(-7) * Molar,
    h_m=exp10(-7.6) * Molar,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    UQ = 3600μM,
    UQH2 = 400μM,
    cytc_ox = 208μM,
    cytc_rd = 325μM - cytc_ox,
    name = :c3_gauthier
    )
	
	@independent_variables t
	D = Differential(t)

    @parameters begin
        rhoC3 = 325μM    ## Complex III activity
        Q_T = 4mM        ## Total CoQ pool
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
        Emcytc = +265mV
        EmO2 = -160mV
        K03_C3 = 1666.63Hz / mM
        KEQ3_C3 = exp(iVT * (EmFeS - EmSQp_QH2p)) ## -10mV
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmbL_bHo - EmQp_SQp)) ## +130mV
        KEQ4_RD_C3 = exp(iVT * (EmbL_bHr - EmQp_SQp)) ## +70mV
        KD_Q = 22000Hz
        K06_C3 = 166.67Hz
        KEQ6_C3 = exp(iVT * (EmbH_bLo - EmbL_bHo)) ## +60mV
        K07_OX_C3 = 13.33Hz / mM
        K07_RD_C3 = 1.67Hz / mM
        KEQ7_OX_C3 = exp(iVT * (EmQn_SQn - EmbH_bLo)) ## +30mV
        KEQ7_RD_C3 = exp(iVT * (EmQn_SQn - EmbH_bLr)) ## +90mV
        K08_OX_C3 = 83.33Hz / mM
        K08_RD_C3 = 8.33Hz / mM
        KEQ8_OX_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLo)) ## +130mV
        KEQ8_RD_C3 = 9.4546   ## +60mV??? should be +190mV?
        K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))  ## -35mV
        K010_C3 = 28.33Hz / mM
        KEQ10_C3 = exp(iVT * (EmO2 - EmQp_SQp)) ## +10mV
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1)) ## +20mV
    end

    ## complex III inhibition by DOX and antimycin
    C3_CONC = rhoC3 * MT_PROT

    @variables begin
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
        fracbLrd(t)
        fracbHrd(t)
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
	## v2: QH2 diffusion (n-side -> p-side)
	v2 = KD_Q * (QH2_n - QH2_p)
    ## v3: QH2 + FeS = Q- + FeS- + 2H+
    Qo_avail = (C3_CONC - SQp) / C3_CONC * (1 - MYXOTHIAZOLE_BLOCK)
    v3 = K03_C3 * (KEQ3_C3 * Qo_avail * fes_ox * QH2_p - fes_rd * SQp * fHi^2)
    ## v4: Q- + bL = Qp + bL-
    el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
    er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
    v4ox = K04_C3 * (KEQ4_OX_C3 * SQp * el4 * cytb_1 - Q_p * er4 * cytb_2)
    v4rd = K04_C3 * (KEQ4_RD_C3 * SQp * el4 * cytb_3 - Q_p * er4 * cytb_4)
    ## v5: Q diffusion (p-side -> n-side)
    v5 = KD_Q * (Q_p - Q_n)
    ## v6: ET from bL to bH
    v6 = K06_C3 * (KEQ6_C3 * cytb_2 * exp(-iVT * β_C3 * δ₂_C3 * dpsi) - cytb_3 * exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi))
    ## v7: ET from bH to Qn
	## v8: ET from bH to SQn
    Qi_avail = (C3_CONC - SQn) / C3_CONC * (1 - ANTIMYCIN_BLOCK)
    el7 = exp(-iVT * γ_C3 * δ₃_C3 * dpsi)
    er7 = exp(iVT * γ_C3 * (1 - δ₃_C3) * dpsi)
    qn = Q_n * Qi_avail
    qh2n = QH2_n * Qi_avail
    v7ox = K07_OX_C3 *  (KEQ7_OX_C3 * cytb_3 * qn * el7 - cytb_1 * SQn * er7)
    v7rd = K07_RD_C3 * (KEQ7_RD_C3 * cytb_4 * qn * el7 - cytb_2 * SQn * er7)
    v8ox = K08_OX_C3 * (KEQ8_OX_C3 * cytb_3 * SQn * fHm^2 * el7 - cytb_1 * qh2n * er7)
    v8rd = K08_RD_C3 * (KEQ8_RD_C3 * cytb_4 * SQn * fHm^2 * el7 - cytb_2 * qh2n * er7)
    ## v9: ET from fes to cytc1
    v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
    ## v10: SQp + O2 = O2- + Q(p)
    v10 = K010_C3 * (KEQ10_C3 * O2 * SQp - sox_m * Q_p)
    ## v33: ET from cytc1 to cytc
    v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

    eqs = [
        C3_CONC ~ cytb_1 + cytb_2 + cytb_3 + cytb_4,
        C3_CONC ~ fes_ox + fes_rd,
        C3_CONC ~ cytc1_ox + cytc1_rd,
        Q_n ~ 0.5 * UQ,
        Q_p ~ 0.5 * UQ,
        QH2_n ~ 0.5 * UQH2,
        QH2_p ~ 0.5 * UQH2,
        fracbLrd ~ (cytb_2 + cytb_4) / C3_CONC,
        fracbHrd ~ (cytb_3 + cytb_4) / C3_CONC,
        ## D(UQH2) ~ dQH2n + dQH2p,
        D(SQn) ~ v7ox + v7rd - v8ox - v8rd,
        D(SQp) ~ v3 - v10 - v4ox - v4rd,
        D(cytb_1) ~ v7ox + v8ox - v4ox,
        D(cytb_2) ~ v4ox + v7rd + v8rd - v6,
        D(cytb_3) ~ v6 - v4rd - v7ox - v8ox,
        ## D(cytb_4) = v4rd - v7rd - v8rd
        D(fes_ox) ~ v9 - v3,
        D(cytc1_ox) ~ v33 - v9,
        vHresC3 ~ v6,
        vROSC3 ~ v10,
    ]
    return System(eqs, t; name)
end

# ╔═╡ e0ca9c48-d9c8-4551-bb67-2bce96e5cf52
md"""
### Semireverse model

Semireverse bc1 complex model adapted from Gauthier, 2013

- Lumped v3 and v4
- SQp is produced by reduction of Q by reduced bL at the Qo site
"""

# ╔═╡ 888d06d6-2203-491f-884d-d0e9f031aa9c
function c3_semireverse(;
    dpsi=150mV,
    MT_PROT=1,
    O2=6μM,
    sox_m=0.001μM,
    h_i=exp10(-7) * Molar,
    h_m=exp10(-7.6) * Molar,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    UQ = 3600μM,
    UQH2 = 400μM,
    cytc_ox = 208μM,
    cytc_rd = 325μM - cytc_ox,
    name = :c3_semireverse)

	@independent_variables t
	D = Differential(t)

    @parameters begin
        rhoC3 = 325μM    ## Complex III activity
        Q_T = 4mM        ## Total CoQ pool
        EmQ_C3 = +65mV   ## Ubiquinone redox potential at complex III Qo
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
		Emcytc = +265mV
        EmO2 = -160mV
        ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmFeS + EmbL_bHo - 2EmQ_C3))
        KEQ4_RD_C3 = exp(iVT * (EmFeS + EmbL_bHr - 2EmQ_C3))
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
        K011_C3 = 1000Hz / mM
        KEQ11_C3 = exp(iVT * (EmO2 - EmQp_SQp))
        ## c1_2+ + c_3+ = c1_3+ + c_2+
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1)) ## +10~20mV
    end

    # QSSA for SQp (SQp proportion is very small)
    # vROS = Qp * ((k11 * O2 * k10 * bL-) - (km11 * sox * km10 * bL)) / (km10 * bL + k11 * O2)
    ## complex III inhibition by DOX and antimycin
    C3_INHIB = 1 - ANTIMYCIN_BLOCK
    C3_CONC = rhoC3 * MT_PROT

    @variables begin
        Q_n(t)
        QH2_n(t)
        QH2_p(t)
        Q_p(t)
        SQn(t) = 0
        SQp(t)
        fes_ox(t) = C3_CONC
        fes_rd(t) ## Conserved
        cytc1_ox(t) = C3_CONC
        cytc1_rd(t) ## Conserved
        blo_bho(t) = C3_CONC
        blr_bho(t) = 0
        blo_bhr(t) = 0
        blr_bhr(t) ## Conserved
        fracbLrd(t)
        fracbHrd(t)
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

    ## QH2 + FeS + bL = Q + FeS- + bL- + 2Ho+
    ## Lumped v3 and v4
    FeS = fes_ox / C3_CONC * (1 - MYXOTHIAZOLE_BLOCK)
    FeSm = fes_rd / C3_CONC * (1 - MYXOTHIAZOLE_BLOCK)
    el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
    er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
    k4ox = K04_C3 * KEQ4_OX_C3 * el4
    k4rd = K04_C3 * KEQ4_RD_C3 * el4
    km4 = K04_C3 * er4 * fHi^2
    v4ox = k4ox * QH2_p * FeS * blo_bho - km4 * Q_p * FeSm * blr_bho
    v4rd = k4rd * QH2_p * FeS * blo_bhr - km4 * Q_p * FeSm * blr_bhr

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
    k11 = K011_C3 * KEQ11_C3
    km11 = K011_C3
    v11 = k11 * SQp * O2 - km11 * Q_p * sox_m

	sqp_ratio = (k10ox * blr_bho + k10rd * blr_bhr + km11 * sox_m) / (km10 * (blo_bho + blo_bhr) + k11 * O2)

    eqs = [
        C3_CONC ~ blo_bho + blr_bho + blo_bhr + blr_bhr,
        C3_CONC ~ fes_ox + fes_rd,
        C3_CONC ~ cytc1_ox + cytc1_rd,
        Q_n ~ 0.5 * UQ,
        Q_p ~ 0.5 * UQ,
        QH2_n ~ 0.5 * UQH2,
        QH2_p ~ 0.5 * UQH2,
        SQp ~ Q_p * sqp_ratio,
        fracbLrd ~ (blr_bho + blr_bhr) / C3_CONC,
        fracbHrd ~ (blo_bhr + blr_bhr) / C3_CONC,
        ## D(UQH2) ~ dQH2n + dQH2p,
        D(SQn) ~ v7ox + v7rd - v8ox - v8rd,
        # D(SQp) ~ v10ox + v10rd - v11, ## =0
        D(blo_bho) ~ v7ox + v8ox - v4ox + v10ox,
        D(blr_bho) ~ v4ox + v7rd + v8rd - v6 - v10ox,
        D(blo_bhr) ~ v6 - v4rd - v7ox - v8ox + v10rd,
        ## D(blr_bhr) = v4rd - v7rd - v8rd - v10rd
        D(fes_ox) ~ v9 - v4ox - v4rd,
        D(cytc1_ox) ~ v33 - v9,
        vHresC3 ~ v6,
        vROSC3 ~ v11,
    ]
    return System(eqs, t; name)
end

# ╔═╡ 5349e449-a9bd-42ce-8114-e1be9500c33e
md"""
### Repulsion model

- Reduced heme bL (bL-) blocks QH2 oxidation at Qo site due to repulsion between bL- and Q-
- SQp less unstable than that in the Gauthier model (Em(Q/SQ) = -300mV)
"""

# ╔═╡ 6fa3d1fb-cf42-4b4d-83fa-724e7f41aae1
function c3_repulsion(;
    dpsi=150mV,
    MT_PROT=1,
    O2=6μM,
    sox_m=0.001μM,
    h_i=exp10(-7) * Molar,
    h_m=exp10(-7.6) * Molar,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    UQ = 3600μM,
    UQH2 = 400μM,
    cytc_ox = 208μM,
    cytc_rd = 325μM - cytc_ox,
    name = :c3_repulse
    )
	@independent_variables t
	D = Differential(t)

    @parameters begin
        rhoC3 = 325μM    ## Complex III activity
        Q_T = 4mM        ## Total CoQ pool
        EmQp = +60mV
        EmSQp_QH2p = +400mV
        EmQp_SQp = 2EmQp - EmSQp_QH2p
        EmQn_SQn = +50mV
        EmSQn_QH2n = +150mV
        EmbL_bHo = -40mV
        EmbL_bHr = EmbL_bHo - 60mV
        EmbH_bLo = +20mV
        EmbH_bLr = EmbH_bLo - 60mV
        EmFeS = +280mV
        Emcytc1 = +245mV
        Emcytc = +265mV
		EmO2 = -160mV
        K03_C3 = 2E5Hz / mM ## 1666.63Hz / mM
        KEQ3_C3 = exp(iVT * (EmFeS - EmSQp_QH2p)) ## ~ 0.01
        K04_C3 = 50.67Hz / mM
        KEQ4_OX_C3 = exp(iVT * (EmbL_bHo - EmQp_SQp))
        KEQ4_RD_C3 = exp(iVT * (EmbL_bHr - EmQp_SQp))
        KD_Q = 22000Hz
        K06_C3 = 10000Hz ## 166.67Hz
        KEQ6_C3 = exp(iVT * (EmbH_bLo - EmbL_bHo))
        K07_OX_C3 = 13.33Hz / mM
        K07_RD_C3 = 1.67Hz / mM
        KEQ7_OX_C3 = exp(iVT * (EmQn_SQn - EmbH_bLo))
        KEQ7_RD_C3 = exp(iVT * (EmQn_SQn - EmbH_bLr))
        K08_OX_C3 = 83.33Hz / mM
        K08_RD_C3 = 8.33Hz / mM
        KEQ8_OX_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLo))
        KEQ8_RD_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLr))
        K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))
        K010_C3 = 28.33Hz / mM
        KEQ10_C3 = exp(iVT * (EmO2 - EmQp_SQp))
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1))
    end

    ## complex III inhibition by DOX and antimycin
    C3_CONC = rhoC3 * MT_PROT

    @variables begin
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
        blo_bho(t) = C3_CONC
        blo_bhr(t) = 0
        blr_bho(t) = 0
        blr_bhr(t) ## Conserved
        fracbLrd(t)
        fracbHrd(t)
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
    ## QH2 + FeS = Q- + FeS- + 2H+
    Qo_avail = (C3_CONC - SQp) / C3_CONC * (1 - fracbLrd) * (1 - MYXOTHIAZOLE_BLOCK)
    v3 = K03_C3 * (KEQ3_C3 * Qo_avail * fes_ox * QH2_p - fes_rd * SQp * fHi^2)
    ## Q- + bL = Qp + bL-
    el4 = exp(-iVT * α_C3 * δ₁_C3 * dpsi)
    er4 = exp(iVT * α_C3 * (1 - δ₁_C3) * dpsi)
    v4ox = K04_C3 * (KEQ4_OX_C3 * SQp * el4 * blo_bho - Q_p * er4 * blr_bho)
    v4rd = K04_C3 * (KEQ4_RD_C3 * SQp * el4 * blo_bhr - Q_p * er4 * blr_bhr)
    ## v5 = Q diffusion (p-side -> n-side)
    v5 = KD_Q * (Q_p - Q_n)
    ## v6 = bL to bH
    el6 = exp(-iVT * β_C3 * δ₂_C3 * dpsi)
    er6 = exp(iVT * β_C3 * (1 - δ₂_C3) * dpsi)
    v6 = K06_C3 * (KEQ6_C3 * blr_bho * el6 - blo_bhr * er6)
    ## v7 = bH to Qn; v8: bH to SQn
    Qi_avail = (C3_CONC - SQn) / C3_CONC * (1 - ANTIMYCIN_BLOCK)
    el7 = exp(-iVT * γ_C3 * δ₃_C3 * dpsi)
    er7 = exp(iVT * γ_C3 * (1 - δ₃_C3) * dpsi)
    qn = Q_n * Qi_avail
    qh2n = QH2_n * Qi_avail
    v7ox = K07_OX_C3 * (KEQ7_OX_C3 * blo_bhr * qn * el7 - blo_bho * SQn * er7)
    v7rd = K07_RD_C3 * (KEQ7_RD_C3 * blr_bhr * qn * el7 - blr_bho * SQn * er7)
    v8ox = K08_OX_C3 * (KEQ8_OX_C3 * blo_bhr * SQn * fHm^2 * el7 - blo_bho * qh2n * er7)
    v8rd = K08_RD_C3 * (KEQ8_RD_C3 * blr_bhr * SQn * fHm^2 * el7 - blr_bho * qh2n * er7)
    ## v9 = fes -> cytc1
    v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
    ## SQp + O2 -> O2- + Q
    v10 = K010_C3 * (KEQ10_C3 * O2 * SQp - sox_m * Q_p)
    ## cytc1_2+  + cytc_3+ = cytc1_3+  + cytc_2+
    v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

    eqs = [
        C3_CONC ~ blo_bho + blr_bho + blo_bhr + blr_bhr,
        C3_CONC ~ fes_ox + fes_rd,
        C3_CONC ~ cytc1_ox + cytc1_rd,
        Q_n ~ 0.5 * UQ,
        Q_p ~ 0.5 * UQ,
        QH2_n ~ 0.5 * UQH2,
        QH2_p ~ 0.5 * UQH2,
        fracbLrd ~ (blr_bho + blr_bhr) / C3_CONC,
        fracbHrd ~ (blo_bhr + blr_bhr) / C3_CONC,
        ## D(UQH2) ~ dQH2n + dQH2p,
        D(SQn) ~ v7ox + v7rd - v8ox - v8rd,
        D(SQp) ~ v3 - v10 - v4ox - v4rd,
        D(blo_bho) ~ v7ox + v8ox - v4ox,
        D(blr_bho) ~ v4ox + v7rd + v8rd - v6,
        D(blo_bhr) ~ v6 - v4rd - v7ox - v8ox,
        ## D(blr_bhr) = v4rd - v7rd - v8rd
        D(fes_ox) ~ v9 - v3,
        D(cytc1_ox) ~ v33 - v9,
        vHresC3 ~ v6,
        vROSC3 ~ v10,
    ]
    return System(eqs, t; name)
end

# ╔═╡ a486b939-b40a-414f-b065-728063523748
md"""
### Rapid equlibrium model

Assuming electron transfer between Qo-bL-bH-Qi and quinone binding/unbinding are fast and there are the following reactions

1. ``\ce{QH2 + CytC^3+ + b_L = Q^{-}b_L + CytC^2+ + 2H^+}``
2. ``\ce{b_H^-Q_n^- + 2H^+ = b_H + QH2}``
3. ``\ce{Q^{-}b_L^* + O2 = Q + b_L^* + O2^-}``
"""

# ╔═╡ 0c97fae0-034c-4d61-b3ea-90bcef259715
function c3_equlibrium(;
    dpsi=150mV,
    MT_PROT=1,
    O2=6μM,
    sox_m=0.001μM,
    h_i=exp10(-7) * Molar,
    h_m=exp10(-7.6) * Molar,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOLE_BLOCK=0,
    UQ = 3600μM,
    UQH2 = 400μM,
    cytc_ox = 208μM,
    cytc_rd = 325μM - cytc_ox,
    name = :c3_equlibrium
    )

	@independent_variables t
	D = Differential(t)

	@parameters begin
		rhoC3 = 325μM    ## Complex III activity
        Q_T = 4mM        ## Total CoQ pool
		## Association constant of Q/QH2 at the Qo site
		KA_Qo = inv(1mM) 
		## Association constant of Q/QH2 at the Qi site
		KA_Qi = inv(400μM)
		## Redox potentials
		EmQp = +60mV
		EmSQp_QH2p = +290mV
		EmQp_SQp = 2EmQp - EmSQp_QH2p ## -170mV
		EmQn_SQn = +50mV
		EmSQn_QH2n = +150mV
		EmbL_bHo = -40mV
        EmbL_bHr = EmbL_bHo - 60mV
        EmbH_bLo = +20mV
        EmbH_bLr = EmbH_bLo - 60mV
		EmFeS = +280mV
        Emcytc1 = +245mV
        Emcytc = +265mV
        EmO2 = -160mV
		## Precalulated fractions from redox potentials
		rQp_SQp = exp(iVT * EmQp_SQp) 
		rQn_SQn = exp(iVT * EmQn_SQn) 
		rbL_bHo = exp(iVT * EmbL_bHo)
        rbL_bHr = exp(iVT * EmbL_bHr)
        rbH_bLo = exp(iVT * EmbH_bLo)
        rbH_bLr = exp(iVT * EmbH_bLr)
		## v3: QH2 + FeS + bL = SQbL + FeS- + 2H+
		K03_C3 = 60Hz / mM
        KEQ3_C3 = exp(iVT * (EmFeS - EmSQp_QH2p)) # -10mV
		## v8: bH-Q- + 2H+ = bH + QH2
        K08_OX_C3 = 83.33Hz / mM
        K08_RD_C3 = 8.33Hz / mM
        KEQ8_OX_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLo)) / KA_Qi # +130mV
        KEQ8_RD_C3 = exp(iVT * (EmSQn_QH2n - EmbH_bLr)) / KA_Qi # +190mV
		K09_C3 = 832.48Hz / mM
        KEQ9_C3 = exp(iVT * (Emcytc1 - EmFeS))  # -35mV
        K010_C3 = 28.33Hz / mM
        KEQ10_C3 = exp(iVT * (EmO2 - EmQp_SQp)) # +10mV
        K33_C3 = 2469.13Hz / mM
        KEQ33_C3 = exp(iVT * (Emcytc - Emcytc1)) # +20mV
	end
	
	@variables begin
		Q_n(t)
        QH2_n(t)
        QH2_p(t)
        Q_p(t)
		## States use the format (Qo, bL, bH, Qi)
		## 0-electron states
		C3_0(t)
		C3_0000(t)
		## 1-electron states
		C3_1(t)
		C3_1000(t)
		C3_0100(t)
		C3_0010(t)
		C3_0001(t)
		## 2-electron states
		C3_2(t)
		C3_1100(t)
		C3_1010(t)
		C3_1001(t)
		C3_0110(t)
		C3_0101(t)
		C3_0011(t)
		## 3-electron states
		C3_3(t)
		C3_1110(t)
		C3_1101(t)
		C3_1011(t)
		C3_0111(t)
		## 4-electron states (?)
		C3_4(t)
		C3_1111(t)
		## Bound ubiquinone in C3
		SQp(t)
		SQn(t)
		fdpsiC3(t) ## Unit of MMP effect
        fes_ox(t)
        fes_rd(t) ## Conserved
        cytc1_ox(t)
        cytc1_rd(t) ## Conserved
        fracbLrd(t)
        fracbHrd(t)
		vQpC3(t)
		vQH2pC3(t)
		vQnC3(t)
		vQH2nC3(t)
        vROSC3(t)
        vHresC3(t)
	end

	## Complex 3 content
	C3_CONC = rhoC3 * MT_PROT

	## Normalized proton concentration
    fHm = h_m * inv(1E-7Molar)
    fHi = h_i * inv(1E-7Molar)
	fQo = Q_p * KA_Qo
	fQi = Q_n * KA_Qi
	## State occupancy weights
	## (Qo, bL, bH, Qi) = (0/1, 0/1, 0/1, 0/1)
	## The dielectric distance is (0.25, 0.5, 0.25)
	w_0000 = 1
	w_1000 = fQo * rQp_SQp * fdpsiC3^2
	w_0100 = rbL_bHo * fdpsiC3
	w_0010 = rbH_bLo / fdpsiC3
	w_0001 = fQi * rQn_SQn / fdpsiC3^2
	w_1100 = w_1000 * w_0100
	w_1010 = w_1000 * w_0010
	w_1001 = w_1000 * w_0001
	w_0110 = rbL_bHr * fdpsiC3 * rbH_bLr / fdpsiC3
	w_0101 = w_0100 * w_0001
	w_0011 = w_0010 * w_0001
	w_1110 = w_1000 * w_0110
	w_1101 = w_1100 * w_0001
	w_1011 = w_1000 * w_0011
	w_0111 = w_0110 * w_0001

	den1 = w_1000 + w_0100 + w_0010 + w_0001
	den2 = w_1100 + w_1010 + w_1001 + w_0110 + w_0101 + w_0011
	den3 = w_1110 + w_1101 + w_1101 + w_0111

	## QH2 + 00xx + FeS = 10xx + FeS- + 2H+
	## Only oxidized bL are eligible for the reaction
	k3 = K03_C3 * KEQ3_C3 * QH2_p * KA_Qo * fes_ox * (1 - MYXOTHIAZOLE_BLOCK)
	km3 = K03_C3 * fes_rd * fHi^2 * (1 - MYXOTHIAZOLE_BLOCK)
	v01 = k3 * C3_0000 - km3 * C3_1000
	v12 = k3 * (C3_0010 + C3_0001) - km3 * (C3_1010 + C3_1001)
	v23 = k3 * C3_0011 - km3 * C3_1011

	## xx11 + 2H+ = xx00 + QH2
	Qi_avail = (C3_CONC - SQn) / C3_CONC * (1 - ANTIMYCIN_BLOCK)
    el7 = exp(-iVT * 0.25 * 0.5 * dpsi)
    er7 = exp(iVT * 0.25 * (1 - 0.5) * dpsi)
	k8ox = K08_OX_C3 * KEQ8_OX_C3 * el7 * fHm^2
	km8ox = K08_OX_C3 * er7 * Qi_avail * QH2_n
	k8rd = K08_RD_C3 * KEQ8_RD_C3 * el7 * fHm^2
	km8rd = K08_RD_C3 * er7 * Qi_avail * QH2_n

	v20 = k8ox * C3_0011 - km8ox * C3_0000
	v31ox = k8ox * C3_1011 - km8ox * C3_1000
	v31rd = k8rd * C3_0111 - km8rd * C3_0100

	## 1xxx + O2 = 0xxx + Q + O2-
	k10 = K010_C3 * KEQ10_C3 * O2
	km10 = K010_C3 * Q_p * KA_Qo * sox_m
	v32 = k10 * (C3_1110 + C3_1101 + C3_1011) - km10 * (C3_0110 + C3_0101 + C3_0011)
	v21 = k10 * (C3_1100 + C3_1010 + C3_1001) - km10 * (C3_0100 + C3_0010 + C3_0001)
	v10 = k10 * C3_1000 - km10 * C3_0000

	## v9: ET from fes to cytc1
    v9 = K09_C3 * (KEQ9_C3 * fes_rd * cytc1_ox - fes_ox * cytc1_rd)
    ## v33: ET from cytc1 to cytc
    v33 = K33_C3 * (KEQ33_C3 * cytc1_rd * cytc_ox - cytc1_ox * cytc_rd)

	eqs = [
		## D(C3_0) ~ -v01 + v20,
		D(C3_1) ~ v01 - v12 + v31ox + v31rd,
		D(C3_2) ~ v12 - v23 - v20,
		D(C3_3) ~ v23 - v31ox - v31rd,
		D(fes_ox) ~ v9 - (v01 + v12 + v23),
        D(cytc1_ox) ~ v33 - v9,
		SQn ~ C3_0001 + C3_1001 + C3_0101 + C3_0011 + C3_1101 + C3_1011 + C3_0111,
		SQp ~ C3_1000 + C3_1100 + C3_1010 + C3_1001 + C3_1110 + C3_1101 + C3_1011,
		C3_CONC ~ C3_0 + C3_1 + C3_2 + C3_3,
		C3_CONC ~ fes_ox + fes_rd,
		C3_CONC ~ cytc1_ox + cytc1_rd,
		Q_n ~ 0.5 * UQ,
		Q_p ~ 0.5 * UQ,
		QH2_n ~ 0.5 * UQH2,
		QH2_p ~ 0.5 * UQH2,
		fdpsiC3 ~ exp(iVT * dpsi / 4),
		C3_0000 ~ C3_0,
		C3_1000 ~ C3_1 * w_1000 / den1,
		C3_0100 ~ C3_1 * w_0100 / den1,
		C3_0010 ~ C3_1 * w_0010 / den1,
		C3_0001 ~ C3_1 * w_0001 / den1,
		C3_1100 ~ C3_2 * w_1100 / den2,
		C3_1010 ~ C3_2 * w_1010 / den2,
		C3_1001 ~ C3_2 * w_1001 / den2,
		C3_0110 ~ C3_2 * w_0110 / den2,
		C3_0101 ~ C3_2 * w_0101 / den2,
		C3_0011 ~ C3_2 * w_0011 / den2,
		C3_1110 ~ C3_3 * w_1110 / den3,
		C3_1101 ~ C3_3 * w_1101 / den3,
		C3_1011 ~ C3_3 * w_1011 / den3,
		C3_0111 ~ C3_3 * w_0111 / den3,
		fracbLrd ~ (C3_0100 + C3_1100 + C3_0110 + C3_0101 + C3_1110 + C3_1101 + C3_0111) / C3_CONC,
		fracbHrd ~ (C3_0010 + C3_1010 + C3_0110 + C3_0011 + C3_1110 + C3_1011 + C3_0111) / C3_CONC,
		## TODO: derive Q binding/unbinding rates
		vQpC3 ~ - (D(C3_1) * (w_1000) / den1 + D(C3_2) * (w_1100 + w_1010 + w_1001) / den2 + D(C3_3) * (w_1110 + w_1101 + w_1011) / den3),
		vQH2pC3 ~ -v01 - v12 - v23,
		vQnC3 ~ - (D(C3_1) * (w_0001) / den1 + D(C3_2) * (w_1001 + w_0101 + w_0011) / den2 + D(C3_3) * (w_1101 + w_1011 + w_0111) / den3),
		vQH2nC3 ~ v20 + v31ox + v31rd,	
		vROSC3 ~ v32 + v21 + v10,
        vHresC3 ~ 2 * vQH2nC3,
	]
	
	return System(eqs, t; name, defaults=[
		C3_1 => 0,
		C3_2 => 0,
		C3_3 => 0, 
		fes_ox=>C3_CONC, 
		cytc1_ox=>C3_CONC
	])
end

# ╔═╡ 621e719a-719c-47c8-b2a5-1debd82bc470
md"""
### Setup the ODE problem
"""

# ╔═╡ 2345cb6e-668b-4f36-82fb-a9fcb72775c7
@parameters begin
    UQ = 3600μM
    UQH2 = 400μM
    dpsi = 150mV
    cytc_ox = 208μM
    cytc_rd = 325μM - cytc_ox
    sox_m = 0.01μM
	O2 = 6μM
end

# ╔═╡ 5eecad75-cb22-4dea-ab84-0bc670d402b8
gsys = c3_gauthier(;dpsi, cytc_ox, cytc_rd, UQ, UQH2, sox_m, O2) |> mtkcompile

# ╔═╡ 21ecac6c-f4a2-4754-b671-fad8f7b88788
prob_g = SteadyStateProblem(gsys, [])

# ╔═╡ 8de17909-61e3-40a3-91d8-7af48a26a1bd
esys = c3_equlibrium(; dpsi, cytc_ox, cytc_rd, UQ, UQH2, sox_m, O2) |> mtkcompile

# ╔═╡ da9d74b7-2b41-4ab0-9cc7-d822c6d12753
observed(esys)

# ╔═╡ fe78a3d5-6d9d-4b51-bf07-c899a6536f66
rsys = c3_repulsion(;dpsi, cytc_ox, cytc_rd, UQ, UQH2, sox_m, O2) |> mtkcompile

# ╔═╡ 9c74298a-80ae-43bb-96cd-d9e86129cb7f
prob_r = SteadyStateProblem(rsys, [rsys.K010_C3 => 33Hz / mM])

# ╔═╡ 8e76757c-6bd4-430b-a15b-430db05f5271
ssys = c3_semireverse(;dpsi, cytc_ox, cytc_rd, UQ, UQH2, sox_m, O2) |> mtkcompile

# ╔═╡ 99a4b73f-a0b2-4125-8da6-daff6f8ad80f
prob_s = SteadyStateProblem(ssys, [ssys.K010_C3 => 130Hz / mM, ssys.K011_C3 => 10000Hz / mM, ssys.K04_C3 => 50Hz / mM])

# ╔═╡ 118cf64a-dcf6-4f3f-9b3a-bd3ab9f88d66
alg = DynamicSS(KenCarp47())

# ╔═╡ 31f75bf3-c8d0-48eb-9091-5a7e4ea0c994
ealg = EnsembleThreads()

# ╔═╡ 44257bde-ba49-48d7-b025-d8ed2a2568e6
# Utility function
extract(sim, k) = map(s -> s[k], sim)

# ╔═╡ 8d5d8eef-f210-4c67-af6e-82e3ef29771a
md"""
## Varying MMP
"""

# ╔═╡ df83f7f4-3ffd-4285-8ef2-94fad068855b
dpsirange = 100:1:200

# ╔═╡ 8b5a16cb-a40f-4115-a476-9474890e053e
alter_dpsi = (prob, i, repeat) -> begin
    prob.ps[dpsi] = dpsirange[i] * mV
    prob
end

# ╔═╡ 64fb2d71-2613-4a9a-b34d-c66f00289e54
eprob_g = EnsembleProblem(prob_g; prob_func=alter_dpsi)

# ╔═╡ d50dfab5-1bb3-49fd-bff7-e3721b42d93c
@time sim_g = solve(eprob_g, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 21b2880b-2724-4d62-9ea4-1ab77797525d
eprob_r = EnsembleProblem(prob_r; prob_func=alter_dpsi)

# ╔═╡ 82ee477e-da1b-4732-abbf-1ec32739fd80
@time sim_r = solve(eprob_r, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 53c8dcfc-bfe1-46e3-871d-50a0af14ebf3
eprob_s = EnsembleProblem(prob_s; prob_func=alter_dpsi)

# ╔═╡ f447d5c4-c9e1-4c19-bd1c-5218ba9810e7
@time sim_s = solve(eprob_s, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 0cd9400e-b458-41a5-86db-ab0687c0b973
prob_e = SteadyStateProblem(esys, [
	esys.K03_C3 => 3900Hz / mM,
	esys.EmQp => 65mV,
	esys.EmSQp_QH2p => +290mV,
	esys.K08_OX_C3 => 83.33Hz / mM,  ## 83.33
    esys.K08_RD_C3 => 8.33Hz / mM,   ## 8.33
	esys.EmbH_bLo => +20mV,
	esys.K010_C3 => 250Hz / mM,
])

# ╔═╡ 5d4ec8b9-9aa7-45e2-a3c5-78d350ab34af
eprob_e = EnsembleProblem(prob_e; prob_func=alter_dpsi)

# ╔═╡ ef63cdf5-0ae6-490f-8f67-2aa64adbeabe
@time sim_e = solve(eprob_e, alg, ealg; trajectories=length(dpsirange), abstol=1e-8, reltol=1e-8)

# ╔═╡ 947c932f-d3aa-4c81-840b-85b73e138980
let
	xs = dpsirange
	ys = [extract(sim_g, gsys.vHresC3) extract(sim_s, ssys.vHresC3) extract(sim_r,rsys.vHresC3) extract(sim_e, esys.vHresC3)]
	plot(xs, ys, xlabel="MMP (mV)", ylabel="Resp. Rate (mM/s)", label=["Semiforward" "Semireverse" "Repulsion" "Eqilibrium"])
end |> PNG

# ╔═╡ 3d48bba5-6f3e-4e47-b2ce-5d7a98dae55d
let
	xs = dpsirange
	ys = [extract(sim_s, ssys.fracbLrd) extract(sim_e, esys.fracbLrd) extract(sim_s, ssys.fracbHrd) extract(sim_e, esys.fracbHrd)]
	plot(xs, ys, xlabel="MMP (mV)", ylabel="Reduced fraction", label=["Semireverse (bL)" "E (bL)" "Semireverse (bH)" "E (bH)"], line=[:solid :dash :solid :dash])	
end |> PNG

# ╔═╡ 8ec0a9f7-f94b-4eb0-a5db-bca292531fa5
let
	xs = dpsirange
	ys = [extract(sim_g, gsys.vROSC3) extract(sim_s, ssys.vROSC3) extract(sim_r, rsys.vROSC3) extract(sim_e, esys.vROSC3)]
	plot(xs, ys, xlabel="MMP (mV)", ylabel="ROS Rate (mM/s)", label=["G" "S" "R" "E"])
end |> PNG

# ╔═╡ 06244f1e-049e-4486-903d-f0bcdbb94dcb
md"""
## Varying UQH2
"""

# ╔═╡ a0a5c6fb-3591-4b3f-b575-5deec1893c2e
qh2range = 10:10:3990

# ╔═╡ 7bdd4656-d836-4cd4-8471-e0755553bb11
alter_qh2 = (prob, i, repeat) -> begin
    prob.ps[UQH2] = qh2range[i] * μM
    prob.ps[UQ] = 4000μM - prob.ps[UQH2]
    prob
end

# ╔═╡ 6ffcafd2-5647-4444-83c3-f8d1b966cc45
@time qh2_g = solve(EnsembleProblem(prob_g; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 439454e4-32e2-4b59-b62a-2bb591ea9ea2
@time qh2_r = solve(EnsembleProblem(prob_r; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 2979d921-a0c0-42c0-a1f0-27ffe9793fb8
@time qh2_s = solve(EnsembleProblem(prob_s; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 9ca0dc68-d53a-4c14-9e0e-1eabe0f7681a
@time qh2_e = solve(EnsembleProblem(prob_e; prob_func=alter_qh2), alg, ealg; trajectories=length(qh2range), abstol=1e-8, reltol=1e-8)

# ╔═╡ 337076c9-3e11-4641-b3b7-d35f18ccb2b0
let
	xs = qh2range ./ 4000 .* 100
	ys = [extract(qh2_g, gsys.vHresC3) extract(qh2_s, ssys.vHresC3) extract(qh2_r, rsys.vHresC3) extract(qh2_e, esys.vHresC3)]
	plot(xs, ys, xlabel="QH2 (%)", ylabel="Resp. Rate (mM/s)", label=["G" "S" "R" "E"])
end |> PNG

# ╔═╡ 24422da1-2f02-4133-abbb-8e693dc4a4f5
let 
	xs = qh2range ./ 4000 .* 100
	ys = [extract(qh2_g, gsys.vROSC3) extract(qh2_s, ssys.vROSC3) extract(qh2_r, rsys.vROSC3) extract(qh2_e, esys.vROSC3)]
	plot(xs, ys, xlabel="QH2 (%)", ylabel="ROS Rate (mM/s)", label=["G" "S" "R" "E"])
end |> PNG

# ╔═╡ 2de8ad2d-21d2-4033-9b4d-870e4d07856f
let
	xs = qh2range ./ 4000 .* 100
	ys = [extract(qh2_s, ssys.fracbLrd) extract(qh2_e, esys.fracbLrd) extract(qh2_s, ssys.fracbHrd) extract(qh2_e, esys.fracbHrd)]
	plot(xs, ys, xlabel="QH2 (%)", ylabel="Reduced fraction", label=["Semireverse (bL)" "E (bL)" "Semireverse (bH)" "E (bH)"], line=[:solid :dash :solid :dash])	
end |> PNG

# ╔═╡ Cell order:
# ╠═13270f16-84b6-11f0-033c-a9652a6e1210
# ╠═acd4eb16-255d-46ee-b498-7e9787837492
# ╠═9061f18a-31ad-46ac-b0fa-aa632d63147c
# ╠═a6e0a7c4-db17-4041-94d3-ec0ae9c31c3c
# ╠═a49251d3-c329-4491-91c2-717f54e85bd8
# ╠═9a284368-7063-4774-a2f1-f3ff892ad6c2
# ╟─e0ca9c48-d9c8-4551-bb67-2bce96e5cf52
# ╟─888d06d6-2203-491f-884d-d0e9f031aa9c
# ╟─5349e449-a9bd-42ce-8114-e1be9500c33e
# ╟─6fa3d1fb-cf42-4b4d-83fa-724e7f41aae1
# ╠═a486b939-b40a-414f-b065-728063523748
# ╠═0c97fae0-034c-4d61-b3ea-90bcef259715
# ╠═621e719a-719c-47c8-b2a5-1debd82bc470
# ╠═2345cb6e-668b-4f36-82fb-a9fcb72775c7
# ╠═5eecad75-cb22-4dea-ab84-0bc670d402b8
# ╠═21ecac6c-f4a2-4754-b671-fad8f7b88788
# ╠═8de17909-61e3-40a3-91d8-7af48a26a1bd
# ╠═da9d74b7-2b41-4ab0-9cc7-d822c6d12753
# ╠═fe78a3d5-6d9d-4b51-bf07-c899a6536f66
# ╠═9c74298a-80ae-43bb-96cd-d9e86129cb7f
# ╠═8e76757c-6bd4-430b-a15b-430db05f5271
# ╠═99a4b73f-a0b2-4125-8da6-daff6f8ad80f
# ╠═118cf64a-dcf6-4f3f-9b3a-bd3ab9f88d66
# ╠═31f75bf3-c8d0-48eb-9091-5a7e4ea0c994
# ╠═44257bde-ba49-48d7-b025-d8ed2a2568e6
# ╠═8d5d8eef-f210-4c67-af6e-82e3ef29771a
# ╠═df83f7f4-3ffd-4285-8ef2-94fad068855b
# ╠═8b5a16cb-a40f-4115-a476-9474890e053e
# ╠═64fb2d71-2613-4a9a-b34d-c66f00289e54
# ╠═d50dfab5-1bb3-49fd-bff7-e3721b42d93c
# ╠═21b2880b-2724-4d62-9ea4-1ab77797525d
# ╠═82ee477e-da1b-4732-abbf-1ec32739fd80
# ╠═53c8dcfc-bfe1-46e3-871d-50a0af14ebf3
# ╠═f447d5c4-c9e1-4c19-bd1c-5218ba9810e7
# ╠═5d4ec8b9-9aa7-45e2-a3c5-78d350ab34af
# ╠═ef63cdf5-0ae6-490f-8f67-2aa64adbeabe
# ╠═947c932f-d3aa-4c81-840b-85b73e138980
# ╠═3d48bba5-6f3e-4e47-b2ce-5d7a98dae55d
# ╠═8ec0a9f7-f94b-4eb0-a5db-bca292531fa5
# ╠═0cd9400e-b458-41a5-86db-ab0687c0b973
# ╠═06244f1e-049e-4486-903d-f0bcdbb94dcb
# ╠═a0a5c6fb-3591-4b3f-b575-5deec1893c2e
# ╠═7bdd4656-d836-4cd4-8471-e0755553bb11
# ╠═6ffcafd2-5647-4444-83c3-f8d1b966cc45
# ╠═439454e4-32e2-4b59-b62a-2bb591ea9ea2
# ╠═2979d921-a0c0-42c0-a1f0-27ffe9793fb8
# ╠═9ca0dc68-d53a-4c14-9e0e-1eabe0f7681a
# ╠═337076c9-3e11-4641-b3b7-d35f18ccb2b0
# ╠═24422da1-2f02-4133-abbb-8e693dc4a4f5
# ╠═2de8ad2d-21d2-4033-9b4d-870e4d07856f
