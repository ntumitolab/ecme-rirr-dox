#===
Coupled LCC-RyR local control model (9 states) by Hinch et al. (2004) https://pmc.ncbi.nlm.nih.gov/articles/PMC1304886/

===#
function get_cicr9_sys(; ca_i, ca_sr, ca_o, vm, A_CAP=1.534E-4cm², V_DS=2e-4μm^3, name=:cicr9)
    @parameters begin
        N_CaRU = 50000      # Number of Ca release units (LCC/RyR pairs)
        V_DS_TOT = V_DS * N_CaRU # Total subspace volume
        P_LCC = 9.13E-13 * cm^3 * Hz # Single LCC conductance
        R_RYR = 19.6 / ms
        R_XFER = 220 / ms   # Ca2+ diffusion rate between dyadic space and cytosol
        V_L = -2mV          # Potential when half LCC open
        kV_L = inv(7mV)     # Steepness of opening potentials
        phi_L = 2.35        # Proportion of time closed in open mode
        t_L = 1ms           # Time switching between C and O states
        it_L = inv(t_L)
        tau_L = 650ms       # LCC inactivation time
        itau_L = inv(tau_L)
        K_L = 0.22μM        # Ca Concentration at LCC inactivation
        iK_L = inv(K_L)
        a_L = 0.0625        # Biasing to make inactivation function of V
        b_L = 14            # Biasing to make inactivation function of_V
        K_RyR = 41μM        # Ca Concentration at RyR activation
        t_R = 1.17t_L       # Time switching between C and O states in RyR
        it_R = inv(t_R)
        phi_R = 0.05        # Proportion of time closed in open mode in RyR
        tau_R = 2.43ms      # RyR inactivation time
        itauR = inv(tau_R)
        θ_R = 0.012         # reciprocal of proportion of time inactivated in open mode
        c_R = 0.01          # Biasing to make inactivation a function of [Ca2+]ds
        d_R = 100           # Biasing to make inactivation a function of [Ca2+]ds
    end

    @variables begin
        # LCC/RyR close/open/inactivated states
        y_oo(t) = 0
        y_oc(t) = 0
        y_oi(t) = 0
        y_co(t) = 0
        y_cc(t) # conserved
        y_ci(t) = 0
        y_io(t) = 0
        y_ic(t) = 0
        y_ii(t) = 0
        # Distribution of LCC/RyR states
        p_cc(t)
        p_co(t)
        p_oc(t)
        p_oo(t)
        # Subspace Ca
        ca_ss_cc(t)
        ca_ss_co(t)
        ca_ss_oc(t)
        ca_ss_oo(t)
        ca_ss(t)
        ICaL(t)
        JRyR(t)
        Jxfer(t)
        JLCC(t)
        α1_L(t)
        αp_L(t)
        αm_L(t)
        ϵm_L(t)
        βm_R(t)
    end

    _ϵp(ca) = ca * itau_L * iK_L * (α1_L + a_L) / (α1_L + 1)
    _βp(ca) = it_R * hil(ca, K_RyR, 2)
    _μp(ca) = itauR * (ca^2 + c_R * K_RyR^2) / (ca^2 + K_RyR^2)
    _μm(ca) = itauR * θ_R * d_R * (ca^2 + c_R * K_RyR^2) / (d_R * ca^2 + c_R * K_RyR^2)

    # State transition rates
    v_oo_oc = y_oo * βm_R - y_oc * _βp(ca_ss_oc)
    v_oc_oi = y_oc * _μp(ca_ss_oc) - y_oi * _μm(ca_ss_oc)
    v_oo_co = y_oo * αm_L - y_co * αp_L
    v_oc_cc = y_oc * αm_L - y_cc * αp_L
    v_oi_ci = y_oi * αm_L - y_ci * αp_L
    v_co_cc = y_co * βm_R - y_cc * _βp(ca_ss_cc)
    v_cc_ci = y_cc * _μp(ca_ss_cc) - y_ci * _μm(ca_ss_cc)
    v_co_io = y_co * _ϵp(ca_ss_co) - y_io * ϵm_L
    v_cc_ic = y_cc * _ϵp(ca_ss_cc) - y_ic * ϵm_L
    v_ci_ii = y_ci * _ϵp(ca_ss_cc) - y_ii * ϵm_L
    v_io_ic = y_io * βm_R - y_ic * _βp(ca_ss_cc)
    v_ic_ii = y_ic * _μp(ca_ss_cc) - y_ii * _μm(ca_ss_cc)

    # calcium flux weights
    wlcc = P_LCC / V_DS * exprel(-2iVT * vm)
    wryr = R_RYR
    wcai = R_XFER

    eqs = [
        Jxfer ~ R_XFER * (ca_ss - ca_i),
        JRyR ~ R_RYR * (p_co * (ca_sr - ca_ss_co) + p_oo * (ca_sr - ca_ss_oo)),
        ICaL ~ NCaRU / A_CAP * P_LCC * (p_oc * ghk(1, vm, ca_ss_oc, ca_o, 2)+ p_oo * ghk(1, vm, ca_ss_oo, ca_o, 2)),
        JLCC ~ - ICaL * A_CAP / (V_SS_SINGLE * Faraday),
        p_cc ~ y_cc + y_ci + y_ic + y_ii,
        p_co ~ y_co + y_io,
        p_oc ~ y_oc + y_oi,
        p_oo ~ y_oo,
        ca_ss_cc ~ ca_i,
        ca_ss_co ~ (wcai * ca_i + wryr * ca_sr) / (wcai + wryr),
        ca_ss_oc ~ (wcai * ca_i + wlcc * cao)/ (wcai + wlcc),
        ca_ss_oo ~ (wcai * ca_i + wryr * ca_sr + wlcc * cao) / (wcai + wryr + wlcc),
        ca_ss ~ ca_ss_oo * p_oo + ca_ss_co * p_co + ca_ss_oc * p_oc + ca_ss_cc * p_cc,
        α1_L ~ exp(kV_L * (vm - V_L)),
        αp_L ~ it_L * hil(α1_L),
        αm_L ~ it_L * phi_L,
        ϵm_L ~ itau_L * b_L * (α1_L + a_L) / (b_L * α1_L + a_L),
        βm_R ~ it_R * phi_R,
        1 ~ y_oo + y_oc + y_oi + y_co + y_cc + y_ci + y_io + y_ic + y_ii,
        D(y_oo) ~ - v_oo_oc - v_oo_co,
        D(y_oc) ~ v_oo_oc - v_oc_oi - v_oc_cc,
        D(y_oi) ~ v_oc_oi - v_oi_ci,
        D(y_co) ~ v_oo_co - v_co_cc - v_co_io,
        D(y_ci) ~ v_cc_ci + v_oi_ci - v_ci_ii,
        D(y_io) ~ v_co_io - v_io_ic,
        D(y_ic) ~ v_io_ic + v_cc_ic - v_ic_ii,
        D(y_ii) ~ v_ci_ii + v_ic_ii
    ]

    return System(eqs, t; name)
end
