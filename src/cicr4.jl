#===
Coupled LCC-RyR local control model (3 states) by Hinch et al. (2004) https://pmc.ncbi.nlm.nih.gov/articles/PMC1304886/

===#
function get_cicr4_sys(ca_i, ca_sr, ca_o, vm, A_CAP, V_DS=2e-4μm^3; name=:cicrsys)
    @parameters begin
        V_L = -2mV # Potential when half LCC open
        kV_L = inv(7mV) # Steepness of opening potentials
        phi_L = 2.35    # Proportion of time closed in open mode
        t_L = 1ms       # Time switching between C and O states
        it_L = inv(t_L)
        tau_L = 650ms   # LCC inactivation time
        itau_L = inv(tau_L)
        K_L = 0.22μM    # Ca Concentration at LCC inactivation
        iK_L = inv(K_L)
        a_L = 0.0625            # Biasing to make inactivation function of V
        b_L = 14                # Biasing to make inactivation function of_V
        J_L = 9.13E-4μm^3/ms    # Permeability of single LCC
        K_RyR = 41μM    # Ca Concentration at RyR activation
        t_R = 1.17t_L   # Time switching between C and O states in RyR
        it_R = inv(t_R)
        phi_R = 0.05    # Proportion of time closed in open mode in RyR
        tau_R = 2.43ms  # RyR inactivation time
        itauR = inv(tau_R)
        θ_R = 0.012     # reciprocal of proportion of time inactivated in open mode
        c_R = 0.01      # Biasing to make inactivation a function of [Ca2+]ds
        d_R = 100       # Biasing to make inactivation a function of [Ca2+]ds
        J_R = 2E-2μm^3/ms # Permeability of single RyR
        N_CaRU = 50000  # Number of Ca release units
        g_D = 0.065μm^3/m # Ca2+ flux rate from dyadic space to cytosol
    end

    @variables begin
        z1(t)
        z2(t)
        z3(t)
        z4(t) # conserved
        ca_ss_cc(t)
        ca_ss_co(t)
        ca_ss_oc(t)
        ca_ss_oo(t)
        JR_co(t)
        JR_oo(t)
        JL_oc(t)
        JL_oo(t)
        α1_L(t)
        αp_L(t)
        αm_L(t)
        ϵm_L(t)
        βm_R(t)
        yoc_z1(t)
        yco_z1(t)
        yoo_z1(t)
        ycc_z1(t)
        yci_z2(t)
        yoi_z2(t)
        yic_z3(t)
        yio_z3(t)
    end

    _ϵp(ca) = ca * itau_L * iK_L * (α1_L + a_L) / (α1_L + 1)
    _βp(ca) = it_R * hil(ca, K_RyR, 2)
    _μp(ca) = itauR * (ca^2 + c_R * K_RyR^2) / (ca^2 + K_RyR^2)
    _μm(ca) = itauR * θ_R * d_R * (ca^2 + c_R * K_RyR^2) / (d_R * ca^2 + c_R * K_RyR^2)


    wr = J_R / g_D
    wl = J_L / g_D * exprel(-2 * iVT * vm)
    cao = ca_o * exp(-2 * iVT * vm)
    woc = αp_L * βm_R * (αp_L + αm_L + βm_R + _βp(ca_ss_cc))
    wco = αm_L * (_βp(ca_ss_cc) * (αm_L + βm_R + _βp(ca_ss_oc)) + _βp(ca_ss_oc) * αp_L)
    woo = αp_L * (_βp(ca_ss_oc) * (αp_L + βm_R + _βp(ca_ss_cc)) + _βp(ca_ss_cc) * αm_L)
    wcc = αm_L * βm_R * (αp_L + αm_L + βm_R + _βp(ca_ss_oc))
    dem = woc + wco + woo + wcc

    r1 = yoc_z1 * _μp(ca_ss_oc) + ycc_z1 * _μp(ca_ss_cc)
    r2 = (αp_L * _μm(ca_ss_oc) + αm_L * _μm(ca_ss_cc)) / (αp_L + αm_L)
    r3 = _μp(ca_ss_cc) * hil(βm_R, _βp(ca_ss_cc))
    r4 = _μm(ca_ss_cc)
    r5 = yco_z1 * _ϵp(ca_ss_co) + ycc_z1 * _ϵp(ca_ss_cc)
    r6 = ϵm_L
    r7 = _ϵp(ca_ss_cc) * hil(αm_L, αp_L)
    r8 = ϵm_L

    eqs = [
        ca_ss_cc ~ ca_i,
        ca_ss_co ~ (ca_i + wr * ca_sr) / (1 + wr),
        ca_ss_oc ~ (ca_i + wl * cao)/ (1 + wl),
        ca_ss_oo ~ (ca_i + wr * ca_sr + wl * cao) / (1 + wr + wl),
        JR_co ~ J_R * (ca_sr - ca_i) / (1 + wr),
        JL_oc ~ J_L * exprel(-2 * iVT * vm) * (cao - ca_i) / (1 + wl),
        JR_oo ~ J_R * (ca_sr - ca_i + wl * (ca_sr - cao)) / (1 + wr + wl),
        JL_oo ~ J_L * exprel(-2 * iVT * vm) * (cao - ca_i + wr * (cao - ca_sr)) / (1 + wr + wl),
        α1_L ~ exp(kV_L * (vm - V_L)),
        αp_L ~ it_L * hil(α1_L),
        αm_L ~ it_L * phi_L,
        ϵm_L ~ itau_L * b_L * (α1_L + a_L) / (b_L * α1_L + a_L),
        βm_R ~ it_R * phi_R,
        1 ~ z1 + z2 + z3 + z4,
        1 ~ yci_z2 + yoi_z2,
        1 ~ yic_z3 + yio_z3,
        1 ~ yoc_z1 + yco_z1 + yoo_z1 + ycc_z1,
        yci_z2 ~ αm_L / (αp_L + αm_L),
        yic_z3 ~ βm_R / (βm_R + _βp(ca_ss_cc)),
        yoc_z1 ~ woc / dem,
        yco_z1 ~ wco / dem,
        yoo_z1 ~ woo / dem,
        D(z1) ~ -(r1 + r5) * z1 + r2 * z2 + r6 * z3,
        D(z2) ~ r1 * z1 - (r2 + r7) * z2 + r8 * z4,
        D(z3) ~ r5 * z1 - (r6 + r3) * z3 + r4 * z4
    ]

    return ODESystem(eqs, t; name)
end
