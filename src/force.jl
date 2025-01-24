"Sarcomere force generation and ATP consumption"
function get_force_sys(atp_i, adp_i, ca_i; name=:forcesys)
    @parameters begin
        ΣLTRPN = 0.07mM             # Total pool of low affinity troponin ca binding sites
        ΣHTRPN = 0.14mM             # Total pool of high affinity troponin ca binding sites
        K_P_LTRPN = 100.0 / (mM * ms)   # Forward constant of low-affinity troponin sites
        K_M_LTRPN = 4e-2 / ms        # Backward constant of low-affinity troponin sites
        K_P_HTRPN = 100.0 / (mM * ms)   # Forward constant of high-affinity troponin sites
        K_M_HTRPN = 3.3e-4 / ms       # Backward constant of high-affinity troponin sites
        K_PN_TROP = 0.04 / ms        # Troponin rate constant (permissive -> nonpermissive)
        SL = 2.15                   # Sarcomere length (μm)
        F_XB = 0.05/ms              # Transition rate from weak to strong cross bridge
        G_MIN = 0.1/ms              # Minimum transition rate from strong to weak cross bridge
        G_OFF = 0.01/ms
        V_MAX_AM = 7.2E-3mM/ms      # Maximal rate of ATP hydrolysis by myofibrils (AM ATPase)
        KM_ATP_AM = 0.03mM          # ATP half saturation constant of AM ATPase
        KI_ADP_AM = 0.26mM          # ADP inhibition constant of AM ATPase
        ζ = 0.1                     # Conversion factor normalizing to  physiological force (N/mm²)
    end

    @variables begin
        p0(t) = 2.601e-5
        p1(t) = 2.248e-5
        p2(t) = 4.199e-5
        p3(t) = 3.657e-5
        n0(t)  # Conserved ~ 1
        n1(t) = 2.243e-5
        force(t)
        force_normal(t)
        Vam(t)
        ltr_ca(t) = 8.949E-3mM
        htr_ca(t) = 0.1321mM
        ltr_free(t) # Conserved
        htr_free(t) # Conserved
        Jtrpn(t)
    end

    G01 = 1 * G_MIN
    G12 = 2 * G_MIN
    G23 = 3 * G_MIN
    F01 = 3 * F_XB
    F12 = 10 * F_XB
    F23 = 7 * F_XB
    ΣPATHS = G01 * G12 * G23 + F01 * G12 * G23 + F01 * F12 * G23 + F01 * F12 * F23
    P1_MAX = F01 * G12 * G23 / ΣPATHS
    P2_MAX = F01 * F12 * G23 / ΣPATHS
    P3_MAX = F01 * F12 * F23 / ΣPATHS
    Φ = 1 + (2.3 - SL) / ((2.3 - 1.7)^1.6)
    G01_SL = Φ * G01
    G12_SL = Φ * G12
    G23_SL = Φ * G23
    G01_OFF = Φ * G_OFF
    N_TROP = 3.5 * SL - 2.0
    K_CA_TRPN = K_M_LTRPN / K_P_LTRPN  # Activation factor for calcium of low affinity troponin sites
    K½_TRPN = hil(1.7μM - 0.8μM * (SL - 1.7) / 0.6, K_CA_TRPN)
    k_np_trop = K_PN_TROP * NaNMath.pow(ltr_ca / (ΣLTRPN * K½_TRPN), N_TROP)

    v01 = F01 * p0 - G01_SL * p1
    v12 = F12 * p1 - G12_SL * p2
    v23 = F23 * p2 - G23_SL * p3
    v04 = K_PN_TROP * p0 - k_np_trop * n0
    v15 = K_PN_TROP * p1 - k_np_trop * n1
    v54 = G01_OFF * n1

    eqs = [
        force ~ ζ * (p1 + n1 + 2 * p2 + 3 * p3) / (P1_MAX + 2P2_MAX + 3P3_MAX),
        force_normal ~ (p1 + p2 + p3 + n1) / (P1_MAX + P2_MAX + P3_MAX),
        Vam ~ V_MAX_AM * hil(atp_i * hil(KI_ADP_AM, adp_i), KM_ATP_AM) / (F01 + F12 + F23) * (F01 * p0 + F12 * p1 + F23 * p2),
        ΣLTRPN ~ ltr_ca + ltr_free,
        ΣHTRPN ~ htr_ca + htr_free,
        D(ltr_ca) ~ K_P_LTRPN * ca_i * ltr_free - K_M_LTRPN * ltr_ca * (1 - 2 / 3 * force_normal),
        D(htr_ca) ~ K_P_HTRPN * ca_i * htr_free - K_M_HTRPN * htr_ca,
        Jtrpn ~ D(ltr_ca) + D(htr_ca),
        1 ~ p0 + p1 + p2 + p3 + n0 + n1,
        D(p0) ~ -v01 - v04,
        D(p1) ~ v01 - v12 - v15,
        D(p2) ~ v12 - v23,
        D(p3) ~ v23,
        D(n1) ~ v15 - v54,
    ]
    return ODESystem(eqs, t; name)
end
