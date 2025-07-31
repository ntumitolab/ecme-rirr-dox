function get_force_eqs(; atp_i=7.9mM, adp_i=0.1mM, ca_i=200nM)
    @independent_variables t
    D = Differential(t)
    @parameters begin
        ΣLTRPN = 70μM               # Total pool of low affinity troponin ca binding sites
        ΣHTRPN = 140μM              # Total pool of high affinity troponin ca binding sites
        K_P_LTRPN = 100 / (mM * ms) # Forward constant of low-affinity troponin sites
        K_M_LTRPN = 40Hz            # Backward constant of low-affinity troponin sites
        K_P_HTRPN = 100 / (mM * ms) # Forward constant of high-affinity troponin sites
        K_M_HTRPN = 0.33Hz          # Backward constant of high-affinity troponin sites
        K_PN_TROP = 40Hz            # Troponin rate constant (permissive -> nonpermissive)
        SL = 2.15                   # Sarcomere length (μm)
        F_XB = 50Hz                 # Transition rate from weak to strong cross bridge
        G_MIN = 100Hz               # Minimum transition rate from strong to weak cross bridge
        G_OFF = 10Hz
        V_MAX_AM = 7.2mM * Hz       # Maximal rate of ATP hydrolysis by myofibrils (AM ATPase)
        KM_ATP_AM = 30μM            # ATP half saturation constant of AM ATPase
        KI_ADP_AM = 260μM           # ADP inhibition constant of AM ATPase
        ζ_AM = 0.1                  # Conversion factor normalizing to  physiological force (N/mm²)

        # Dependent parameters
        G01_AM = 1G_MIN
        G12_AM = 2G_MIN
        G23_AM = 3G_MIN
        F01_AM = 3F_XB
        F12_AM = 10F_XB
        F23_AM = 7F_XB
        ΣPATHS_AM = G01_AM * G12_AM * G23_AM + F01_AM * G12_AM * G23_AM + F01_AM * F12_AM * G23_AM + F01_AM * F12_AM * F23_AM
        P1_MAX_AM = F01_AM * G12_AM * G23_AM / ΣPATHS_AM
        P2_MAX_AM = F01_AM * F12_AM * G23_AM / ΣPATHS_AM
        P3_MAX_AM = F01_AM * F12_AM * F23_AM / ΣPATHS_AM
        iDEN_FORCE = inv(P1_MAX_AM + 2P2_MAX_AM + 3P3_MAX_AM)
        iDEN_FORCE_N = inv(P1_MAX_AM + P2_MAX_AM + P3_MAX_AM)
        Φ_AM = 1 + (2.3 - SL) / ((2.3 - 1.7)^1.6) # ~ 1.34
        G01_SL_AM = Φ_AM * G01_AM
        G12_SL_AM = Φ_AM * G12_AM
        G23_SL_AM = Φ_AM * G23_AM
        G01_OFF_AM = Φ_AM * G_OFF
        N_TROP = 3.5 * SL - 2.0   ## Hill coefficient ~ 5.525
        K_CA_TRPN = K_M_LTRPN / K_P_LTRPN  # Activation factor for calcium of low affinity troponin sites
        K½_TRPN = hil((1.7 - 0.8 * (SL - 1.7) / 0.6), K_CA_TRPN)
    end

    @variables begin
        x_p0(t) = 0
        x_p1(t) = 0
        x_p2(t) = 0
        x_p3(t) = 0
        x_n0(t) # Conserved ~ 1
        x_n1(t) = 0
        force(t)
        force_normal(t)
        vAm(t)
        ltr_ca(t) = 9μM
        htr_ca(t) = 132μM
        ltr_free(t) # Conserved = ΣLTRPN - ltr_ca
        htr_free(t) # Conserved = ΣHTRPN - htr_free
        Jtrpn(t)
    end

    k_np_trop = K_PN_TROP * NaNMath.pow(ltr_ca / (ΣLTRPN * K½_TRPN), N_TROP)

    rs = Dict()
    add_rate!(rs, F01_AM, x_p0, G01_SL_AM, x_p1)
    add_rate!(rs, F12_AM, x_p1, G12_SL_AM, x_p2)
    add_rate!(rs, F23_AM, x_p2, G23_SL_AM, x_p3)
    add_rate!(rs, K_PN_TROP, x_p0, k_np_trop, x_n0)
    add_rate!(rs, K_PN_TROP, x_p1, k_np_trop, x_n1)
    add_rate!(rs, G01_OFF_AM, x_n1, 0, x_n0)

    deq = [D(x) ~ rs[x] for x in (x_p0,x_p1,x_p2,x_n1)]
    eqs = [
        force ~ ζ_AM * (x_p1 + x_n1 + 2 * x_p2 + 3 * x_p3) * iDEN_FORCE,
        force_normal ~ (x_p1 + x_p2 + x_p3 + x_n1) * iDEN_FORCE_N,
        vAm ~ V_MAX_AM / (F01_AM + F12_AM + F23_AM) * hil(atp_i * hil(KI_ADP_AM, adp_i), KM_ATP_AM) * (F01_AM * x_p0 + F12_AM * x_p1 + F23_AM * x_p2),
        ΣLTRPN ~ ltr_ca + ltr_free,
        ΣHTRPN ~ htr_ca + htr_free,
        D(ltr_ca) ~ K_P_LTRPN * ca_i * ltr_free - K_M_LTRPN * ltr_ca * (1 - 2 / 3 * force_normal),
        D(htr_ca) ~ K_P_HTRPN * ca_i * htr_free - K_M_HTRPN * htr_ca,
        Jtrpn ~ D(ltr_ca) + D(htr_ca),
        1 ~ x_p0 + x_p1 + x_p2 + x_p3 + x_n0 + x_n1,
    ]
    eqs_force = [deq; eqs]
    return (;eqs_force, vAm, Jtrpn)
end

"Sarcomere force generation and ATP consumption"
function get_force_sys(; atp_i, adp_i, ca_i, name=:forcesys)
    eqs_force, _, _ = get_force_eqs(; atp_i, adp_i, ca_i)
    return System(eqs_force, t; name)
end
