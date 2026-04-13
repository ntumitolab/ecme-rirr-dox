# OXPHOS

function get_etc_eqs(;
    nad_m,                  ## NAD concentration
    nadh_m,                 ## NADH concentration
    dpsi,                   ## Mitochondrial membrane potential
    sox_m,                  ## Superoxide concentration in the matrix
    succinate,              ## Succinate concentration
    fumarate,               ## Fumarate concentration
    oaa,                    ## Oxaloacetate concentration
    MT_PROT=1,                  ## OXPHOS protein content scale factor
    O2=6μM,                     ## Oxygen concentration
    h_i=exp10(-7) * Molar,      ## IMS proton concentration
    h_m=exp10(-7.6) * Molar,    ## Matrix proton concentration
    DOX = 0μM,                  ## Doxorubicin concentration
    ROTENONE_BLOCK=0,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOL_BLOCK=0,
    STIGMATELLIN_BLOCK=0,
    CYANIDE_BLOCK=0,)

    @parameters begin
        Q_T = 4mM        ## CoQ pool
        cytc_T = 325μM   ## Ctyc pool
    end

    @variables begin
        Q_n(t) = 1805μM
        QH2_n(t) = 123μM
        QH2_p(t) = 123μM
        Q_p(t) ## Q_T - Qn - QH2n - QH2p - SQn - SQp
        vSDH(t)
        cytc_rd(t)=0
        cytc_ox(t) ## Conserved
        vHres(t)
        vROS(t)
    end

    @unpack eqs_c1, vQC1, vQH2C1, vNADHC1, vROSC1, vHresC1 = get_c1_eqs(; nad_m, nadh_m, dpsi, sox_m, O2, h_i, h_m, C1_INHIB_DOX=_redoxcycling_c1(DOX), ROTENONE_BLOCK, MT_PROT)

    eqs_c2 = get_c2_eqs(; vSDH, Q_n, QH2_n, succinate, fumarate, oaa, C2_INHIB_DOX=_inhibit_c2(DOX), MT_PROT)

    @unpack eqs_c3, vHresC3, vROSC3 = get_eqs_c3(; vQH2C1, vSDH, cytc_ox, cytc_rd, dpsi, sox_m, C3_INHIB_DOX=_inhibit_c3(DOX), MT_PROT, O2, h_i, h_m, ANTIMYCIN_BLOCK, MYXOTHIAZOL_BLOCK, STIGMATELLIN_BLOCK)

    @unpack eqs_c4, vO2, vHresC4 = get_c4_eqs(; dpsi, cytc_ox, cytc_rd, O2, MT_PROT, C4_INHIB_DOX=_inhibit_c4(DOX), CYANIDE_BLOCK)

    eqs = [
        vROS ~ vROSC1 + vROSC3,
        vHres ~ vHresC1 + vHresC3 + vHresC4,
    ]

    eqs_etc = [eqs_c1; eqs_c2; eqs_c3; eqs_c4; eqs]
    return (; eqs_etc, vHres, vROS, vSDH, vNADHC1, vO2)
end

"Electron transport chain (ETC)"
function get_etc_sys(;
    nad_m,
    nadh_m,
    dpsi,                       # Mitochondrial membrane potential
    sox_m,                      # Superoxide concentration in the matrix
    succinate,                        # Succinate concentration
    fumarate,                        # Fumarate concentration
    oaa,                        # Oxaloacetate concentration
    DOX=0μM,                    # Doxorubicin concentration
    MT_PROT=1,                  # OXPHOS protein content scale factor
    O2=6μM,                     # Oxygen concentration
    h_i=exp10(-7) * Molar,      # IMS proton concentration
    h_m=exp10(-7.6) * Molar,    # Matrix proton concentration
    ROTENONE_BLOCK=0,
    ANTIMYCIN_BLOCK=0,
    MYXOTHIAZOL_BLOCK=0,
    STIGMATELLIN_BLOCK=0,
    CYANIDE_BLOCK=0,
    name=:etcsys)

    @unpack eqs_etc = get_etc_eqs(;
        DOX, MT_PROT, O2, nad_m, nadh_m,
        dpsi, sox_m, succinate, fumarate, oaa, h_i, h_m,
        ROTENONE_BLOCK,
        ANTIMYCIN_BLOCK,
        MYXOTHIAZOL_BLOCK,
        STIGMATELLIN_BLOCK,
        CYANIDE_BLOCK)

    return System(eqs_etc, t; name)
end

function get_c5_eqs(; dpsi, atp_i, adp_i, atp_m, adp_m,
    h_i=exp10(-7) * Molar, h_m=exp10(-7.6) * Molar,
    pi_m=8.6512mM, MT_PROT=1, C5_INHIB=1,
    use_mg=false, mg_i=1mM, mg_m=0.4mM)
    @parameters begin
        ρF1 = 5.0mM                 # Concentration of ATP synthase
        P1_C5 = 1.346E-8
        P2_C5 = 7.739E-7
        P3_C5 = 6.65E-15
        PA_C5 = 1.656E-5Hz
        PB_C5 = 3.373E-7Hz
        PC1_C5 = 9.651E-14Hz
        PC2_C5 = 4.585E-19Hz        ## Magnus model
        ## Equilibrium constant of ATP synthase (ΔG ~ -30kJ/mol)
        ## 1.71E6 * mM in Magnus model was different from Caplan's model (1.71E6 Molar)
        KEQ_C5 = 2E5Molar
        G_H_MITO = 2E-6mM / ms / mV # Proton leak rate constant
        VMAX_ANT = 5E-3mM / ms      # Max rate of ANT, (Wei, 2011)
        H_ANT = 0.5  # Voltage steepness
    end

    @variables begin
        AF1(t)      # Relative electrochemical activity of ATP + H2O <-> ADP + Pi
        vC5(t)      # ATP synthesis rate
        vHu(t)      # Porton flux via ATP synthase
        vANT(t)     # ANT reaction rate
        vHleak(t)   # Proton leak rate
        Δp(t)      # proton motive force
        E_PHOS(t)   # phosphorylation potential
    end

    # F1-Fo ATPase
    if use_mg
        adp3_m, hadp_m, _, poly_adp_m = breakdown_adp(adp_m, h_m, mg_m)
        _, _, mgatp_m, poly_atp_m = breakdown_atp(atp_m, h_m, mg_m)
        poly_pi_m = pipoly(h_m)
        ke_app = KEQ_C5 * poly_atp_m / (poly_adp_m * poly_pi_m)
        v_af1 = ke_app * mgatp_m / (pi_m * (hadp_m + adp3_m))
    else
        v_af1 = KEQ_C5 * atp_m / (pi_m * adp_m)
    end

    vb = exp(iVT * 3 * 50mV) # Boundary potential
    vh = exp(iVT * 3 * dpsi)
    fh = (h_m / h_i)^3
    common = -ρF1 * MT_PROT * C5_INHIB / ((1 + P1_C5 * AF1) * vb + (P2_C5 + P3_C5 * AF1) * vh)
    v_c5 = common * ((PA_C5 * fh + PC1_C5 * vb) * AF1 - (PA_C5 + PC2_C5 * AF1) * vh)
    v_hu = 3 * common * (PA_C5 * fh * AF1 - vh * (PA_C5 + PB_C5))

    # Adenine nucleotide translocator (ANT)
    # Free adenylates
    if use_mg
        atp4_i = atp_i / (1 + h_i / KA_ATP + mg_i / KMG_ATP)
        atp4_m = atp_m / (1 + h_m / KA_ATP + mg_m / KMG_ATP)
        adp3_i = adp_i / (1 + h_i / KA_ADP + mg_i / KMG_ADP)
        adp3_m = adp_m / (1 + h_m / KA_ADP + mg_m / KMG_ADP)
    else
        atp4_i = 0.25 * atp_i
        atp4_m = 0.025 * atp_m
        adp3_i = 0.45 * adp_i
        adp3_m = 0.17 * adp_m
    end

    f_i = atp4_i / adp3_i
    f_m = adp3_m / atp4_m
    v_ant = VMAX_ANT * (1 - f_i * f_m * exp(-iVT * dpsi)) / ((1 + f_i * exp(-iVT * H_ANT * dpsi)) * (1 + f_m))

    eqs_c5 = [
        Δp ~ dpsi + nernst(h_i, h_m, 1),
        E_PHOS ~ VT * NaNMath.log(AF1),
        vANT ~ v_ant,
        AF1 ~ v_af1,
        vHleak ~ G_H_MITO * Δp,
        vC5 ~ v_c5,
        vHu ~ v_hu,
    ]
    return (; eqs_c5, vANT, vC5, vHu, vHleak)
end

function get_c5_sys(; dpsi, h_i, h_m, atp_i, adp_i, atp_m, adp_m, pi_m=8.6512mM, MT_PROT=1, C5_INHIB=1, use_mg=false, mg_i=1mM, mg_m=0.4mM, name=:c5sys)
    @unpack eqs_c5 = get_c5_eqs(;
        dpsi, h_i, h_m, atp_i, adp_i, atp_m, adp_m, pi_m,
        MT_PROT, C5_INHIB, use_mg, mg_i, mg_m)
    return System(eqs_c5, t; name)
end
