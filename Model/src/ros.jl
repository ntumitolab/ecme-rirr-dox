"""
Rate of superoxide dismutase. Based on (McAdam, 1977)
The activity will be capped at high concentrations of SOX.
"""
function _vsod(sox, h2o2, K1, K3, K5, KI_H2O2, E0)
    k3′ = K3 * (1 + h2o2 / KI_H2O2)
    denom = K5 * (2K1 + k3′) + k3′ * K1 * sox
    return 2 * E0 * K5 * K1 * sox * (K1 + k3′) / denom
end

function get_ros_eqs(; dpsi, sox_m, nadph_i=75μM, V_MITO_V_MYO=0.615)
    @parameters begin
        # superoxide dismutase (SOD)
        K1_SOD = 1200 / mM / ms     # Reaction rate constant of SOD (Zhou, 2009)
        K3_SOD = 24 / mM / ms       # Inhibition rate constant of SOD (Zhou, 2009)
        K5_SOD = 0.24Hz             # Recovery rate constant of SOD (Zhou, 2009)
        KI_H2O2_SOD = 500μM         # H2O2 inhibition constant of SOD (Zhou, 2009)
        ET_SOD_I = 3μM              # Cytosolic SOD concentration # 1.43μM in (Zhou, 2009)
        ET_SOD_M = 0.3μM            # Mitochondrial SOD concentration # (Kembro, 2013)
        # glutathione peroxidase (GPX)
        𝚽1_GPX = 5μM * ms  # Rate constant of GPX
        𝚽2_GPX = 750μM * ms  # Rate constant of GPX
        ET_GPX = 10μM  # GPX concentration ()
        # thioredoxin peroxidase (TPX)
        𝚽1_TPX = 3.83mM * ms  # Rate constant of TPX
        𝚽2_TPX = 1.85mM * ms  # Rate constant of TPX
        ET_TPX = 3μM  # TPX concentration
        # GR (glutathion reductace) parameters
        K1_GR = 5Hz  # Catalytic constant of GR
        ET_GR = 10μM  # Mitochondrial concentration of GR
        KM_NADPH_GR = 15μM  # Michaelis constant for NADPH of GR
        KM_GSSG_GR = 60μM  # Michaelis constant for oxidized GSH of GR
        ΣGSH_i = 1mM  # Cytosolic glutathione pool (mM)
        ΣGSH_m = 1mM  # Mitochondrial glutathione pool (mM)
        # Thioredoxin reductase (TR) parameters
        K1_TR = 22.75Hz         # Catalytic constant of TR
        ET_TR = 0.35μM          # Enzyme concentration
        KM_NADPH_TR = 65μM      # Michaelis constant for NADPH
        KM_TRXSS_TR = 35μM      # Michaelis constant for oxidized thioredoxin of TR
        TRX_T_TR = 25μM         # Thioredoxin pool(mM)
        # Catalase
        K1_CAT = 17 / (mM * ms)  # Catalytic constant of Catalase
        ET_CAT = 10μM            # Catalase concentration
        FR_CAT = 0.05 / mM       # H2O2 inhibition factor of CAT
        # IMAC (Inner mitochondrial anion channel) from Cortassa et al. (2004)
        A_IMAC = 0.001                  # Basal IMAC conductance factor
        B_IMAC = 10000                  # Activated IMAC conductance factor by cytoplasmic superoxide
        KCC_SOX_IMAC = 10μM             # Activation constant by cytoplasmic superoxide of IMAC
        GL_IMAC = 3.5E-8mM / ms / mV    # Leak conductance of IMAC (Zhou, 2009)
        G_MAX_IMAC = GL_IMAC * 100      # Maximal conductance of IMAC (Zhou, 2009)
        k_IMAC = 0.07 / mV             # Steepness factor of IMAC (some papers say it's +0.07/mV)
        DPSI_OFFSET_IMAC = 4mV          # Potential at half saturation
        J_IMAC = 0.5                    # Fraction of superoxide in IMAC conductance
    end

    @variables begin
        sox_i(t) = 1nM
        h2o2_i(t) = 0.76nM
        # h2o2_m(t)
        gsh_i(t) # Conserved
        # gsh_m(t)
        gssg_i(t) = 2.22μM
        # gssg_m(t)
        # nadph_m(t)
        vSOD_i(t)
        vGPX_i(t)
        vGR_i(t)
        vCAT(t)     ## Catalase flux
        vTrROS(t)   ## SOX flux via IMAC
        vIMAC(t)    ## IMAC ion flux
        gIMAC(t)    ## IMAC conductance
        fvIMAC(t)   ## IMAC activated by voltage
        faIMAC(t)   ## IMAC activated by SOX
        ΔVROS(t)    ## Nernst potential of SOX
    end
    eqs_ros = [
        ΔVROS ~ nernst(sox_i, sox_m, -1),
        vTrROS ~ J_IMAC * gIMAC * (dpsi + ΔVROS),
        vIMAC ~ gIMAC * dpsi,
        gIMAC ~ fvIMAC * faIMAC,
        fvIMAC ~ GL_IMAC + G_MAX_IMAC / (1 + exp(k_IMAC * (DPSI_OFFSET_IMAC + dpsi))),
        faIMAC ~ A_IMAC + B_IMAC * hil(sox_i, KCC_SOX_IMAC),
        vGR_i ~ ET_GR * K1_GR * hil(nadph_i, KM_NADPH_GR) * hil(gssg_i, KM_GSSG_GR),
        vGPX_i ~ ET_GPX * h2o2_i * gsh_i / (𝚽1_GPX * gsh_i + 𝚽2_GPX * h2o2_i),
        vCAT ~ 2K1_CAT * ET_CAT * h2o2_i * exp(-FR_CAT * h2o2_i),
        vSOD_i ~ _vsod(sox_i, h2o2_i, K1_SOD, K3_SOD, K5_SOD, KI_H2O2_SOD, ET_SOD_I),
        ΣGSH_i ~ gsh_i + 2gssg_i,
        D(sox_i) ~ V_MITO_V_MYO * vTrROS - vSOD_i,
        D(h2o2_i) ~ 0.5vSOD_i - vGPX_i - vCAT,
        D(gssg_i) ~ -0.5 * (vGR_i - vGPX_i)
    ]
    return (; eqs_ros, vTrROS, vIMAC, vCAT, vGPX_i, vGR_i, vSOD_i)
end

"ROS diffusion and detox system"
function get_ros_sys(; dpsi, sox_m, nadph_i=75μM, V_MITO_V_MYO=0.615, name=:rossys)
    @unpack eqs_ros = get_ros_eqs(; dpsi, sox_m, nadph_i, V_MITO_V_MYO)
    return System(eqs_ros, t; name)
end
