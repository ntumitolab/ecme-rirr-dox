# ROS transport and scavenging systems

"""
Rate of superoxide dismutase. Based on McAdam, 1977 and Zhou, 2009.
The reaction rate is zero-order due to low apparent Km.
"""
function _vsod(sox, h2o2, K1, K3, K5, KI_H2O2, E0)
    k3′ = K3 * (1 + h2o2 / KI_H2O2)
    denom = K5 * (2K1 + k3′) + k3′ * K1 * sox
    return 2 * E0 * K5 * K1 * sox * (K1 + k3′) / denom
end

get_default_sod_params() = ComponentArray(;
    K1_SOD = 1200 / mM / ms,     # Reaction rate constant of SOD (Zhou, 2009)
    K3_SOD = 24 / mM / ms,       # Inhibition rate constant of SOD (Zhou, 2009)
    K5_SOD = 0.24Hz,             # Recovery rate constant of SOD (Zhou, 2009)
    KI_H2O2_SOD = 500μM,         # H2O2 inhibition constant of SOD (Zhou, 2009)
    ET_SOD_I = 3μM,              # Cytosolic SOD concentration # 1.43μM in (Zhou, 2009)
    ET_SOD_M = 0.3μM             # Mitochondrial SOD concentration # (Kembro, 2013)
)

function vSODi(sox_i, h2o2_i, p=get_default_sod_params())
    @unpack K1_SOD, K3_SOD, K5_SOD, KI_H2O2_SOD, ET_SOD_I = p
    return _vsod(sox_i, h2o2_i, K1_SOD, K3_SOD, K5_SOD, KI_H2O2_SOD, ET_SOD_I)
end

function vSODm(sox_m, h2o2_m, p=get_default_sod_params())
    @unpack K1_SOD, K3_SOD, K5_SOD, KI_H2O2_SOD, ET_SOD_M = p
    return _vsod(sox_m, h2o2_m, K1_SOD, K3_SOD, K5_SOD, KI_H2O2_SOD, ET_SOD_M)
end

"Dalziel type Ping-pong mechanism"
dalziel(A, B, E0, PHI1, PHI2) = E0 * A * B / (PHI1 * A + PHI2 * B)

get_default_gpx_params() = ComponentArray(;
    PHI1_GPX = 5μM * ms,     # Dalziel constant of GPX
    PHI2_GPX = 750μM * ms,   # Dalziel constant of GPX
    ET_GPXI = 10μM,          # GPX concentration (cytosolic)
    ET_GPXM = 10μM,          # GPX concentration (mitochondrial)

)

"Glutathione peroxidase rate (cytosolic)"
vGPXi(h2o2_i, gsh_i, p=get_default_gpx_params()) = dalziel(h2o2_i, gsh_i, p.ET_GPXI, p.PHI2_GPX, p.PHI1_GPX)

"Glutathione peroxidase rate (mitochondrial)"
vGPXm(h2o2_m, gsh_m, p=get_default_gpx_params()) = dalziel(h2o2_m, gsh_m, p.ET_GPXM, p.PHI2_GPX, p.PHI1_GPX)

get_default_tpx_params() = ComponentArray(;
    PHI1_TPX = 3.83mM * ms,  # Dalziel constant of TPX
    PHI2_TPX = 1.85mM * ms,  # Dalziel constant of TPX
    ET_TPXI = 3μM,           # TPX concentration (cytosolic)
    ET_TPXM = 3μM,           # TPX concentration (mitochondrial)
)

"Thioredoxin peroxidase (TPX) rate (cytosolic)"
vTPXi(h2o2_i, trx_i, p=get_default_tpx_params()) = dalziel(h2o2_i, trx_i, p.ET_TPXI, p.PHI2_TPX, p.PHI1_TPX)

"Thioredoxin peroxidase (TPX) rate (mitochondrial)"
vTPXm(h2o2_m, trx_m, p=get_default_tpx_params()) = dalziel(h2o2_m, trx_m, p.ET_TPXM, p.PHI2_TPX, p.PHI1_TPX)

# GR (glutathion reductace) parameters
get_default_gr_params() = ComponentArray(;
    K1_GR = 5Hz,              # Catalytic constant of GR
    ET_GRm = 10μM,            # Mitochondrial concentration of GR
    ET_GRi = 10μM,            # Cytosolic concentration of GR
    KM_NADPH_GR = 15μM,       # Michaelis constant for NADPH of GR
    KM_GSSG_GR = 60μM,        # Michaelis constant for oxidized GSH of GR
    GSHi_Total = 1mM,         # Total glutathione pool (cytosolic)
    GSHm_Total = 1mM,         # Total glutathione pool (mitochondrial)
    nadph_i=75μM              # Cytosolic NADPH concentration
)

"GR (glutathion reductace) rate (cytosolic)"
vGRi(nadph_i, gssg_i, p=get_default_gr_params()) = p.ET_GRi * p.K1_GR * hil(nadph_i, p.KM_NADPH_GR) * hil(gssg_i, p.KM_GSSG_GR)

"GR (glutathion reductace) rate (mitochondrial)"
vGRm(nadph_m, gssg_m, p=get_default_gr_params()) = p.ET_GRm * p.K1_GR * hil(nadph_m, p.KM_NADPH_GR) * hil(gssg_m, p.KM_GSSG_GR)

get_default_tr_params() = ComponentArray(;
    K1_TR = 22.75Hz,          # Catalytic constant of TR
    ET_TRm = 0.35μM,          # Enzyme concentration (mitochondrial)
    ET_TRi = 0.35μM,          # Enzyme concentration (cytosolic)
    KM_NADPH_TR = 65μM,       # Michaelis constant for NADPH
    KM_TRXSS_TR = 35μM,       # Michaelis constant for oxidized thioredoxin of TR
    TRX_T_I = 25μM,           # Thioredoxin pool (cytosolic)
    TRX_T_M = 25μM,           # Thioredoxin pool (mitochondrial)
)

"TR (thioredoxin reductase) rate (cytosolic)"
vTRi(nadph_i, trxss_i, p=get_default_tr_params()) = p.ET_TRi * p.K1_TR * hil(nadph_i, p.KM_NADPH_TR) * hil(trxss_i, p.KM_TRXSS_TR)

"TR (thioredoxin reductase) rate (mitochondrial)"
vTRm(nadph_m, trxss_m, p=get_default_tr_params()) = p.ET_TRm * p.K1_TR * hil(nadph_m, p.KM_NADPH_TR) * hil(trxss_m, p.KM_TRXSS_TR)

get_default_cat_params() = ComponentArray(;
    K1_CAT = 17 / (mM * ms),  # Catalytic constant of Catalase
    ET_CATi = 10μM,           # Cytosolic Catalase concentration
    FR_CAT = 0.05 / mM        # H2O2 inhibition factor of CAT
)

"Catalase rate (cytosolic)"
vCATi(h2o2_i, p=get_default_cat_params()) = 2 * p.K1_CAT * p.ET_CATi * h2o2_i * exp(-p.FR_CAT * h2o2_i)

get_default_imac_params() = ComponentArray(;
    A_IMAC = 0.001,                  # Basal IMAC conductance factor
    B_IMAC = 10000,                  # Activated IMAC conductance factor by cytoplasmic superoxide
    KCC_SOX_IMAC = 10μM,             # Activation constant by cytoplasmic superoxide of IMAC
    GL_IMAC = 3.5E-8mM / ms / mV,    # Leak conductance of IMAC (Zhou, 2009)
    G_MAX_IMAC = GL_IMAC * 100,      # Maximal conductance of IMAC (Zhou, 2009)
    k_IMAC = 0.07 / mV,              # Steepness factor of IMAC
    DPSI_OFFSET_IMAC = 4mV,          # Potential at half saturation
    J_IMAC = 0.5                     # Fraction of superoxide in IMAC conductance
)

"Redox-sensitive IMAC (Inner mitochondrial anion channel) from Cortassa et al. (2004)"
function vIMAC(dpsi, sox_i, sox_m, p=get_default_imac_params())
    @unpack A_IMAC, B_IMAC, KCC_SOX_IMAC, GL_IMAC, G_MAX_IMAC, k_IMAC, DPSI_OFFSET_IMAC, J_IMAC = p
    ΔVROS = nernst(sox_i, sox_m, -1)
    faIMAC = A_IMAC + B_IMAC * hil(sox_i, KCC_SOX_IMAC)
    fvIMAC = GL_IMAC + G_MAX_IMAC / (1 + exp(k_IMAC * (DPSI_OFFSET_IMAC + dpsi)))
    gIMAC = fvIMAC * faIMAC
    return (; ΔVROS=ΔVROS, vTrROS=J_IMAC * gIMAC * (dpsi + ΔVROS), vIMAC=gIMAC * dpsi)
end

#===
# ODEs

V_MITO_V_MYO=0.615
D(sox_i) ~ V_MITO_V_MYO * vTrROS - vSOD_i,
D(h2o2_i) ~ 0.5vSOD_i - vGPX_i - vCAT,
D(gssg_i) ~ -0.5 * (vGR_i - vGPX_i)
ΣGSH_i ~ gsh_i + 2gssg_i,
===#
