using Parameters
# Concentration parameters (in mM)
@with_kw struct Concentrations
    k_o = 5.4  # Extracellular potassium
    na_o = 140.0  # Extracellular sodium
    ca_o = 2.0  # Extracellular calcium
    mg_m = 0.4  # Mitochondrial magnesium (li-2015)
    mg_i = 1.0  # Cytosolic magnesium (Gauthier-2013)
    # mg_i = 3.1  # Cytosolic magnesium (Li-2015)
    pi_i = 3.0  # Cytosolic inorganic phosphate (Gauthier-2013)
    # pi_i = 2.0  # Cytosolic inorganic phosphate (Li-2015)
    PH_I = 7.0  # Cytosolic pH
    h_i = exp10(-PH_I) * 1E3 # Cytosolic proton concentration
    PH_M = 7.6  # Mitocondrial pH (From the Matlab code)
    h_m = exp10(-PH_M) * 1E3
    ΔpH = PH_I - PH_M
    FAC_PH = h_i * 1E4
    AXP_M_T = 1.01  # Mitochondrial ATP + ADP pool (Gauthier-2013)
    AXP_I_T = 8.0  # Cytosolic ATP + ADP pool (Li-2015)
    NAD_T = 1.0  # Mitochondrial NAD + NADH pool (Gauthier-2013)
    # NAD_T = 10.0  # Mitochondrial NAD + NADH pool (Li-2015)
    NADP_M_T = 0.1  # Mitochondrial NADP + NADPH pool (Gauthier-2013)
    nadph_i = 0.075  # Cytosolic NADPH (Gauthier-2013)
    # nadph_i = 1.0  # Cytosolic NADPH (Li-2015)
    O2 = 6E-3 # Oxygen concentration
    pi_m = 8.6512  # Gauthier-2013 and Kembro-2013
end

# Conservation equations
_atp_m(adp_m, c::Concentrations) = c.AXP_M_T - adp_m
_adp_i(atp_i, c::Concentrations) = c.AXP_I_T - atp_i
_nad(nadh, c::Concentrations) = c.NAD_T - nadh
_nadp_m(nadph_m, c::Concentrations) = c.NADP_M_T - nadph_m
