# Calcium buffer by CMDN and CSQN
β_ca(ca, KM, ET) = hil((ca + KM)^2, KM * ET)
β_ca_cmdn(ca, KM_CA_CMDN=2.38μM, ET_CMDN=50μM) = β_ca(ca, KM_CA_CMDN, ET_CMDN)
β_ca_csqn(ca, KM_CA_CSQN=0.8mM, ET_CSQN=5.0mM) = β_ca(ca, KM_CA_CSQN, ET_CSQN)

"Calcium recycle currents"
function get_jca_sys(atp_i, adp_i, ca_i, ca_nsr, ca_jsr, ca_ss, ca_o, na_i, na_o, vm; name=:caresys)
    @variables begin
        jUp(t)
        iPMCA(t)
        iNaCa(t)
        jTr(t)
        jXfer(t)
    end

    # Diffusion rate between JSR and NSR
    @parameters R_TR = inv(9.09ms)
    # Diffusion rate between subspace and cytosol
    @parameters R_XFER= inv(0.5747ms)
    # Plasma membrane calcium pump (PMCA)
    @parameters IMAX_PMCA = 0.575μA / cm² [description = "Max PMCA current"]
    @parameters KM_CA_PMCA = 0.5μM    [description = "Ca half-saturation constant"]
    @parameters KM1_ATP_PMCA = 12μM [description = "ATP 1st half-saturation constant"]
    @parameters KM2_ATP_PMCA = 230μM  [description = "ATP 2nd half-saturation constant"]
    @parameters KI_ADP_PMCA = 1.0mM [description = "ADP half-inhibition constant"]
    ipca = let
        f_atp = hil(atp_i * hilr(adp_i, KI_ADP_PMCA), KM1_ATP_PMCA) + hil(atp_i, KM2_ATP_PMCA)
        f_ca = hil(ca_i, KM_CA_PMCA)
        IMAX_PMCA * f_atp * f_ca
    end

    # SER Calcium pump (SERCA)
    @parameters VF_SERCA =  2.989E-4mM * kHz [description = "Max forward SERCA rate"]
    @parameters VR_SERCA =  3.179E-4mM * kHz [description = "Max reverse SERCA rate"]
    @parameters KMF_CA_SERCA = 0.24μM [description = "Michaelis constant for Ca of forward SERCA reaction"]
    @parameters KMR_CA_SERCA = 1.64269mM [description =  "Michaelis constant for Ca of reverse SERCA reaction"]
    @parameters NFB_SERCA = 1.4 [description = "Cooperativity of forward SERCA reaction"]
    @parameters NRB_SERCA = 1.0 [description = "Cooperativity of reverse SERCA reaction"]
    @parameters KM_ATP_SERCA = 10μM [description = "Michaelis constant for ATP in SERCA"]
    @parameters KI1_ADP_SERCA = 140μM [description = "1st Michaelis constant for ADP in SERCA"]
    @parameters KI2_ADP_SERCA = 5.1mM [description = "2nd Michaelis constant for ADP in SERCA"]
    jup = let
        fb = NaNMath.pow(ca_i / KMF_CA_SERCA, NFB_SERCA)
        rb = NaNMath.pow(ca_nsr / KMR_CA_SERCA, NRB_SERCA)
        f_atp_inv = KM_ATP_SERCA / atp_i * (adp_i / KI1_ADP_SERCA + 1) + (adp_i / KI2_ADP_SERCA + 1)
        (VF_SERCA * fb - VR_SERCA * rb) / ((fb + rb + 1) * f_atp_inv)
    end

    # Sodium-Calcium exchanger (NCX)
    @parameters begin
        K_NCX = 9000μA / cm²   # Scaling factor of Na/Ca exchange
        KM_NA_NCX = 87.5mM         # Na half saturating concentration
        KM_CA_NCX = 1.38mM         # Ca half saturating concentration
        K_SAT_NCX = 0.1            # Steepness factor
        η_NCX = 0.35               # Voltage dependence factor
    end

    inaca = let
        vmax_ncx = K_NCX / (KM_NA_NCX^3 + na_o^3) / (KM_CA_NCX + ca_o)
        a_eta = exp(η_NCX * vm * iVT)
        a_etam1 = exp((η_NCX - 1) * vm * iVT)
        num = a_eta * na_i^3 * ca_o - a_etam1 * na_o^3 * ca_i
        vmax_ncx * num / (1 + K_SAT_NCX * a_etam1)
    end

    eqs = [
        iPMCA ~ ipca,
        jUp ~ jup,
        iNaCa ~ inaca,
        jTr ~ R_TR * (ca_nsr - ca_jsr),
        jXfer ~ R_XFER * (ca_ss - ca_i),
    ]
    return ODESystem(eqs, t; name)
end
