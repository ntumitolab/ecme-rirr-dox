# Calcium buffering by CMDN or CSQN
β_ca(ca, KM, ET) = hil((ca + KM)^2, KM * ET)

"Calcium recycle currents"
function get_jca_sys(; atp_i, adp_i, ca_i, ca_nsr, ca_jsr, ca_ss, ca_o, na_i, na_o, vm, name=:caresys)
    @variables begin
        Jup(t)
        IPMCA(t)
        INaCa(t)
        Jtr(t)
        Jxfer(t)
    end

    @parameters begin
        R_TR = inv(9.09ms)          # Diffusion rate between JSR and NSR
        R_XFER= inv(0.5747ms)       # Diffusion rate between subspace and cytosol
        IMAX_PMCA = 0.575μAcm⁻²     # Max PMCA current
        KM_CA_PMCA = 0.5μM          # Ca half-activation constant
        KM1_ATP_PMCA = 12μM         # ATP 1st half-activation constant
        KM2_ATP_PMCA = 230μM        # ATP 2nd half-activation constant
        KI_ADP_PMCA = 1mM           # ADP half-inhibition constant
        VF_SERCA =  0.2989μM/ms     # Max forward SERCA rate
        VR_SERCA =  0.3179μM/ms     # Max reverse SERCA rate
        iKMF_CA_SERCA = inv(0.24μM) # Inverse of Michaelis constant for Ca of forward SERCA reaction
        iKMR_CA_SERCA = inv(1.64269mM) # nverse of Michaelis constant for Ca of reverse SERCA reaction
        NFB_SERCA = 1.4             # Cooperativity of forward SERCA reaction
        NRB_SERCA = 1.0             # Cooperativity of reverse SERCA reaction
        KM_ATP_SERCA = 10μM         # Michaelis constant for ATP in SERCA
        KI1_ADP_SERCA = 140μM       # 1st Michaelis constant for ADP in SERCA
        KI2_ADP_SERCA = 5.1mM       # 2nd Michaelis constant for ADP in SERCA
        K_NCX = 9000μAcm⁻²          # Rate of Na/Ca exchanger (NCX)
        KM_NA_NCX = 87.5mM          # Na half saturating concentration for NCX
        KM_CA_NCX = 1.38mM          # Ca half saturating concentration for NCX
        K_SAT_NCX = 0.1             # Steepness factor for NCX
        η_NCX = 0.35                # Voltage dependence factor for NCX
        IMAX_NCX = K_NCX / (KM_NA_NCX^3 + na_o^3) / (KM_CA_NCX + ca_o) # Precalculated NCX rate
    end

    f_atp_ipca = hil(atp_i * hil(KI_ADP_PMCA, adp_i), KM1_ATP_PMCA) + hil(atp_i, KM2_ATP_PMCA)
    f_ca_ipca = hil(ca_i, KM_CA_PMCA)

    j_up = let
        fb = NaNMath.pow(ca_i * iKMF_CA_SERCA, NFB_SERCA)
        rb = NaNMath.pow(ca_nsr * iKMR_CA_SERCA, NRB_SERCA)
        f_atp = atp_i / (KM_ATP_SERCA * (adp_i / KI1_ADP_SERCA + 1) + atp_i * (adp_i / KI2_ADP_SERCA + 1))
        j_up = f_atp * (VF_SERCA * fb - VR_SERCA * rb) / (fb + rb + 1)
    end

    i_naca = let
        a_eta = exp(η_NCX * iVT * vm)
        a_etam1 = exp((η_NCX - 1) * iVT * vm)
        i_naca = IMAX_NCX * (a_eta * na_i^3 * ca_o - a_etam1 * na_o^3 * ca_i) / (1 + K_SAT_NCX * a_etam1)
    end

    eqs = [
        IPMCA ~ IMAX_PMCA * f_atp_ipca * f_ca_ipca,
        Jup ~ j_up ,
        INaCa ~ i_naca,
        Jtr ~ R_TR * (ca_nsr - ca_jsr),
        Jxfer ~ R_XFER * (ca_ss - ca_i),
    ]
    return System(eqs, t; name)
end
