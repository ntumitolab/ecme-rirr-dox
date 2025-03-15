"Na-K pump"
function get_inak_sys(atp_i, adp_i, vm, na_i, na_o, k_o; name=:naksys)
    @parameters begin
        KM_NA_NAK = 10.0mM  # Na half-saturate constant of Na-K ATPase
        KM_K_NAK = 1.5mM    # K half-saturate constant of Na-K ATPase
        IMAX_NAK = 3.147μAcm⁻² # Max Na-K ATPase current
        KM_ATP_NAK = 8μM    # ATP half-saturate constant of Na-K ATPase
        KI_ADP_NAK = 100μM  # ADP half-saturate constant of Na-K ATPase
    end

    @variables INaK(t)
    fko = hil(k_o, KM_K_NAK)
    fnao = expm1(na_o / 67.3mM) / 7
    fnai = hil(na_i, KM_NA_NAK, 1.5)
    f_nak = inv(1.0 + 0.1245 * exp(-0.1iVT * vm) + 0.0365 * fnao * exp(-iVT * vm))
    f_atp = hil(atp_i * hil(KI_ADP_NAK, adp_i), KM_ATP_NAK)
    eqs = [INaK ~ IMAX_NAK * fko * fnai * f_atp * f_nak]
    return ODESystem(eqs, t; name)
end
