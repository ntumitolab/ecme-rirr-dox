get_nak_default_params() = ComponentArray(
    KM_NA_NAK = 10mM,  # Na half-saturate concentration of Na-K ATPase
    KM_K_NAK = 1.5mM,    # K half-saturate concentration of Na-K ATPase
    IMAX_NAK = 4.5μAcm⁻², # Max Na-K ATPase current # 3.147μAcm⁻² in (Zhou, 2009)
    KM_ATP_NAK = 8μM,    # ATP half-saturate concentration of Na-K ATPase
    KI_ADP_NAK = 100μM,  # ADP half-saturate concentration of Na-K ATPase
)

"Na-K pump rate"
function vNaK(; na_i, na_o, k_o, atp_i, adp_i, vm, p=get_nak_default_params())
    @unpack KM_NA_NAK, KM_K_NAK, IMAX_NAK, KM_ATP_NAK, KI_ADP_NAK = p
    fnao_NAK = expm1(na_o / 67.3mM) / 7
    fnai = hil(na_i, KM_NA_NAK, 1.5)
    fko_NAK = hil(k_o, KM_K_NAK)
    f_nak = inv(1.0 + 0.1245 * exp(-0.1iVT * vm) + 0.0365 * fnao_NAK * exp(-iVT * vm))
    f_atp = hil(atp_i * hil(KI_ADP_NAK, adp_i), KM_ATP_NAK)
    return IMAX_NAK * fko_NAK * fnai * f_atp * f_nak
end
