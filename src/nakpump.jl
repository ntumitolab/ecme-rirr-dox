"Na-K pump"
function get_inak_sys(atp_i, adp_i, vm, na_i, na_o, k_o; name=:naksys)
    @parameters KM_NA_NAK = 10.0mM [description = "Na half-saturate constant of Na-K ATPase"]
    @parameters KM_K_NAK = 1.5mM [description = "K half-saturate constant of Na-K ATPase"]
    @parameters IMAX_NAK = 3.147μA / cm²  [description = "Max Na-K ATPase current"]
    @parameters KM_ATP_NAK = 8μM [description = "ATP half-saturate constant of Na-K ATPase"]
    @parameters KI_ADP_NAK = 100μM [description = "ADP half-saturate constant of Na-K ATPase"]
    @variables iNaK(t)
    fko = hil(k_o, KM_K_NAK)
    fnao = expm1(na_o / 67.3mM) / 7
    fnai = hil(na_i, KM_NA_NAK, 1.5)
    f_nak = 1.0 + 0.1245 * exp(-0.1iVT * vm) + 0.0365 * fnao *exp(-iVT * vm)
    f_atp = hil(atp_i * hilr(adp_i, KI_ADP_NAK), KM_ATP_NAK)
    return ODESystem([iNaK ~ IMAX_NAK * fko * fnai * f_atp / f_nak], t; name)
end
