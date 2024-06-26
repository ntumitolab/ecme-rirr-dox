"Sodium and background currents"
function get_ina_sys(na_i, na_o, ca_i, ca_o, vm; name=:inasys)
    @parameters begin
        G_NA = 12.8mS / cm²
        P_NSNA = 1.75E-7cm * Hz
        KM_CA_NSNA = 1.2μM
        G_NAB = 3.22E-3mS / cm²
    end
    @variables begin
        m_na(t) = 0.0327    # Fast Na gating (activation)
        h_na(t) = 0.9909    # Fast Na gating (inactivation)
        j_na(t) = 0.9941    # Fast Na gating (slow inactivation)
        ΔVNa(t)             # Reversal potential for Na
        ΔVCa(t)             # Reversal potential for Ca
        iNa(t)              # Fast Na current
        iNsNa(t)            # Non-specific Na current
        iNaB(t)             # Background Na current
        iCaB(t)             # Background Ca current
    end
    v = vm / mV
    mα = 0.32 / 0.1 * exprel(-0.1 * (v + 47.13))
    mβ = 0.08 * exp(-v / 11)

    ishigh = (v >= -40)

    hαhi = 0
    hβhi = 7.6923 * expit((v + 10.66) / 11.1)
    hαlo = 0.135 * exp(-(v + 80) / 6.8)
    hβlo = 3.56 * exp(0.079v) + 3.1E5 * exp(0.35v)
    hα = ishigh * hαhi + (1 - ishigh) * hαlo
    hβ = ishigh * hβhi + (1 - ishigh) * hβlo

    jαhi = 0
    jβhi = 0.3 * exp(-2.535e-7v) * expit(0.1 * (v + 32))
    jαlo = (-127140 * exp(0.2444v) - 3.474e-5 * exp(-0.04391v)) * (v + 37.78) * expit(-0.311 * (v + 79.23))
    jβlo = 0.212 * exp(-0.01052v) * expit(0.1378 * (v + 40.14))
    jα = ishigh * jαhi + (1 - ishigh) * jαlo
    jβ = ishigh * jβhi + (1 - ishigh) * jβlo

    # Background Ca current
    @parameters G_CAB = 5.45E-4mS / cm²

    eqs = [
        D(m_na) ~ inv(ms) * (mα - m_na * (mα + mβ)),
        D(h_na) ~ inv(ms) * (hα - h_na * (hα + hβ)),
        D(j_na) ~ inv(ms) * (jα - j_na * (jα + jβ)),
        ΔVNa ~ vm - nernst(na_o, na_i),
        iNa ~ G_NA * m_na^3 * h_na * j_na * ΔVNa,
        iNsNa ~ 0.75 * hil(ca_i^3, KM_CA_NSNA^3) * ghk(P_NSNA, vm, na_i, na_o, 1),
        iNaB ~ G_NAB * ΔVNa,
        ΔVCa ~ vm - nernst(ca_o, ca_i, 2),
        iCaB ~ G_CAB * ΔVCa
    ]
    return ODESystem(eqs, t; name)
end
