"Sodium and background currents"
function get_ina_sys(na_i, na_o, ca_i, ca_o, vm; name=:inasys)
    @parameters begin
        G_NA = 12.8mS / cm²
        P_NSNA = 1.75E-7cm * Hz
        KM_CA_NSNA = 1.2μM
        G_NAB = 3.22E-3mS / cm²
        G_CAB = 5.45E-4mS / cm²
    end
    @variables begin
        m_na(t) = 0.0327    # Fast Na gating (activation)
        h_na(t) = 0.9909    # Fast Na gating (inactivation)
        j_na(t) = 0.9941    # Fast Na gating (slow inactivation)
        ΔVNa(t)             # Reversal potential for Na
        ΔVCa(t)             # Reversal potential for Ca
        INa(t)              # Fast Na current
        INsNa(t)            # Non-specific Na current
        INaB(t)             # Background Na current
        ICaB(t)             # Background Ca current
    end
    v = vm / mV
    mα = 0.32 / 0.1 * exprel(-0.1 * (v + 47.13))
    mβ = 0.08 * exp(-v / 11)
    hβ = 7.6923 * expit((v + 10.66) / 11.1)
    hα = 0.135 * exp(-(v + 80) / 6.8)
    hα = ishigh * hαhi + (1 - ishigh) * hαlo
    hβ = ishigh * hβhi + (1 - ishigh) * hβlo
    jβ = 0.3 * exp(-2.535e-7v) * expit(0.1 * (v + 32))
    jα = max((-127140 * exp(0.2444v) - 3.474e-5 * exp(-0.04391v)) * (v + 37.78) * expit(-0.311 * (v + 79.23)), 0)

    eqs = [
        D(m_na) ~ inv(ms) * (mα - m_na * (mα + mβ)),
        D(h_na) ~ inv(ms) * (hα - h_na * (hα + hβ)),
        D(j_na) ~ inv(ms) * (jα - j_na * (jα + jβ)),
        ΔVNa ~ vm - nernst(na_o, na_i),
        INa ~ G_NA * m_na^3 * h_na * j_na * ΔVNa,
        INsNa ~ 0.75 * hil(ca_i^3, KM_CA_NSNA^3) * ghk(P_NSNA, vm, na_i, na_o, 1),
        INaB ~ G_NAB * ΔVNa,
        ΔVCa ~ vm - nernst(ca_o, ca_i, 2),
        ICaB ~ G_CAB * ΔVCa
    ]
    return ODESystem(eqs, t; name)
end
