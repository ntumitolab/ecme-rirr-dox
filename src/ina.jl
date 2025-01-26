"Sodium and background currents"
function get_ina_sys(na_i, na_o, ca_i, ca_o, vm; name=:inasys)
    @parameters begin
        G_NA = 12.8milliseimens / cm²
        P_NSNA = 1.75E-7cm * Hz
        KM_CA_NSNA = 1.2μM
        G_NAB = 3.22E-3milliseimens / cm²
        G_CAB = 5.45E-4milliseimens / cm²
    end
    @variables begin
        m_na(t) = 0.0327    # Fast Na gating (activation)
        h_na(t) = 0.9909    # Fast Na gating (inactivation)
        j_na(t) = 0.9941    # Fast Na gating (slow inactivation)
        ENa(t)              # Reversal potential for Na
        ECa(t)              # Reversal potential for Ca
        ΔVNa(t)
        ΔVCa(t)
        INa(t)              # Fast Na current
        INsNa(t)            # Non-specific Na current
        INaB(t)             # Background Na current
        ICaB(t)             # Background Ca current
    end
    v = vm / mV
    mα = 0.32 / 0.1 * exprel(-0.1 * (v + 47.13))
    mβ = 0.08 * exp(-v / 11)
    hα = 0.135 * exp(-(v + 80) / 6.8)
    hβ = 7.6923 * expit((v + 10.66) / 11.1)
    jα = max((-127140 * exp(0.2444v) - 3.474e-5 * exp(-0.04391v)) * (v + 37.78) * expit(-0.311 * (v + 79.23)), 0)
    jβ = 0.3 * exp(-2.535e-7v) * expit(0.1 * (v + 32))

    eqs = [
        INa ~ G_NA * m_na^3 * h_na * j_na * ΔVNa,
        INsNa ~ 0.75 * hil(ca_i^3, KM_CA_NSNA^3) * ghk(P_NSNA, vm, na_i, na_o, 1),
        INaB ~ G_NAB * ΔVNa,
        ICaB ~ G_CAB * ΔVCa,
        D(m_na) ~ mα - m_na * (mα + mβ),
        D(h_na) ~ hα - h_na * (hα + hβ),
        D(j_na) ~ jα - j_na * (jα + jβ),
        ENa ~ nernst(na_o, na_i),
        ΔVNa ~ vm - ENa,
        ECa ~ nernst(ca_o, ca_i, 2),
        ΔVCa ~ vm - ECa,
    ]
    return ODESystem(eqs, t; name)
end
