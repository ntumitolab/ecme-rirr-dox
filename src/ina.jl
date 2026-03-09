function get_ina_eqs(; na_i, na_o, ca_i, ca_o, vm)
    @parameters begin
        G_NA = 12.8mScm⁻²
        P_NSNA = 1.75E-7cm * Hz
        KM_CA_NSNA = 1.2μM
        G_NAB = 3.22E-3mScm⁻²
        G_CAB = 5.45E-4mScm⁻²
    end
    @variables begin
        m_na(t) = 0.0327    # Fast Na gating (activation)
        h_na(t) = 1         # Fast Na gating (inactivation)
        j_na(t) = 1         # Fast Na gating (slow inactivation)
        ENa(t)              # Reversal potential for Na
        ECa(t)              # Reversal potential for Ca
        INa(t)              # Fast Na current
        INsNa(t)            # Non-specific Na current
        INaB(t)             # Background Na current
        ICaB(t)             # Background Ca current
        mα(t)               # Gating variable for INa
        mβ(t)               # Gating variable for INa
        hα(t)               # Gating variable for INa
        hβ(t)               # Gating variable for INa
        jα(t)               # Gating variable for INa
        jβ(t)               # Gating variable for INa
    end

    eqs_ina = [
        INa ~ G_NA * m_na^3 * h_na * j_na * (vm - ENa),
        INsNa ~ 0.75 * hil(ca_i, KM_CA_NSNA, 3) * ghk(P_NSNA, vm, na_i, na_o, 1),
        INaB ~ G_NAB * (vm - ENa),
        ICaB ~ G_CAB * (vm - ECa),
        D(m_na) ~ mα - m_na * (mα + mβ),
        D(h_na) ~ hα - h_na * (hα + hβ),
        D(j_na) ~ jα - j_na * (jα + jβ),
        ENa ~ nernst(na_o, na_i),
        ECa ~ nernst(ca_o, ca_i, 2),
        mα ~ 0.32/ms / 0.1 * exprel(-0.1/mV * (vm + 47.13mV)),
        mβ ~ 0.08/ms * exp(-vm * inv(11mV)),
        hα ~ 0.135 * exp(-(vm + 80mV) * inv(6.8mV)),
        hβ ~ inv(0.13ms) * expit((vm + 10.66mV) * inv(11.1mV)),
        jα ~ max((-127140 * exp(0.2444/mV * vm) - 3.474e-5 * exp(-0.04391/mV * vm)) * (vm + 37.78mV) / (ms * mV) * expit(-0.311/mV * (vm + 79.23mV)), 0),
        jβ ~ 0.3/ms * exp(-2.535e-7/mV * vm) * expit((vm + 32mV) * inv(10mV))
    ]
    return (; eqs_ina, INa, INsNa, INaB, ICaB)
end

"Sodium and background currents"
function get_ina_sys(; na_i, na_o, ca_i, ca_o, vm, name=:inasys)
    @unpack eqs_ina = get_ina_eqs(; na_i, na_o, ca_i, ca_o, vm)
    return System(eqs_ina, t; name)
end
