function get_ik_eqs(; na_i, na_o, k_i, k_o, vm, atp_i, adp_i, mg_i=1mM)
    fko = sqrt(k_o / 5.4mM)
    @parameters begin
        G_K1 = 0.75mScm⁻² * fko  # Time-independent
        G_K = 0.282mScm⁻² * fko  # Time-dependent
        P_NA_K = 0.01833         # Permeability ratio of Na to K in IKs
        G_KP = 8.28E-3mScm⁻²   # Plateau potassium current
        # ATP-inhibited K channel (Zhou, 2009, adapted from Ferrero et al.)
        G0_KATP = 1.8mScm⁻²   # KATP channel conductance
    end
    @variables begin
        EK(t)               # Reversal potential of K
        IK1(t)              # Time-independent K current
        IK(t)               # Time-dependent delayed rectifier K current
        IKp(t)              # Plateau K current
        IKatp(t)            # ATP-dependent K channel (KATP) current
        x_k(t) = 0          # Time-dependent K channel gating variable
    end
    # Time-independent potassium current
    ΔVK = vm - EK
    vk = ΔVK / mV
    α1 = 1.02 / (1 + exp(0.2385/mV * (vk - 59.215mV)))
    β1 = (0.4912 * exp(0.08032/mV * (vk + 5.476mV)) + exp(0.06175/mV * (vk - 594.31mV))) / (1 + exp(-0.5143/mV * (vk + 4.753mV)))
    # Time-dependent delayed rectifier K current system
    α = 7.19e-5 / ms / (0.148) * exprel(-0.148/mV * (vm + 30mV))
    β = 1.31e-4 / ms / (0.0687) * exprel(0.0687/mV * (vm + 30mV))
    E_KNa = nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i)
    x1 = expit(-inv(40mV) * (vm - 40mV))
    # ATP-dependent K channel (KATP) current
    hatp = 1.3 + 0.74 * exp(-0.09adp_i / μM)   # Hill factor (Ferrero)
    km_atp = 35.8μM + 17.9μM * NaNMath.pow(adp_i / μM, 0.256) # fixed, it's μM rather than mM
    f_atp = hil(km_atp, atp_i, hatp)  # Inhibition by ATP

    eqs_ik = [
        D(x_k) ~ α - x_k * (α + β),
        EK ~ nernst(k_o, k_i),
        IK1 ~ G_K1 * (α1 / (α1 + β1)) * ΔVK,
        IK ~ G_K * x1 * x_k^2 * (vm - E_KNa),
        IKp ~ G_KP * ΔVK / (1 + exp(-(vm - 7.488) / 5.98)),
        IKatp ~ G0_KATP * f_atp * ΔVK
    ]

    return (; eqs_ik, IK1, IK, IKp, IKatp)
end

"Potassium currents, Zhou 2009 version"
function get_ik_sys(; na_i, na_o, k_i, k_o, mg_i, vm, atp_i, adp_i, name=:iksys)
    return create_system(get_ik_eqs; na_i, na_o, k_i, k_o, mg_i, vm, atp_i, adp_i, name)
end
