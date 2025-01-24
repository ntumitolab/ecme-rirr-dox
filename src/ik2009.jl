"Potassium currents"
function get_ik_sys(na_i, na_o, k_i, k_o, mg_i, vm, atp_i, adp_i; name=:iksys)
    @parameters begin
        G_K1 = 0.75mS / cm² * sqrt(k_o / 5.4mM)  # Time-independent
        G_K = 0.282mS / cm² * sqrt(k_o / 5.4mM)  # Time-dependent
        P_NA_K = 0.01833         # Permeability ratio of Na to K in IKs
        G_KP = 8.28E-3mS / cm²   # Plateau potassium current
        # ATP-inhibited K channel (Zhou, 2009, adapted from Ferrero et al.)
        G0_KATP = 1.8mS / cm²   # KATP channel conductance
    end
    @variables begin
        ΔVK(t)              # Reversal potential of K
        IK1(t)              # Time-independent K current
        IK(t)               # Time-dependent delayed rectifier K current
        IKp(t)              # Plateau K current
        IKatp(t)            # ATP-dependent K channel (KATP) current
        x_k(t) = 1.1212e-4  # Time-dependent K channel gating variable
    end

    # Time-independent potassium current
    vk = ΔVK / mV
    α1 = expit(-0.2385 * (vk - 59.215), 1.02)
    β1 = (0.4912 * exp(0.08032 * (vk + 5.476)) + exp(0.06175 * (vk - 594.31))) * expit(0.5143 * (vk + 4.753))
    # Time-dependent delayed rectifier K current system
    α = 7.19e-5 / ms / (0.148) * exprel(-0.148 * (v + 30))
    β = 1.31e-4 / ms / (0.0687) * exprel(0.0687* (v + 30))
    EKNa = nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i)
    x1 = expit(-(vm - 40) / 40)
    # ATP-dependent K channel (KATP) current
    hatp = 1.3 + 0.74 * exp(-0.09 * adp_i / μM)   # Hill factor (Ferrero)
    km_atp = 35.8μM + 17.9μM * NaNMath.pow(adp_i / μM, 0.56) # fixed, it's μM rather than mM
    f_atp = hil(km_atp, atp_i, hatp)  # Inhibition by ATP

    eqs = [
        D(x_k) ~ α - x_k * (α + β),
        ΔVK ~ vm - nernst(k_o, k_i),
        IK1 ~ G_K1 * (α1 / (α1 + β1)) * ΔVK,
        IK ~ G_K * x1 * x_k^2 * (vm - EKNa),
        IKp ~ G_KP * ΔVK * expit((vm - 7.488mV) / 5.98mV),
        IKatp ~ G0_KATP * f_atp * ΔVK
    ]
    return ODESystem(eqs, t; name)
end
