using Parameters
# Potassium channels
@with_kw struct IKParams
    k_o = 5.4
    G_KP = 8.28E-3
    P_NA_K = 0.0183
    G_K = 0.282 * sqrt(k_o / 5.4)
    G_K1 = 0.748 * sqrt(k_o / 5.4)
end

# Rate of change of Xk gating variable
function dxk(x_k, vm)
    Δv = vm + 30
    α = -7.9e-5 * Δv / expm1(-0.148 * Δv)
    β = 1.31e-4 * Δv / expm1(0.0687 * Δv)
    d_x_k = α - x_k * (α + β)
end

# Time-dependent K current (iK) (μA/cm²)
function ik(vm, eNaK, x_k, G_K)
    x1_inv = 1 + exp((vm - 40) / 40)
    i_k = G_K * x_k^2 / x1_inv * (vm - eNaK)
end

# Plataeu K current (iKp) (uA/cm^2)
ikp(vm, eK, G_KP) = G_KP * (vm - eK) / (1 + exp((7.488 - vm) / 5.98))

# Time-independent K current (iK1) (uA/cm^2)
function ik1(vm, eK, G_K1)
    Δv = vm - eK
    α = 1.02 / (1 + exp(0.2385 * (Δv - 59.215)))
    β = (0.4912 * exp(0.08032 * (Δv + 5.476)) + exp(0.06175 * (Δv - 594.31))) / (1 + exp(-0.5143 * (Δv + 4.753)))
    ik1 = G_K1 * Δv * α / (α + β)
end
