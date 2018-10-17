# Fast sodium current

# Rate of change of gating variables
function d_mna(m_na, vm)
    mα = 0.32 * (vm + 47.13) / (-expm1(-0.1 * (vm + 47.13)))
    mβ = 0.08 * exp(-vm / 11.0)
    d_m_na = mα - m_na * (mα + mβ)
end

function d_hna(h_na, vm)
    if vm >= -40
        hα = 0.0
        hβ = inv(0.13 * (1 + exp((vm + 10.66) / -11.1)))
    else
        hα = 0.135 * exp(-(vm + 80.0) / 6.8)
        hβ = 3.56 * exp(0.079 * vm) + 3.1e5 * exp(0.35 * vm)
    end
    d_h_na = hα - h_na * (hα + hβ)
end

function d_jna(j_na, vm)
    if vm >= -40
        jα = 0.0
        jβ = 0.3 * exp(-2.535e-7 * vm) / (1 + exp(-0.1 * (vm + 32)))
    else
        jα = ((-127140 * exp(0.2444 * vm) - 3.471e-5 * exp(-0.04391 * vm))
              * (vm + 37.78) / (1 + exp(0.311 * (vm + 79.23))))
        jβ = 0.1212 * exp(-0.01052 * vm) / (1 + exp(-0.1378 * (vm + 40.14)))
    end
    d_j_na = jα - j_na * (jα + jβ)
end

# Fast sodium current (iNa)
ina(vm, eNa, m_na, h_na, j_na, G_NA) = G_NA * (vm - eNa) * m_na^3 * h_na * j_na
