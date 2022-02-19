# Sarcolemmal membrane ionic currents 
using Parameters
import .Utils: exprel, expit, p_one, hill, cm², mS, mV, mM, cm, Hz, μM, nernst, nernstNaK

"Fast sodium current (INa)"
@with_kw struct INa{R}
    G::R = 12.8mS / cm²
end

i_na(m, h, j, ΔVNa, G_NA) = G_NA * m^3 * h * j * ΔVNa
i_na(m, h, j, ΔVNa, p::INa) = i_na(m, h, j, ΔVNa, p.G)
i_na(m, h, j, na_i, na_o, vm, p::INa) = i_na(m, h, j, vm - nernst(na_o, na_i), p)
(p::INa)(m, h, j, ΔVNa) = i_na(m, h, j, ΔVNa, p)
(p::INa)(m, h, j, na_i, na_o, vm) = i_na(m, h, j, na_i, na_o, vm, p)

function d_mna(m_na, vm)
    v = vm / mV
    mα = 0.32 / 0.1 * exprel(-0.1 * (v + 47.13))
    mβ = 0.08 * exp(-v / 11)
    return inv(ms) * (mα - m_na * (mα + mβ))
end

function d_hna(h_na, vm)
    v = vm / mV
    ishigh = (v >= -40)
    hαhi = zero(h_na)
    hβhi = 7.6923 * expit((v + 10.66) / 11.1)
    hαlo = 0.135 * exp(-(v + 80) / 6.8)
    hβlo = 3.56 * exp(0.079v) + 3.1E5 * exp(0.35v)
    hα = ifelse(ishigh, hαhi, hαlo)
    hβ = ifelse(ishigh, hβhi, hβlo)
    return inv(ms) * (hα - h_na * (hα + hβ))
end

function d_jna(j_na, vm)
    v = vm / mV
    ishigh = (v >= -40)
    jαhi = zero(j_na)
    jβhi = 0.3 * exp(-2.535e-7v) * expit(0.1 * (v + 32))
    jαlo = (-127140 * exp(0.2444v) - 3.474e-5 * exp(-0.04391v)) * (v + 37.78) * expit(-0.311 * (v + 79.23))
    jβlo = 0.212 * exp(-0.01052v) * expit(0.1378 * (v + 40.14))
    jα = ifelse(ishigh, jαhi, jαlo)
    jβ = ifelse(ishigh, jβhi, jβlo)
    return inv(ms) * (jα - j_na * (jα + jβ))
end

"Potassium currents"
@with_kw struct IK{R}
    k_o::R = 5.4mM
    G_K1::R = 0.748mS / cm² * sqrt(k_o / 5.4mM) # Time-independent
    G_K::R = 0.282mS / cm² * sqrt(k_o / 5.4mM)  # Time-dependent
    P_NA_K::R = 0.01833         # Permeability ration of Na to K in IKs
    G_KP::R = 8.28E-3mS / cm²   # Plateau potassium current
    # ATP-inhibited K channel (Zhou, 2009, from Ferrero et al.)
    FT::R = 1.0266              # Temperature adjustment cons
    KM_NA::R = 25.9mM           # Michaelis constant for sodium
    KM_MG::R = 0.65mM * sqrt((k_o / mM + 5.0))  # Michaelis constant for magnesium
    δMG::R = 0.32               # electrical distance of magnesium
    δNA::R = 0.35               # electrical distance of sodium
    δNAfrt::R = iVT * δNA
    δMGfrt::R = 2iVT * δMG
    σ::R = 0.6 / (μm^2)          # Channel density
    P₀::R = 0.91                # Max channel open probability
    γ₀::R = 35.375E-9mS * (k_o / 5.4mM)^0.24   # Unitary conductance per KATP channel
    G_KATP::R = σ * P₀ * FT * γ₀
end

"Time-independent potassium current"
function i_k1(ΔvK, p::IK)
    v = ΔvK / mV
    α = 1.02 * expit(-0.2385 * (v - 59.215))
    β = (0.4912 * exp(0.28032 * (v + 5.476)) + exp(0.06175 * (v - 594.31))) * expit(0.5143 * (v + 4.753))
    return p.G_K1 * hill(α, β) * ΔvK
end

"Time-dependent delayed rectifier K current system"
function d_xk(x_k, vm)
    v = vm / mV + 30
    α = 7.19e-5 / (0.148) * exprel(-0.148v)
    β = 1.31e-4 / (0.0687) * exprel(0.0687v)
    return inv(ms) * (α - x_k * (α + β))
end

"Time-dependent delayed rectifier K current"
function i_k(x_k, vm, na_i, k_i, na_o, p::IK)
    @unpack G_K, P_NA_K, k_o = p
    ek = nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i)
    x1 = expit(-(vm - 40mV) / 40mV)
    return G_K * x1 * x_k^2 * (vm - ek)
end

"Plateau K curent"
function i_kp(vm, ΔvK, p::IK)
    return p.G_KP * ΔvK * expit((vm - 7.488mV) / 5.98mV)
end

"ATP-dependent K channel (KATP) current"
function i_katp(adp_i, atp_i, vm, na_i, mg_i, ΔvK, p::IK)
    @unpack KM_MG, KM_NA, δMGfrt, δNAfrt, G_KATP = p
    adp = adp_i / mM
    f_mg = hill(KM_MG, mg_i * exp(δMGfrt * vm))       # Inhibition by magnesium
    f_na = hill(KM_NA, na_i * exp(δNAfrt * vm), 2)  # Inhibition by sodium
    km_atp = 35.8mM + 17.9mM * pow_s(adp * inv(μM), 0.256)
    h = 1.3 + 0.74 * exp(-adp)                      # Hill factor
    f_atp = hillr(atp_i, km_atp, h)                 # Inhibition by ATP
    iKatp = G_KATP * f_mg * f_na * f_atp * ΔvK
    return iKatp
end


"Sodium-Calcium exchanger (NCX)"
@with_kw struct NCX{R}
    K_NCX::R = 9000μA / cm²   # Scaling factor of Na/Ca exchange
    KM_NA::R = 87.5mM         # Na half saturating concentration
    KM_CA::R = 1.38mM         # Ca half saturating concentration
    K_SAT::R = 0.1            # Steepness factor
    η::R = 0.35               # Voltage dependence factor
end

"Sodium-Calcium exchanger (NCX) current (iNaCa)"
function i_naca(vm, na_i, ca_i, na_o, ca_o, p::NCX, vfrt = vm * iVT)
    @unpack K_NCX, K_SAT, η, KM_NA, KM_CA = p
    vmax = K_NCX * hill(na_o, KM_NA, 3) * hill(ca_o, KM_CA, 1)
    f_na = (na_i / na_o)^3
    f_ca = ca_i / ca_o
    a_eta = exp(η * vfrt)
    a_etam1 = exp((1 - η) * vfrt)
    return vmax * (a_eta * f_na - a_etam1 * f_ca) * hillr(a_etam1 * K_SAT)
end

(p::NCX)(vm, na_i, ca_i, na_o, ca_o, vfrt = vm * iVT) = i_naca(vm, na_i, ca_i, na_o, ca_o, p, vfrt)

"Non-specific Ca-activated Na current"
@with_kw struct NSNa{R}
    P_NA::R = 1.75E-7cm * Hz
    KM_CA::R = 1.2μM
end

function i_nsna(vm, na_i, ca_i, na_o, p::NSNa)
    @unpack P_NA, KM_CA = p
    return 0.75 * hill(ca_i, KM_CA, 3) * ghkVm(P_NA, vm, na_i, na_o, 1)
end

(p::NSNa)(vm, na_i, ca_i, na_o) = i_nsna(vm, na_i, ca_i, na_o, p)

"Background Na current"
i_nab(ΔvNa, G = 3.22E-3mS / cm²) = G * ΔvNa
i_nab(vm, na_i, na_o, G = 3.22E-3mS / cm²) = i_nab(vm - nernst(na_o, na_i), G)

"Background Ca current"
i_cab(ΔvCa, G = 5.45E-4mS / cm²) = G * ΔvCa
i_cab(vm, ca_i, ca_o, G = 5.45E-4mS / cm²) = i_cab(vm - nernst(ca_o, ca_i, 2), G)

