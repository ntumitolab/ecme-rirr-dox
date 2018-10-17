using Parameters
include("common.jl")
# ATP-dependent K channel (KATP) parameters
@with_kw struct KATPParams
    k_o = 5.4
    # Temperature adjustment
    FT = 1.0266
    # Default Michaelis constant for sodium
    KM_NA = 25.9
    # Default Michaelis constant for sodium
    KM_MG = 0.65 * sqrt(k_o + 5.0)
    # electrical distance of magnesium
    δMG = 0.32
    # electrical distance of sodium
    δNA = 0.35
    # Channel density (/cm²)
    σ = 6.0E8
    # Max channel open probability
    P₀ = 0.91
    # Unitary conductance of KATP channel
    γ = 35.375E-9 * (k_o / 5.4)^2.4
    # Mg factor of KATP channel
    IMAX = σ * P₀ * FT * γ
end

#=
ATP-inhibited K current from Gauthier et. al, 2012, based on Ferrero et al. 1996
Scalar version
=#
function ikatp(vm, eK, na_i, atp_i, adp_i, mg_i, KM_MG, KM_NA, δMG, δNA, IMAX)
    f_mg = _mm(KM_MG, mg_i * _ra(2 * δMG * vm))
    f_na = _hills(KM_NA, na_i * _ra(δNA * vm), 2)
    km_atp = 35.8e-3 + 17.9 * 1000^(-0.44) * adp_i^0.56
    h = 1.9 + 0.74 * exp(-adp_i)
    f_atp = _hills(km_atp, atp_i, h)
    iKatp = IMAX * f_mg * f_na * f_atp * (vm - eK)
end

function ikatp(vm, eK, na_i, atp_i, adp_i, mg_i, pKATP::KATPParams)
    @unpack KM_MG, KM_NA, δMG, δNA, IMAX = pKATP
    ikatp(vm, eK, na_i, atp_i, adp_i, mg_i, KM_MG, KM_NA, δMG, δNA, IMAX)
end
