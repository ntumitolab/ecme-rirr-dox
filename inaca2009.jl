using Parameters
include("common.jl")
# Cytosolic Na-Ca exchanger (NCX) parameters
@with_kw struct NCXParams
    na_o = 140.0
    ca_o = 2.0
    # Scaling factor (μA/uF = mV/ms)
    K_NCX = 9000.0
    # Na half saturating concentration (mM)
    KM_NA = 87.5
    # Ca half saturating concentration (mM)
    KM_CA = 1.38
    # Steepness factor
    K_SAT = 0.1
    # Voltage dependence factor
    η = 0.35
    na_o3 = na_o^3
    SCALE = K_NCX / ((KM_NA^3 + na_o3) * (KM_CA + ca_o))
end

# Sodium-Calcium exchanger NCX current (iNaCa) (uA/cm^2)
function inaca(na_i, ca_i, vm, evfrtm1, ca_o, KM_NA, KM_CA, K_NCX, K_SAT, η, SCALE, na_o3)
    eη = _ra(η * vm)
    eηm1 = eη / (evfrtm1 + 1)
    iNaCa = SCALE * (eη * na_i^3 * ca_o - eηm1 * na_o3 * ca_i) / (1 + K_SAT * eηm1)
end

function inaca(na_i, ca_i, vfrt, evfrtm1, ca_o, pNCX::NCXParams)
    @unpack KM_NA, KM_CA, K_NCX, K_SAT, η, SCALE, na_o3 = pNCX
    iNaCa = inaca(na_i, ca_i, vfrt, evfrtm1, ca_o, KM_NA, KM_CA, K_NCX, K_SAT, η, SCALE, na_o3)
end
