using Parameters
include("common.jl")

# Background current parameters
@with_kw struct BGParams
    G_CAB = 3.22E-3  # Conductance of Background Calcium current
    G_NAB = 5.45E-4  # Conductance of Background Sodium current
end

# Background Ca current (iCaB) (μA/cm²)
icab(vm, eCa, G_CAB) = G_CAB * (vm - eCa)

# Background Na current (iNaB) (μA/cm²)
inab(vm, eNa, G_NAB) = G_NAB * (vm - eNa)

# Non-specific Ca-activated current parameters
@with_kw struct NSParams
    # Na permebility of non-specific current (cm/s)
    P_NA = 1.75E-7
    # K permebility of non-specific current (cm/s)
    P_K = 0.
    # Ca Michaelis constant for non-specific current (mM)
    KM_CA = 1.2E-3
end

# Non-specific Ca-activated Na current (iNsNa) (uA/cm^2)"
function insna(na_i, ca_i, vfrt, evfrtm1, na_o, KM_CA, P_NA)
    imax = 0.75 * P_NA * F * vfrt * (na_i * (evfrtm1 + 1) - na_o) / evfrtm1
    insna = imax * _hills(ca_i, KM_CA, 3)
end

function insna(na_i, ca_i, vfrt, evfrtm1, na_o, pNS::NSParams)
    @unpack KM_CA, P_NA = pNS
    insna(na_i, ca_i, vfrt, evfrtm1, na_o, KM_CA, P_NA)
end
