# Calcium buffering dynamics
include("common.jl")

# Calcium buffering parameters
@with_kw struct CABUFParams
    KM_CA
    ET
end

_β_ca(ca, KM, ET) = _mm((ca + KM)^2, KM * ET)
_β_ca(ca, p::CABUFParams) = _β_ca(ca, p.KM_CA, p.ET)
