using Symbolics
using Groebner

#===
## SOD model

Based on the `McAdam, 1977` model
===#
@variables k1 k3 k5 EA EB EC sox

eqs_sod = let
    vAB = k1 * sox * EA - k1 * sox * EB
    vBC = k3 * sox * EB
    vCA = k5 * EC
    dA = -vAB + vCA
    dB = vAB - vBC
    dC = vBC - vCA
    @assert dA + dB + dC == 0
    [dA, dB, EA + EB + EC - 1]
end

#---
@time sol_sod = Symbolics.symbolic_solve(eqs_sod, [EA, EB, EC]) |> only

# Superoxide consumption rate
@time vSOD = sox * (k1 * (sol_sod[EA] + sol_sod[EB]) + k3 * sol_sod[EB]) |> simplify

#===
## Complex I Gauthier model

Based on the `Gauthier, 2013` model
===#
using Symbolics
using Groebner
using Random

@variables a12 a21 a65 a56 a61 a16 a23 a32 a34 a43 a47 a74 a57 a75 a42 a24

@variables I1 I2 I3 I4 I5 I6 I7

vars_c1g = [I1, I2, I3, I4, I5, I6, I7]

eqs_c1g = let I7 = 1
    v12 = I1 * a12 - I2 * a21
    v23 = I2 * a23 - I3 * a32
    v34 = I3 * a34 - I4 * a43
    v47 = I4 * a47 - I7 * a74
    v75 = I7 * a75 - I5 * a57
    v56 = I5 * a56 - I6 * a65
    v61 = I6 * a61 - I1 * a16
    v42 = I4 * a42 - I2 * a24
    d1 = -v12 + v61
    d2 = v12 - v23 + v42
    d3 = v23 - v34
    d4 = v34 - v47 - v42
    d5 = v75 - v56
    d6 = v56 - v61
    d7 = v47 - v75
    @assert d1 + d2 + d3 + d4 + d5 + d6 + d7 == 0
    [d1, d2, d3, d4, d5, d6, sum(vars_c1g)-1]
end

@time sol_c1g = Symbolics.symbolic_solve(eqs_c1g, vars_c1g)[1]
