# # Enzyme steady-state flux
# Quasi-steady state solutions of enzyme states from transition rates.
using Symbolics
using Groebner

# ## SOD model
# Based on the `McAdam, 1977` model
@variables k1 k3 k5 EA EB EC sox

eqs = let
    vAB = k1 * sox * EA - k1 * sox * EB
    vBC = k3 * sox * EB
    vCA = k5 * EC
    dA = -vAB + vCA
    dB = vAB - vBC
    [dA, dB, EA + EB + EC - 1]
end

@time sol = Symbolics.symbolic_solve(eqs, [EA, EB, EC])[1]

# Superoxide consumption rate
vSOD = sox * (k1 * (sol[EA] + sol[EB]) + k3 * sol[EB])

simplify(vSOD)

# ## Complex I Gauthier model
@variables a12 a21 a65 a56 a61 a16 a23 a32 a34 a43 a47 a74 a57 a75 a42 a24
@variables I1 I2 I3 I4 I5 I6 I7

eqs = let
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
    @assert isequal(d1 + d2 + d3 + d4 + d5 + d6 + d7, 0)
    [d1, d2, d3, d4, d5, d6, I1 + I2 + I3 + I4 + I5 + I6 + I7 - 1]
end

@time sol = Symbolics.symbolic_solve(eqs, [I1, I2, I3, I4, I5, I6, I7])[1]

# Weights of states
for k in (I1, I2, I3, I4, I5, I6, I7)
    println(k, " = ", numerator(sol[k]))
end

# ## Complex I simplified Markevich model

@variables b12 b21 b23 b32 b34 b43 b41 b14
@variables C1 C1_Q C1_SQ C1_QH2

eqs = let
    v12 = C1 * b12 - C1_Q * b21
    v23 = C1_Q * b23 - C1_SQ * b32
    v34 = C1_SQ * b34 - C1_QH2 * b43
    v41 = C1_QH2 * b41 - C1 * b14

    d1 = -v12 + v41
    d2 = v12 - v23
    d3 = v23 - v34
    d4 = v34 - v41

    @assert isequal(sum([d1, d2, d3, d4]), 0)

    [d1, d2, d3, d4, sum([C1, C1_Q, C1_SQ, C1_QH2]) - 1]
end

@time sol = Symbolics.symbolic_solve(eqs, [C1, C1_Q, C1_SQ, C1_QH2])[1]

# Weights of all 7 states
for k in [C1, C1_Q, C1_SQ, C1_QH2]
    println(k, " = ", numerator(sol[k]))
end
