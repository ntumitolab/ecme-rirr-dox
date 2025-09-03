### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 36058b88-7a73-11f0-393d-bdc74750b193
begin
	import Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
	using Symbolics
    using Groebner
	using PlutoUI
end

# ╔═╡ ba4ad82e-b2a2-4a0d-83f5-b496e2dfa856
PlutoUI.TableOfContents()

# ╔═╡ 85ba8dcd-aba1-4c02-acef-9742a1b6a881
md"""
## SOD model

Based on the `McAdam, 1977` model
"""

# ╔═╡ 61ee1a32-c186-4bf8-99ae-d335621c5d37
@variables k1 k3 k5 EA EB EC sox

# ╔═╡ 717ddb6e-eb89-42ba-bb5a-2344f2b86a9f
eqs_sod = let
    vAB = k1 * sox * EA - k1 * sox * EB
    vBC = k3 * sox * EB
    vCA = k5 * EC
    dA = -vAB + vCA
    dB = vAB - vBC
    [dA, dB, EA + EB + EC - 1]
end

# ╔═╡ fb81a368-ba2d-49ec-b9eb-44e6607cc47f
@time sol_sod = Symbolics.symbolic_solve(eqs_sod, [EA, EB, EC])[1]

# ╔═╡ 8ddb88ca-260b-4b6b-b602-3e60138c808c
# Superoxide consumption rate
vSOD = sox * (k1 * (sol_sod[EA] + sol_sod[EB]) + k3 * sol_sod[EB]) |> simplify

# ╔═╡ 37e15cfc-e18c-4f2d-818d-4625ebb3d1e1
md"""
## Complex I Gauthier model
"""

# ╔═╡ d4e1b99b-b3fc-4487-b04a-2bdea532334c
@variables a12 a21 a65 a56 a61 a16 a23 a32 a34 a43 a47 a74 a57 a75 a42 a24

# ╔═╡ 791ba1c0-b92e-4bad-826f-e68689a1216d
@variables I1 I2 I3 I4 I5 I6 I7

# ╔═╡ 4240a4ba-5585-48c8-a167-1d6759e8168b
eqs_c1g = let
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

# ╔═╡ 6dad921f-5c40-4bd7-b88e-b11b88136e7c
@time sol_c1g = Symbolics.symbolic_solve(eqs_c1g, [I1, I2, I3, I4, I5, I6, I7])[1]

# ╔═╡ cdda990f-b0ca-4bb8-b57c-648b66d62538
# Weights of states
for k in (I1, I2, I3, I4, I5, I6, I7)
    println(k, " = ", numerator(sol_c1g[k]))
end

# ╔═╡ 7553c72f-7667-480d-8709-b188d1a4a251
md"""
## Complex I simplified Markevich model
"""

# ╔═╡ aa1ca2c3-7b60-40cb-a570-c4d4bf88833a
@variables b12 b21 b23 b32 b34 b43 b41 b14

# ╔═╡ 8733db27-b0ca-4d8b-a0ff-5808fbbf22f4
@variables C1 C1_Q C1_SQ C1_QH2

# ╔═╡ ecf33f17-49b8-4b8a-8f3a-0ace9f039dd3
eqs_c1m = let
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

# ╔═╡ 9fcb379c-4b45-4985-8890-3d13628db560
@time sol_c1m = Symbolics.symbolic_solve(eqs_c1m, [C1, C1_Q, C1_SQ, C1_QH2])[1]

# ╔═╡ ddc34d66-4f65-495a-899e-ae23fbc3018c
for k in [C1, C1_Q, C1_SQ, C1_QH2]
    println(k, " = ", numerator(sol_c1m[k]))
end

# ╔═╡ Cell order:
# ╠═36058b88-7a73-11f0-393d-bdc74750b193
# ╠═ba4ad82e-b2a2-4a0d-83f5-b496e2dfa856
# ╠═85ba8dcd-aba1-4c02-acef-9742a1b6a881
# ╠═61ee1a32-c186-4bf8-99ae-d335621c5d37
# ╠═717ddb6e-eb89-42ba-bb5a-2344f2b86a9f
# ╠═fb81a368-ba2d-49ec-b9eb-44e6607cc47f
# ╠═8ddb88ca-260b-4b6b-b602-3e60138c808c
# ╠═37e15cfc-e18c-4f2d-818d-4625ebb3d1e1
# ╠═d4e1b99b-b3fc-4487-b04a-2bdea532334c
# ╠═791ba1c0-b92e-4bad-826f-e68689a1216d
# ╠═4240a4ba-5585-48c8-a167-1d6759e8168b
# ╠═6dad921f-5c40-4bd7-b88e-b11b88136e7c
# ╠═cdda990f-b0ca-4bb8-b57c-648b66d62538
# ╠═7553c72f-7667-480d-8709-b188d1a4a251
# ╠═aa1ca2c3-7b60-40cb-a570-c4d4bf88833a
# ╠═8733db27-b0ca-4d8b-a0ff-5808fbbf22f4
# ╠═ecf33f17-49b8-4b8a-8f3a-0ace9f039dd3
# ╠═9fcb379c-4b45-4985-8890-3d13628db560
# ╠═ddc34d66-4f65-495a-899e-ae23fbc3018c
