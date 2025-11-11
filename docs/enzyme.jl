#===
# King-Altman method for enzyme kinetics

===#
using Symbolics
using LinearAlgebra

# Make a (n-1)x(n-1) matrix by skipping the i-th row and i-th column
function skip_colrow(mat, i::Int)
    rows = collect(1:size(mat, 1))
    cols = collect(1:size(mat, 2))
    return mat[rows.!=i, cols.!=i]
end

# Accumulate rates into the transition rate matrix
function accumulate_rate!(mat, rate, src::Int, dst::Int)
	mat[dst, src] += rate
	mat[src, src] -= rate
	return mat
end

#===
## SOD model McAdam model

From `McAdam, 1977`
===#

@variables k1 k3 k5 sox

@time sol_sod = let
	k12 = k1 * sox
    k21 = k1 * sox
    k23 = k3 * sox
    k31 = k5
	mat = fill(Num(0), 3, 3)
	accumulate_rate!(mat, k12, 1, 2)
	accumulate_rate!(mat, k21, 2, 1)
	accumulate_rate!(mat, k23, 2, 3)
	accumulate_rate!(mat, k31, 3, 1)
	fA = det(skip_colrow(mat, 1))
    fB = det(skip_colrow(mat, 2))
    fC = det(skip_colrow(mat, 3))

	sol_sod = Dict(
        :A => fA |> expand |> simplify,
        :B => fB |> expand |> simplify,
        :C => fC |> expand |> simplify,
		:DEN => (fA + fB + fC) |> expand |> simplify
    )
end

# Superoxide consumption rate
@time vSOD = sox * (k1 * (sol_sod[:A] + sol_sod[:B]) + k3 * sol_sod[:B]) / (sol_sod[:DEN]) |> simplify

#===
## Complex I simplified Markevich model

Four-state Q-site reaction cycle.
===#

mat_c1q = let
	@variables b12 b21 b23 b32 b34 b43 b41 b14
	mat = fill(Num(0), 4, 4)
	accumulate_rate!(mat, b12, 1, 2)
	accumulate_rate!(mat, b21, 2, 1)
	accumulate_rate!(mat, b23, 2, 3)
	accumulate_rate!(mat, b32, 3, 2)
	accumulate_rate!(mat, b34, 3, 4)
	accumulate_rate!(mat, b43, 4, 3)
	accumulate_rate!(mat, b41, 4, 1)
	accumulate_rate!(mat, b14, 1, 4)
end

#---
@time weights_c1q = [(-1)^(4-1) * det(skip_colrow(mat_c1q, i)) |> expand |> simplify for i in 1:4]
# Total weight
@time total_weight_c1q = sum(weights_c1q) |> expand |> simplify

#===
## Complex I Gauthier model

Based on the `Gauthier, 2013` model
===#

mat_c1g = let
	@variables a12 a21 a65 a56 a61 a16 a23 a32 a34 a43 a47 a74 a57 a75 a42 a24
	mat = fill(Num(0), 7, 7)
	accumulate_rate!(mat, a12, 1, 2)
	accumulate_rate!(mat, a21, 2, 1)
	accumulate_rate!(mat, a65, 6, 5)
	accumulate_rate!(mat, a56, 5, 6)
	accumulate_rate!(mat, a61, 6, 1)
	accumulate_rate!(mat, a16, 1, 6)
	accumulate_rate!(mat, a23, 2, 3)
	accumulate_rate!(mat, a32, 3, 2)
	accumulate_rate!(mat, a34, 3, 4)
	accumulate_rate!(mat, a43, 4, 3)
	accumulate_rate!(mat, a47, 4, 7)
	accumulate_rate!(mat, a74, 7, 4)
	accumulate_rate!(mat, a57, 5, 7)
	accumulate_rate!(mat, a75, 7, 5)
	accumulate_rate!(mat, a42, 4, 2)
	accumulate_rate!(mat, a24, 2, 4)
end

#---
@time weights_c1g = [(-1)^(7-1) * det(skip_colrow(mat_c1g, i)) |> expand |> simplify for i in 1:7]

# Total weight
@time total_weight_c1g = sum(weights_c1g)
