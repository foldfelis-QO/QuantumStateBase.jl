### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 120da254-1c28-11ec-1938-a55c8cde5727
using Pkg; Pkg.develop(path=".."); Pkg.activate("..")

# ╔═╡ c99f9fc5-0ab1-444a-80f6-b5241edf2f8f
begin
	using QuantumStateBase
	using Plots
	using LinearAlgebra
end

# ╔═╡ db69bedc-2803-4268-9f56-34e4a88a324e
state = SqueezedThermalState(ComplexF, ξ(3, π/4), 2, )

# ╔═╡ cd19d00d-bf2b-479d-8c88-493311d97303
d = gaussian_state_sampler(state, 4096)

# ╔═╡ b9ea40f4-1ed1-4862-8858-be07f2c04633
scatter(d[1, :], d[2, :])

# ╔═╡ ad5fa0c5-0eef-4357-9e14-2970fb046cef
exp(Complex{BigFloat}(5))

# ╔═╡ 7ea5a7c7-cbd9-4ac1-8ea0-f3eb1209e6f9
exp(rand(Complex{BigFloat}, 5, 5))

# ╔═╡ 72fc005e-544c-4dab-9d83-6226248c79e8
?exp

# ╔═╡ 520fef4d-c509-4ff2-9505-a00997a7db2d
LinearAlgebra.BlasFloat

# ╔═╡ 7d5d1fe8-d2a5-4e41-85ad-62b5e918ec27
Complex{BigFloat} isa LinearAlgebra.BlasFloat

# ╔═╡ 572d1bed-bd04-4380-b8df-f9e72c34a206
BigFloat isa LinearAlgebra.BlasFloat

# ╔═╡ Cell order:
# ╟─120da254-1c28-11ec-1938-a55c8cde5727
# ╠═c99f9fc5-0ab1-444a-80f6-b5241edf2f8f
# ╠═db69bedc-2803-4268-9f56-34e4a88a324e
# ╠═cd19d00d-bf2b-479d-8c88-493311d97303
# ╠═b9ea40f4-1ed1-4862-8858-be07f2c04633
# ╠═ad5fa0c5-0eef-4357-9e14-2970fb046cef
# ╠═7ea5a7c7-cbd9-4ac1-8ea0-f3eb1209e6f9
# ╠═72fc005e-544c-4dab-9d83-6226248c79e8
# ╠═520fef4d-c509-4ff2-9505-a00997a7db2d
# ╠═7d5d1fe8-d2a5-4e41-85ad-62b5e918ec27
# ╠═572d1bed-bd04-4380-b8df-f9e72c34a206
