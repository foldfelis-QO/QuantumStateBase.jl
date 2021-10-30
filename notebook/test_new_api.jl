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
state = SqueezedThermalState(Complex{Float64}, ξ(0.5, π/4), 0.5, )

# ╔═╡ cd19d00d-bf2b-479d-8c88-493311d97303
d = gaussian_state_sampler(Float64, state, 4096)

# ╔═╡ b9ea40f4-1ed1-4862-8858-be07f2c04633
scatter(d[1, :], d[2, :])

# ╔═╡ Cell order:
# ╟─120da254-1c28-11ec-1938-a55c8cde5727
# ╠═c99f9fc5-0ab1-444a-80f6-b5241edf2f8f
# ╠═db69bedc-2803-4268-9f56-34e4a88a324e
# ╠═cd19d00d-bf2b-479d-8c88-493311d97303
# ╠═b9ea40f4-1ed1-4862-8858-be07f2c04633
