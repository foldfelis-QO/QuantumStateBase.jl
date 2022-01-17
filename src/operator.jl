export
    create,
    annihilate,
    displace,
    squeeze

function check_finite_size(s::AbstractArray)
    ℵ₀ in size(s) && throw(BoundsError("infinite size"))

    return true
end

############
# a† and a #
############

Creation(T::Type{<:Number}) = .√(Diagonal(T.(0:∞))[:, 2:end])
Creation(T, dim::Integer) = view(Creation(T), 1:dim, 1:dim)

function create(v::AbstractVector{T}) where {T}
    dim = length(v)

    check_finite_size(v) && (return Creation(T, dim) * v)
end

function create(ρ::AbstractMatrix{T}) where {T}
    dim = size(ρ, 1)
    c = Creation(T, dim)

    check_finite_size(ρ) && (return c * ρ * c')
end

Annihilation(T::Type{<:Number}) = .√(Diagonal(T.(0:∞))[2:end, :])
Annihilation(T, dim::Integer) = view(Annihilation(T), 1:dim, 1:dim)

function annihilate(v::AbstractVector{T}) where {T}
    dim = length(v)

    check_finite_size(v) && (return Annihilation(T, dim) * v)
end

function annihilate(ρ::AbstractMatrix{T}) where {T}
    dim = size(ρ, 1)
    a = Annihilation(T, dim)

    check_finite_size(ρ) && (return a * ρ * a')
end

###########
# α and ξ #
###########

struct ComplexVec{T<:Real}
    r::T
    θ::T
end

ComplexVec(r::T, θ::T) where {T} = ComplexVec{T}(r, θ)

function Base.show(io::IO, complexvec::ComplexVec{T}) where {T}
    print(io, "ComplexVec{$T}($(complexvec.r)exp(-$(complexvec.θ)im))")
end

z(complexvec::ComplexVec) = complexvec.r * exp(-im * complexvec.θ)

################
# displacement #
################

function Displacement(T::Type{<:Complex}, r::Real, θ::Real, dim)
    U = real(T)
    α = ComplexVec(U(r), U(θ))

    return exp(collect( # `collect` applied due to limitation of `exp(A::AbstractMatrix)`
        z(α) * Creation(U, dim) - z(α)' * Annihilation(U, dim)
    ))
end

function displace(v::AbstractVector{T}, r, θ) where {T}
    dim = length(v)

    check_finite_size(v) && (return Displacement(complex(T), r, θ, dim) * v)
end

function displace(ρ::AbstractMatrix{T}, r, θ) where {T}
    dim = size(ρ, 1)
    d = Displacement(complex(T), r, θ, dim)

    check_finite_size(ρ) && (return d * ρ * d')
end

#############
# squeezing #
#############

function Squeezing(T::Type{<:Complex}, r::Real, θ::Real, dim)
    U = real(T)
    ξ = ComplexVec(U(r), U(θ))

    return exp(collect( # `collect` applied due to limitation of `exp(A::AbstractMatrix)`
        (z(ξ)' * Annihilation(U, dim)^2)/2 - (z(ξ) * Creation(U, dim)^2)/2
    ))
end

function squeeze(v::AbstractVector{T}, r, θ) where {T}
    dim = length(v)

    check_finite_size(v) && (return Squeezing(complex(T), r, θ, dim) * v)
end

function squeeze(ρ::AbstractMatrix{T}, r, θ) where {T}
    dim = size(ρ, 1)
    s = Squeezing(complex(T), r, θ, dim)

    check_finite_size(s) && (return s * ρ * s')
end

###################
# BHD measurement #
###################

# ##### for arb. state in intensity-to-measurement-phase quadrature coordinate #####

# # |θ, x⟩ = ∑ₙ |n⟩ ⟨n|θ, x⟩ = ∑ₙ ψₙ(θ, x) |n⟩
# # ⟨n|θ, x⟩ = ψₙ(θ, x) = exp(im n θ) (2/π)^(1/4) exp(-x^2) Hₙ(√2 x)/√(2^n n!)
# function ψₙ(n::Integer, θ::Real, x::Real)
#     return (2/π)^(1/4) * exp(im*n*θ - x^2) * hermiteh(n, sqrt(2)x) / sqrt(2^n * factorial(n))
# end

# function 𝛑̂!(result::Matrix{<:Complex}, θ::Real, x::Real; dim=DIM)
#     view(result, :, 1) .= ψₙ.(big(0):big(dim-1), θ, x)
#     result .= view(result, :, 1) * view(result, :, 1)'

#     return result
# end

# function 𝛑̂(T::Type{<:Complex}, θ::Real, x::Real; dim=DIM)
#     result = Matrix{T}(undef, dim, dim)

#     return 𝛑̂!(result, θ, x, dim=dim)
# end

# 𝛑̂(θ::Real, x::Real; dim=DIM) = 𝛑̂(ComplexF64, θ, x, dim=dim)

# ##### for Gaussian state in intensity-to-measurement-phase quadrature coordinate #####

# # π̂ₓ = (â exp(-im θ) + â† exp(im θ)) / 2

# tr_mul(𝐚, 𝐛) = sum(𝐚[i, :]' * 𝐛[:, i] for i in 1:size(𝐚, 1))
# create_μ(state::StateMatrix{T}) where {T} = tr_mul(Creation(T, dim=state.dim), state.𝛒)
# create²_μ(state::StateMatrix{T}) where {T} = tr_mul(Creation(T, dim=state.dim)^2, state.𝛒)
# annihilate_μ(state::StateMatrix{T}) where {T} = tr_mul(Annihilation(T, dim=state.dim), state.𝛒)
# annihilate²_μ(state::StateMatrix{T}) where {T} = tr_mul(Annihilation(T, dim=state.dim)^2, state.𝛒)
# create_annihilate_μ(state::StateMatrix{T}) where {T} = tr_mul(
#     Creation(T, dim=state.dim) * Annihilation(T, dim=state.dim),
#     state.𝛒
# )

# # ⟨π̂ₓ²⟩ = ⟨ââ exp(-2im θ) + â†â† exp(2im θ) + ââ† + â†â⟩ / 4
# # ⟨π̂ₓ²⟩ = (exp(-2im θ)⟨â²⟩ + exp(2im θ)⟨â†²⟩ + 1 + 2⟨ââ†⟩) / 4
# # here, ⟨ââ† + â†â⟩ = 1 + 2⟨ââ†⟩ due to the commutation relation
# function π̂ₓ²_μ(θs::AbstractVector{<:Number}, state::StateMatrix)
#     return (
#         exp.(-2im*θs) .* annihilate²_μ(state) .+
#         exp.(2im*θs) .* create²_μ(state) .+
#         1 .+ 2create_annihilate_μ(state)
#     ) ./ 4
# end

# # ⟨π̂ₓ⟩ = ⟨â exp(-im θ) + â† exp(im θ)⟩ / 2
# # ⟨π̂ₓ⟩ = (exp(-im θ)⟨â⟩ + exp(im θ)⟨â†⟩) / 2
# function π̂ₓ_μ(θs::AbstractVector{<:Number}, state::StateMatrix)
#     return (
#         exp.(-im*θs) .* annihilate_μ(state) .+
#         exp.(im*θs) .* create_μ(state)
#     ) ./ 2
# end
