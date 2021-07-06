export
    Creation,
    create!,
    create,
    Annihilation,
    annihilate!,
    annihilate,

    Arg,
    α,
    ξ,

    Displacement,
    displace!,

    Squeezing,
    squeeze!

############
# a† and a #
############

Creation(; dim=DIM) = diagm(-1 => sqrt.(1:dim-1))

function create!(state::StateVector{<:Number})
    dim = state.dim
    state.v = Creation(dim=dim) * state.v

    return state
end

function create!(state::StateMatrix{<:Number})
    dim = state.dim
    𝐜 = Creation(dim=dim)
    state.𝛒 = 𝐜 * state.𝛒 * 𝐜'

    return state
end

create(state::AbstractState) = create!(copy(state))

Annihilation(; dim=DIM) = diagm(1 => sqrt.(1:dim-1))

function annihilate!(state::StateVector{<:Number})
    dim = state.dim
    state.v = Annihilation(dim=dim) * state.v

    return state
end

function annihilate!(state::StateMatrix{<:Number})
    dim = state.dim
    𝐚 = Annihilation(dim=dim)
    state.𝛒 = 𝐚 * state.𝛒 * 𝐚'

    return state
end

annihilate(state::AbstractState) = annihilate!(copy(state))

###########
# α and ξ #
###########

struct Arg{T <: Real}
    r::T
    θ::T
end

Base.show(io::IO, arg::Arg{T}) where {T} = print(io, "Arg{$T}($(arg.r)exp($(arg.θ)im))")

z(arg::Arg{<:Real}) = arg.r * exp(im * arg.θ)

α(r::T, θ::T) where {T} = Arg{T}(r, θ)
const ξ = α

################
# displacement #
################

function Displacement(α::Arg{<:Real}; dim=DIM)
    return exp(z(α) * Creation(dim=dim) - z(α)' * Annihilation(dim=dim))
end

function displace!(state::StateVector{<:Number}, α::Arg{<:Real})
    dim = state.dim
    state.v = Displacement(α, dim=dim) * state.v

    return state
end

function displace!(state::StateMatrix{<:Number}, α::Arg{<:Real})
    dim = state.dim
    𝐝 = Displacement(α, dim=dim)
    state.𝛒 = 𝐝 * state.𝛒 * 𝐝'

    return state
end

#############
# squeezing #
#############

function Squeezing(ξ::Arg{<:Real}; dim=DIM)
    return exp(0.5 * z(ξ)' * Annihilation(dim=dim)^2 - 0.5 * z(ξ) * Creation(dim=dim)^2)
end

function squeeze!(state::StateVector{<:Number}, ξ::Arg{<:Real})
    dim = state.dim
    state.v = Squeezing(ξ, dim=dim) * state.v

    return state
end

function squeeze!(state::StateMatrix{<:Number}, ξ::Arg{<:Real})
    dim = state.dim
    𝐬 = Squeezing(ξ, dim=dim)
    state.𝛒 = 𝐬 * state.𝛒 * 𝐬'

    return state
end

###############
# measurement #
###############

# ##### for arb. statein θ-x quadrature coordinate #####

# |θ, x⟩ = ∑ₙ |n⟩ ⟨n|θ, x⟩ = ∑ₙ ψₙ(θ, x) |n⟩
# ⟨n|θ, x⟩ = ψₙ(θ, x) = exp(im n θ) (2/π)^(1/4) exp(-x^2) Hₙ(√2 x)/√(2^n n!)
# coeff_ψₓ = (2/π)^(1/4)/√(2^n n!)
# ψₓ = coeff_ψₓ(n) exp(im n θ) exp(-x^2) Hₙ(√2 x)
calc_coeff_ψₙ(n::BigInt) = (n/π)^(1/4) / sqrt(2^n * factorial(n))
COEFF_ψₓ = [calc_coeff_ψₙ(big(n)) for n in 0:499]

function coeff_ψₙ(n::Integer)
    (n < 500) && (return COEFF_ψₓ[n+1])
    return calc_coeff_ψₙ(n)
end

function ψₙ(n::Integer, θ::Real, x::Real)
    return coeff_ψₙ(n) * exp(im * n * θ - x^2) * hermiteh(n, sqrt(2)x)
end

function 𝛑̂!(result::Matrix{<:Complex}, θ::Real, x::Real; dim=DIM)
    view(result, :, 1) .= ψₙ.(0:dim-1, θ, x)
    result .= view(result, :, 1) * view(result, :, 1)'

    return result
end

function 𝛑̂(θ::Real, x::Real; dim=DIM, T=ComplexF64)
    result = Matrix{T}(undef, dim, dim)
    U = T.parameters[1]

    return 𝛑̂!(result, U(θ), U(x), dim=dim)
end

# ##### for Gaussian state in θ-x quadrature coordinate #####

# π̂ₓ = (â exp(-im θ) + â† exp(im θ)) / 2

tr_mul(𝐚, 𝐛) = sum(𝐚[i, :]' * 𝐛[:, i] for i in 1:size(𝐚, 1))
create_μ(state::StateMatrix) = tr_mul(Creation(dim=state.dim), state.𝛒)
create²_μ(state::StateMatrix) = tr_mul(Creation(dim=state.dim)^2, state.𝛒)
annihilate_μ(state::StateMatrix) = tr_mul(Annihilation(dim=state.dim), state.𝛒)
annihilate²_μ(state::StateMatrix) = tr_mul(Annihilation(dim=state.dim)^2, state.𝛒)
create_annihilate_μ(state::StateMatrix) = tr_mul(
    Creation(dim=state.dim) * Annihilation(dim=state.dim),
    state.𝛒
)

# ⟨π̂ₓ²⟩ = ⟨ââ exp(-2im θ) + â†â† exp(2im θ) + ââ† + â†â⟩ / 4
# ⟨π̂ₓ²⟩ = (exp(-2im θ)⟨â²⟩ + exp(2im θ)⟨â†²⟩ + 1 + 2⟨ââ†⟩) / 4
# here, ⟨ââ† + â†â⟩ = 1 + 2⟨ââ†⟩ due to the commutation relation
function π̂ₓ²_μ(θs::AbstractVector{<:Number}, state::StateMatrix)
    return (
        exp.(-2im*θs) .* annihilate²_μ(state) .+
        exp.(2im*θs) .* create²_μ(state) .+
        1 .+ 2create_annihilate_μ(state)
    ) ./ 4
end

# ⟨π̂ₓ⟩ = ⟨â exp(-im θ) + â† exp(im θ)⟩ / 2
# ⟨π̂ₓ⟩ = (exp(-im θ)⟨â⟩ + exp(im θ)⟨â†⟩) / 2
function π̂ₓ_μ(θs::AbstractVector{<:Number}, state::StateMatrix)
    return (
        exp.(-im*θs) .* annihilate_μ(state) .+
        exp.(im*θs) .* create_μ(state)
    ) ./ 2
end
