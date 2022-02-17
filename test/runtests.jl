using QuantumStateBase
using Test
using LinearAlgebra
using ClassicalOrthogonalPolynomials

const QSB = QuantumStateBase
const DIM = 100

@testset "QuantumStateBase.jl" begin
    # polynomial
    include("polynomial.jl")

    # state
    include("operator.jl")
    include("state.jl")

    # wigner
    # include("wigner_util.jl")
    # include("wigner.jl")
end
