using QuantumStateBase
using Test
using InfiniteArrays
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

    # pdf of quadrature
    # include("quadrature_pdf.jl")
    # include("sampler.jl")
end
