module QuantumStateBase
    using DataDeps
    using LinearAlgebra
    using ClassicalOrthogonalPolynomials

    const DIM = 100

    function __init__()
        register(DataDep(
            "QuantumStateBase",
            """Data for QuantumStateBase.""",
            ""
        ))

        mkpath(joinpath(DataDeps.standard_loadpath[1], "QuantumStateBase"))
    end

    datadep_root() = datadep"QuantumStateBase"

    # polynomial
    include("polynomial.jl")

    # state
    include("representation.jl")
    include("basis.jl")
    include("operator.jl")
    include("state.jl")

    # # wigner
    include("wigner_util.jl")
    include("wigner.jl")

    # pdf of quadrature
    include("quadrature_pdf.jl")
    include("sampler.jl")
end
