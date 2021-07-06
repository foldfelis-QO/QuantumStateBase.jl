module QuantumStateBase
    using DataDeps

    const DIM = 70

    function __init__()
        register(DataDep(
            "QuantumStateBase",
            """
            Data for QuantumStateBase.
            """,
            ""
        ))

        mkpath(joinpath(DataDeps.standard_loadpath[1], "QuantumStateBase"))
    end

    datadep_root() = datadep"QuantumStateBase"

    include("state/state.jl")
end
