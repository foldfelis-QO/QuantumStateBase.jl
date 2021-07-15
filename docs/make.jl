using QuantumStateBase
using Documenter

DocMeta.setdocmeta!(QuantumStateBase, :DocTestSetup, :(using QuantumStateBase); recursive=true)

withenv("LINES"=>8, "COLUMNS"=>80) do
    makedocs(;
        modules=[QuantumStateBase],
        authors="JingYu Ning <foldfelis@gmail.com> and contributors",
        repo="https://github.com/foldfelis-QO/QuantumStateBase.jl/blob/{commit}{path}#{line}",
        sitename="QuantumStateBase.jl",
        format=Documenter.HTML(;
            prettyurls=get(ENV, "CI", "false") == "true",
            canonical="https://foldfelis-qo.github.io/QuantumStateBase.jl",
            assets=String[],
        ),
        pages=[
            "Home" => "index.md",
        ],
    )
end

deploydocs(;
    repo="github.com/foldfelis-QO/QuantumStateBase.jl",
)
