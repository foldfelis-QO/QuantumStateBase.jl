export
    FockState,
    NumberState,
    VacuumState,
    SinglePhotonState

function FockState(T::Type{<:Number}, n::Integer, dim::Integer, ::Type{StateVector})
    v = zeros(T, dim)
    v[n+1] = 1

    return StateVector{T}(v, dim)
end

function FockState(T::Type{<:Number}, n::Integer, dim::Integer, ::Type{StateMatrix})
    𝛒 = zeros(T, dim, dim)
    𝛒[n+1, n+1] = 1

    return StateMatrix{T}(𝛒, dim)
end

FockState(n; T=ComplexF64, dim=DIM, rep=StateVector) = FockState(T, n, dim, rep)

NumberState(n; T=ComplexF64, dim=DIM, rep=StateVector) = FockState(T, n, dim, rep)

VacuumState(; T=ComplexF64, dim=DIM, rep=StateVector) = FockState(T, 0, dim, rep)

SinglePhotonState(; T=ComplexF64, dim=DIM, rep=StateVector) = FockState(T, 1, dim, rep)
