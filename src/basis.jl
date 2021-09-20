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

"""
    FockState(n; T=ComplexF64, dim=DIM, rep=StateVector)

Construct a Fock state in `rep` representation.

* `n`: Photon number of Fock state.
* `T`: Numeric data types, default is `ComplexF64`
* `dim`: Maximum photon number for truncate, default is $DIM
* `rep`: In which representation, default is `StateVector`.

# Examples
```jldoctest
julia> state = FockState(1);

julia> state = FockState(1, dim=100);

julia> state = FockState(1, rep=StateMatrix);
```
"""
FockState(T::Type{<:Number}, n; dim=DIM, rep=StateVector) = FockState(T, n, dim, rep)
FockState(n; dim=DIM, rep=StateVector) = FockState(ComplexF64, n, dim, rep)

"""
    NumberState(n; T=ComplexF64, dim=DIM, rep=StateVector)

Exactly the same with `FockState`.
"""
NumberState(T::Type{<:Number}, n; dim=DIM, rep=StateVector) = FockState(T, n, dim, rep)
NumberState(n; dim=DIM, rep=StateVector) = FockState(ComplexF64, n, dim, rep)

"""
    VacuumState(; T=ComplexF64, dim=DIM, rep=StateVector)

Construct a vacuum state in `rep` representation.

* `T`: Numeric data types, default is `ComplexF64`
* `dim`: Maximum photon number for truncate, default is $DIM
* `rep`: In which representation, default is `StateVector`.

# Examples
```jldoctest
julia> state = VacuumState();

julia> state = VacuumState(dim=100);

julia> state = VacuumState(rep=StateMatrix);
```
"""
VacuumState(T::Type{<:Number}; dim=DIM, rep=StateVector) = FockState(T, 0, dim, rep)
VacuumState(; dim=DIM, rep=StateVector) = FockState(ComplexF64, 0, dim, rep)

"""
    SinglePhotonState(; T=ComplexF64, dim=DIM, rep=StateVector)

Construct a single photon state in `rep` representation.

* `T`: Numeric data types, default is `ComplexF64`
* `dim`: Maximum photon number for truncate, default is $DIM
* `rep`: In which representation, default is `StateVector`.

# Examples
```jldoctest
julia> state = SinglePhotonState();

julia> state = SinglePhotonState(dim=100);

julia> state = SinglePhotonState(rep=StateMatrix);
```
"""
SinglePhotonState(T::Type{<:Number}; dim=DIM, rep=StateVector) = FockState(T, 1, dim, rep)
SinglePhotonState(; dim=DIM, rep=StateVector) = FockState(ComplexF64, 1, dim, rep)
