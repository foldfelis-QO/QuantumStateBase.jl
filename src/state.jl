export
    # Fock basis
    FockState,
    NumberState,
    VacuumState,
    SinglePhotonState,

    # pure state
    CoherentState,
    SqueezedState,

    # mixed state
    ThermalState,
    SqueezedThermalState

################
# constructors #
################

##### Fock basis #####

"""
    FockState([T::Type{<:Number}], n::Integer, [rep::Type]; dim::Integer)

Construct a Fock state in `rep`  representation.

## Arguments

* `T`: Numerical type. Default as `Float64` if not assigned.
* `n`: Photon number of Fock state.
* `rep`: In which representation, default as `Vector` if not assigned.
* `dim`: Maximum photon number for truncate.
"""
function FockState(T::Type{<:Number}, n::Integer, ::Type{Vector}; dim)
    v = zeros(T, dim); v[n+1] = one(T)

    return v
end

function FockState(T::Type{<:Number}, n::Integer, ::Type{Matrix}; dim)
    ρ = zeros(T, dim, dim); ρ[n+1, n+1] = one(T)

    return ρ
end

FockState(n, rep; dim) = FockState(Float64, n, rep, dim=dim)

FockState(n; dim) = FockState(n, Vector, dim=dim)

const NumberState = FockState

"""
    VacuumState([T::Type{<:Number}], [rep::Type]; dim::Integer)

Construct a vacuum state in `rep` representation.

## Arguments

* `T`: Numerical type. Default as `Float64` if not assigned.
* `rep`:In which representation, default as `Vector` if not assigned.
* `dim`: Maximum photon number for truncate.
"""
VacuumState(T, rep; dim) = FockState(T, 0, rep, dim=dim)

VacuumState(rep; dim) = FockState(Float64, 0, rep, dim=dim)

VacuumState(; dim) = FockState(0, Vector, dim=dim)

"""
    SinglePhotonState([T::Type{<:Number}], [rep::Type]; dim::Integer)

Construct a single photon state in `rep` representation.

## Arguments

* `T`: Numerical type. Default as `Float64` if not assigned.
* `rep`:In which representation, default as `Vector` if not assigned.
* `dim`: Maximum photon number for truncate.
"""
SinglePhotonState(T, rep; dim) = FockState(T, 1, rep, dim=dim)

SinglePhotonState(rep; dim) = FockState(Float64, 1, rep, dim=dim)

SinglePhotonState(; dim) = FockState(1, Vector, dim=dim)

##### pure state #####

### coherent state

"""
    CoherentState([T::Type{<:Complex}], r::Real, θ::Real, [rep::Type]; dim::Integer)

Construct a coherent state in `rep` representation.

* `T`: Numerical type. Default as `ComplexF64` if not assigned.
* `r`: The radius in phase space.
* `θ`: The phase angle in phase space.
* `rep`:In which representation, default as `Vector` if not assigned.
* `dim`: Maximum photon number for truncate.
"""
function CoherentState(T::Type{<:Complex}, r, θ, rep; dim)
    return displace(VacuumState(real(T), rep, dim=dim), r, θ)
end

CoherentState(r, θ, rep; dim) = CoherentState(ComplexF64, r, θ, rep, dim=dim)

CoherentState(r, θ; dim) = CoherentState(r, θ, Vector, dim=dim)

### squeezed state

"""
    SqueezedState([T::Type{<:Complex}], r::Real, θ::Real, [rep::Type]; dim::Integer)

Construct a squeezed state in `rep` representation.

Squeezed state is defined if its electric field strength for some phases has a quantum
uncertainty smaller than that of a coherent state.

## Arguments

* `T`: Numerical type. Default as `ComplexF64` if not assigned.
* `r`: The squeezing level.
* `θ`: The squeezing phase.
* `rep`:In which representation, default as `Vector` if not assigned.
* `dim`: Maximum photon number for truncate.

This constructor will construct ``| \\xi \\rangle = \\hat{S}(\\xi) | 0 \\rangle``
"""
function SqueezedState(T::Type{<:Complex}, r, θ, rep; dim)
    return squeeze(VacuumState(real(T), rep, dim=dim), r, θ)
end

SqueezedState(r, θ, rep; dim) = SqueezedState(ComplexF64, r, θ, rep, dim=dim)

SqueezedState(r, θ; dim) = SqueezedState(r, θ, Vector, dim=dim)

##### mixed state #####

### thermal state

bose_einstein(n::Number, n̄::Real) = n̄^n / (1 + n̄)^(n+1)

bose_einstein(n̄) = n -> bose_einstein(n, n̄)

"""
    ThermalState([T::Type{<:Number}], n̄::Real; dim::Integer)

Construct a thermal state.

Thermal state is a mixed state with photon number distribution described by
Bose-Einstein distribution.

## Arguments

* `T`: Numerical type. Default as `Float64` if not assigned.
* `n̄`: Average photon number at a certain temperature.
* `dim`: Maximum photon number for truncate.
"""
ThermalState(T::Type{<:Number}, n̄; dim) = diagm(bose_einstein(n̄).(T.(0:dim)))

ThermalState(n̄; dim) = ThermalState(Float64, n̄, dim=dim)

### squeezed thermal state

"""
    SqueezedThermalState([T::Type{<:Complex}], r::Real, θ::Real, n̄::Real; dim::Integer)

Squeezed thermal state is defined if its electric field strength for some phases has a
quantum uncertainty smaller than that of a coherent state.

## Arguments

* `T`: Numeric type. Default as `ComplexF64` if not assigned.
* `r`: Squeezing level.
* `θ`: Squeezing phase.
* `n̄`: Average photon number at a certain temperature.
* `dim`: Maximum photon number for truncate.

This constructor will construct ``\\rho = \\hat{S}(\\xi) \\rho_{th} \\hat{S}(\\xi)^{T}``
"""
function SqueezedThermalState(T::Type{<:Complex}, r, θ, n̄; dim)
    return squeeze(ThermalState(real(T), n̄, dim=dim), r, θ)
end

SqueezedThermalState(r, θ, n̄; dim) = SqueezedThermalState(ComplexF64, r, θ, n̄, dim=dim)

###########
# methods #
###########

# TODO: ρ(v)
# TODO: purity(ρ)
# TODO: fidelity(ρ1, ρ2)
# TODO: fidelity(v1, v2) ???
