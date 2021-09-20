export
    CoherentState,
    SqueezedState,
    ThermalState,
    SqueezedThermalState

##############
# pure state #
##############

"""
    CoherentState(α::ComplexVec{<:Real}; dim=DIM, rep=StateVector)

Coherent state is defined as the eigenstate of annihilation operator.

* `α`: Eigenvalue of annihilation operator.
* `dim`: Maximum photon number for truncate, default is $DIM.
* `rep`: In which representation, default is `StateVector`.

``\\hat{a} | \\alpha \\rangle = \\alpha | \\alpha \\rangle``

This constructor will construct ``| \\alpha \\rangle = \\hat{D}(\\alpha) | 0 \\rangle``
"""
function CoherentState(T::Type{<:Number}, α::ComplexVec; dim=DIM, rep=StateVector)
    return displace!(VacuumState(T, dim=dim, rep=rep), α)
end

CoherentState(α::ComplexVec; dim=DIM, rep=StateVector) = CoherentState(ComplexF64, α, dim=dim, rep=rep)

"""
    SqueezedState(ξ::ComplexVec{<:Real}; dim=DIM, rep=StateVector)

Squeezed state is defined if its electric field strength for some phases
has a quantum uncertainty smaller than that of a coherent state.

* `ξ`: Squeezing factor
* `dim`: Maximum photon number for truncate, default is $DIM.
* `rep`: In which representation, default is `StateVector`.

This constructor will construct ``| \\xi \\rangle = \\hat{S}(\\xi) | 0 \\rangle``
"""
function SqueezedState(T::Type{<:Number}, ξ::ComplexVec; dim=DIM, rep=StateVector)
    return squeeze!(VacuumState(T, dim=dim, rep=rep), ξ)
end

SqueezedState(ξ::ComplexVec; dim=DIM, rep=StateVector) = SqueezedState(ComplexF64, ξ, dim=dim, rep=rep)

###############
# mixed state #
###############

bose_einstein(n::Number, n̄::Real) = n̄^n / (1 + n̄)^(n+1)

bose_einstein(n̄::Real) = n -> bose_einstein(n, n̄)

"""
    ThermalState(n̄::Real; dim=DIM)

Thermal state is a mixed state with photon number distribution described by Bose-Einstein distribution.

* `n̄`: Average photon number at temperature T.
* `dim`: Maximum photon number for truncate, default is $DIM.
"""
ThermalState(T::Type{<:Number}, n̄::Real; dim=DIM) = StateMatrix(diagm(bose_einstein(n̄).(T.(0:dim-1))), dim)
ThermalState(n̄::Real; dim=DIM) = ThermalState(ComplexF64, n̄, dim=dim)

"""
    SqueezedThermalState(ξ::ComplexVec{<:Real}, n̄::Real; dim=DIM)

Squeezed state is defined if its electric field strength for some phases
has a quantum uncertainty smaller than that of a coherent state.

* `ξ`: Squeezing factor
* `n̄`: Average photon number at temperature T.
* `dim`: Maximum photon number for truncate, default is $DIM.

This constructor will construct ``\\rho = \\hat{S}(\\xi) \\rho_{th} \\hat{S}(\\xi)^{T}``
"""
function SqueezedThermalState(T::Type{<:Number}, ξ::ComplexVec, n̄::Real; dim=DIM)
    return squeeze!(ThermalState(T, n̄, dim=dim), ξ)
end

SqueezedThermalState(ξ::ComplexVec, n̄::Real; dim=DIM) = SqueezedThermalState(ComplexF64, ξ, n̄; dim=dim)
