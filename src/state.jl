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
function CoherentState(α::ComplexVec{<:Real}; dim=DIM, rep=StateVector)
    return displace!(VacuumState(dim=dim, rep=rep), α)
end

"""
    SqueezedState(ξ::ComplexVec{<:Real}; dim=DIM, rep=StateVector)

Squeezed state is defined if its electric field strength for some phases
has a quantum uncertainty smaller than that of a coherent state.

* `ξ`: Squeezing factor
* `dim`: Maximum photon number for truncate, default is $DIM.
* `rep`: In which representation, default is `StateVector`.

This constructor will construct ``| \\xi \\rangle = \\hat{S}(\\xi) | 0 \\rangle``
"""
function SqueezedState(ξ::ComplexVec{<:Real}; dim=DIM, rep=StateVector)
    return squeeze!(VacuumState(dim=dim, rep=rep), ξ)
end

###############
# mixed state #
###############

bose_einstein(n::Integer, n̄::Real) = n̄^n / (1 + n̄)^(n+1)

bose_einstein(n̄::Real) = n -> bose_einstein(n, n̄)

"""
    ThermalState(n̄::Real; dim=DIM)

Thermal state is a mixed state with photon number distribution described by Bose-Einstein distribution.

* `n̄`: Average photon number at temperature T.
* `dim`: Maximum photon number for truncate, default is $DIM.
"""
ThermalState(n̄::Real; dim=DIM) = StateMatrix(diagm(ComplexF64.(bose_einstein(n̄).(0:dim-1))), dim)

"""
    SqueezedThermalState(ξ::ComplexVec{<:Real}, n̄::Real; dim=DIM)

Squeezed state is defined if its electric field strength for some phases
has a quantum uncertainty smaller than that of a coherent state.

* `ξ`: Squeezing factor
* `n̄`: Average photon number at temperature T.
* `dim`: Maximum photon number for truncate, default is $DIM.

This constructor will construct ``\\rho = \\hat{S}(\\xi) \\rho_{th} \\hat{S}(\\xi)^{T}``
"""
function SqueezedThermalState(ξ::ComplexVec{<:Real}, n̄::Real; dim=DIM)
    return squeeze!(ThermalState(n̄, dim=dim), ξ)
end
