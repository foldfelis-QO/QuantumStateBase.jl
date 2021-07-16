export
    CoherentState,
    SqueezedState,
    ThermalState,
    SqueezedThermalState

##############
# pure state #
##############

function CoherentState(α::ComplexVec{<:Real}; dim=DIM, rep=StateVector)
    return displace!(VacuumState(dim=dim, rep=rep), α)
end

function SqueezedState(ξ::ComplexVec{<:Real}; dim=DIM, rep=StateVector)
    return squeeze!(VacuumState(dim=dim, rep=rep), ξ)
end

###############
# mixed state #
###############

bose_einstein(n::Integer, n̄::Real) = n̄^n / (1 + n̄)^(n+1)

bose_einstein(n̄::Real) = n -> bose_einstein(n, n̄)

ThermalState(n̄::Real; dim=DIM) = StateMatrix(diagm(ComplexF64.(bose_einstein(n̄).(0:dim-1))), dim)

function SqueezedThermalState(ξ::ComplexVec{<:Real}, n̄::Real; dim=DIM)
    return squeeze!(ThermalState(n̄, dim=dim), ξ)
end
