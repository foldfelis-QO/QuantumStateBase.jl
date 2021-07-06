@testset "pure state" begin
    dim = 70

    @test CoherentState(α(2., π/4), dim=dim) == displace!(VacuumState(dim=dim), α(2., π/4))
    @test CoherentState(α(2., π/4), dim=dim, rep=StateMatrix) ==
        displace!(VacuumState(dim=dim, rep=StateMatrix), α(2., π/4))
    @test SqueezedState(ξ(2., π/4), dim=dim) == squeeze!(VacuumState(dim=dim), ξ(2., π/4))
    @test SqueezedState(ξ(2., π/4), dim=dim, rep=StateMatrix) ==
        squeeze!(VacuumState(dim=dim, rep=StateMatrix), ξ(2., π/4))
end

@testset "mixed state" begin
    dim = 70

    n̄ = 0.5
    n = 5
    @test QuantumStateBase.bose_einstein(n̄)(n) == n̄^n / (1 + n̄)^(n+1)

    @test ThermalState(n̄, dim=dim) ==
        StateMatrix(diagm(QuantumStateBase.bose_einstein(n̄).(0:dim-1)), dim)
    @test SqueezedThermalState(ξ(1., π/4), n̄, dim=dim) ==
        squeeze!(ThermalState(n̄, dim=dim), ξ(1., π/4))
end
