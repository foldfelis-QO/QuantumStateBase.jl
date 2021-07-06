@testset "pdf of quadrature" begin
    @show 𝐩 = q_pdf(
        VacuumState(rep=StateMatrix),
        LinRange(0, 2π, 2),
        LinRange(-1, 1, 100),
    )
end
