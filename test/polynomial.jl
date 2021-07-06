@testset "laguerre" begin
    @test all(
        laguerrel(n, α, x) ≈ laguerre(n, α)(x)
        for n in 0:100, α in 0:100, x in -20:0.1:20
    )
end
