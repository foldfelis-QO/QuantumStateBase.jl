@testset "factorial_ij" begin
    @test all(
        QSB.factorial_div(i, j) == factorial(big(i)) / factorial(big(j))
        for i in 0:40, j in 0:40
    )
end

@testset "z" begin
    x_range = -10:0.01:10
    p_range = -6:0.01:6
    @test all(
        QSB.z(x, p) == sqrt(2.)*(x + p*im)
        for x in x_range, p in p_range
    )

    xs = collect(x_range)
    ps = collect(p_range)
    @test QSB.z(xs, ps) == sqrt(2.).*(xs .+ ps' .* im)
end

@testset "gaussian_function" begin
    tol = 1e-14
    x_range = -10:0.01:10
    p_range = -6:0.01:6
    @test all(
        isapprox(QSB.gaussian_function(x, p), exp(-(x^2 + p^2)) / π, atol=tol)
        for x in x_range, p in p_range
    )

    xs = collect(x_range)
    ps = collect(p_range)
    @test isapprox(
        QSB.gaussian_function(xs, ps),
        exp.(-(xs.^2 .+ (ps').^2)) ./ π,
        atol=tol
    )
end

@testset "(-1)^i" begin
    @test all(
        QSB.neg_one_to_power_of(i) == (-1)^i
        for i in 1:50
    )
end

@testset "coefficient_of_wave_function" begin
    tol = 1e-14
    @test all(
        isapprox(
            QSB.coefficient_of_wave_function(m, n),
            ((n ≥ m) ? (-1)^(m-1) : (-1)^(n-1)) *
            sqrt(factorial(big(min(m, n))) / factorial(big(max(m, n)))),
            atol=tol
        )
        for m in 1:70, n in 1:70
    )
end

@testset "z_to_power" begin
    tol = 1e-14

    m_range = n_range = 1:70
    x_range = -1:0.1:1
    p_range = -0.6:0.1:0.6
    @test all(
        QSB.z_to_power(m, n, x, p) == (
            (n >= m) ? sqrt(2.)*(x - p*im) : sqrt(2.)*(x + p*im)
        )^(max(m, n) - min(m, n))
        for m in m_range, n in n_range, x in x_range, p in p_range
    )

    xs = collect(x_range)
    ps = collect(p_range)
    @test all(
        QSB.z_to_power(m, n, xs, ps) == QSB.z_to_power.(m, n, xs, ps')
        for m in m_range, n in n_range
    )
end

@testset "laguerre" begin
    m_range = n_range = 1:70
    x_range = -1:0.1:1
    p_range = -0.6:0.1:0.6

    result = Bool[]
    # n >= m
    for m in m_range, n in m:n_range.stop, x in x_range, p in p_range
        push!(
            result,
            abs(laguerre(m, n, x, p) - laguerre(m-1, n-m, abs2(sqrt(2.)*(x+p*im)))) < 1e-7
        )
    end
    # n < m
    for n in n_range, m in n:m_range.stop, x in x_range, p in p_range
        push!(
            result,
            abs(laguerre(m, n, x, p) - laguerre(n-1, m-n, abs2(sqrt(2.)*(x+p*im)))) < 1e-7
        )
    end
    @test all(result)

    xs = collect(x_range)
    ps = collect(p_range)
    @test all(
        laguerre(m, n, xs, ps) == laguerre.(m, n, xs, ps')
        for m in m_range, n in n_range
    )
end
