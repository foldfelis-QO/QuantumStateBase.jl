########################
# fractions of formula #
########################

factorial_div(i::Integer, j::Integer) = factorial(big(i)) / factorial(big(j))

z(x::Real, p::Real) = sqrt(2.)*(x + p*im)
z(x::Vector{<:Real}, p::Vector{<:Real}) = z.(x, p')

gaussian_function(x::Real, p::Real) = exp(-0.5 * abs2(z(x,p))) / π
gaussian_function(x::Vector{<:Real}, p::Vector{<:Real}) = gaussian_function.(x, p')

neg_one_to_power_of(i::Integer) = (i % 2 == 0) ? 1 : -1

function coefficient_of_wave_function(m::Integer, n::Integer)
    if n ≥ m
        # adjust index bases for number state on `m`
        return neg_one_to_power_of(m-1) * sqrt(factorial_div(m, n))
    else
        # adjust index bases for number state on `n`
        return neg_one_to_power_of(n-1) * sqrt(factorial_div(n, m))
    end
end

function z_to_power(m::Integer, n::Integer, x::Real, p::Real)
    if n ≥ m
        return conj(z(x, p'))^(n - m)
    else
        return z(x, p')^(m - n)
    end
end

function z_to_power(m::Integer, n::Integer, x::Vector{<:Real}, p::Vector{<:Real})
    z_to_power.(m, n, x, p')
end

function laguerre(m::Integer, n::Integer, x::Real, p::Real)
    if n ≥ m
        # adjust index bases for number state on `m`
        return laguerre(m-1, n - m, abs2(z(x, p)))
    else
        # adjust index bases for number state on `n`
        return laguerre(n-1, m - n, abs2(z(x, p)))
    end
end

function laguerre(m::Integer, n::Integer, x::Vector{<:Real}, p::Vector{<:Real})
    return laguerre.(m, n, x, p')
end

#########
# utils #
#########

function save_𝐰(bin_path::String, 𝐰::Array{ComplexF64,4})
    @info "Save W_{m,n,x,p} to $bin_path"
    mem = open(bin_path, "w+")
    write(mem, 𝐰)
    close(mem)
end

function load_𝐰(
    m_dim::Integer,
    n_dim::Integer,
    x_range::AbstractRange,
    p_range::AbstractRange,
    bin_path::String
)
    @info "Load W_{m,n,x,p} from $bin_path"
    mem = open(bin_path)
    𝐰 = Mmap.mmap(
        mem,
        Array{ComplexF64,4},
        (m_dim, n_dim, length(x_range), length(p_range))
    )
    close(mem)

    return 𝐰
end

function gen_wigner_bin_path(
    m_dim::Integer,
    n_dim::Integer,
    x_range::AbstractRange,
    p_range::AbstractRange,
)
    bin_path = joinpath(
        mkpath(joinpath(datadep_root(), "wigner_function")),
        "W " *
        "m=$(m_dim) n=$(n_dim) " *
        "x=$(range2str(x_range)) p=$(range2str(p_range)).bin"
    )

    return bin_path
end

range2str(range::AbstractRange) = replace(string(range), r":|," => "_")

check_zero(m_dim, n_dim) = !iszero(m_dim) && !iszero(n_dim)

check_empty(x_range, p_range) = !isempty(x_range) && !isempty(p_range)

function check_argv(m_dim, n_dim, x_range, p_range)
    return check_zero(m_dim, n_dim) && check_empty(x_range, p_range)
end
