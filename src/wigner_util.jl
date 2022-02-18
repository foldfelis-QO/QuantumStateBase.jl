using Mmap

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

function save_𝐰(bin_name::String, 𝐰::Array{ComplexF64,4})
    global my_artifacts
    hash = create_my_artifact() do artifact_dir
        path = joinpath(artifact_dir, bin_name)
        mem = open(path, "w+")
        write(mem, 𝐰)
        close(mem)

        return path
    end
    bind_my_artifact!(my_artifacts[], bin_name, hash)

    @info "Save Wₘₙₓₚ to $bin_name"

    return my_artifact_path(hash)
end

function load_𝐰(
    bin_path::String,
    x_range::AbstractRange,
    p_range::AbstractRange,
    dim::Integer
)
    mem = open(bin_path)
    𝐰 = Mmap.mmap(
        mem,
        Array{ComplexF64,4},
        (dim, dim, length(x_range), length(p_range))
    )
    close(mem)

    @info "Load Wₘₙₓₚ from $bin_path"

    return 𝐰
end

function gen_wigner_bin_name(
    x_range::AbstractRange,
    p_range::AbstractRange,
    dim::Integer,
)
    return string(hash((x_range, p_range, dim)))
end
