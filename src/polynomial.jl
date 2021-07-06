export
    laguerre

function laguerre(n::Integer, α::Integer, x::T) where {T<:Real}
    if n == 0
        return one(T)
    elseif n == 1
        return one(T) + α - x
    else
        L_prev = one(T)
        L_k = one(T) + α - x
        for k = 1:(n-1)
            L_next = ((2k + one(T) + α - x) * L_k - (k+α) * L_prev) / (k+one(T))
            L_prev, L_k = L_k, L_next
        end
        return L_k
    end
end

laguerre(n::Integer, α::Integer) = x -> laguerre(n, α, x)

# function hermite(n::T, x::Real) where {T <: Integer}
#     result = 0
#     for k in 0:floor(T, n/2)
#         coeff = (-1) ^ k *
#             2 ^ (n-2k) *
#             factorial(n) / (factorial(k)*factorial(n-2k))

#         x_power = x^(n-2k)

#         result += coeff * x_power
#     end

#     return result
# end
