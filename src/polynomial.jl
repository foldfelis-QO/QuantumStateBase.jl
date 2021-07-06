using ClassicalOrthogonalPolynomials

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
