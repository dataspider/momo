module MDL

setprecision(BigFloat, 64)

"""
scalable binomial coefficient
"""
function choose(n::Int64, k::Int64)
    if 0 < k <= n
        p = BigInt(1)
        for t=0:(min(k,n-k)-1)
            p = floor((p * (n - t)) / (t + 1))
        end
        return p
    else
        return 0.
    end
end

"""
binary logarithm
"""
function log2(n::Real)
    return log(2, n)
end

"""
binary logarithm returning 0. at zero
"""
function log2_zero(n::Real)
    if n == 0
        return 0.
    else
        return log2(n)
    end
end

"""
binary logarithm of binomial coefficient
"""
function log2_choose(n::Int64, k::Int64)
    if n == 0 || k == 0
        return 0.
    else
        return log2_zero(choose(n,k))
    end
end

"""
iterated logarithm (log-star)
"""
function log2_star(n::Union{Int64,Float64,BigInt,BigFloat})
    if n <= 1
        return 0.
    else
        return 1. + log2_star(log2(n))
    end
end

"""
:param z: an int >= 1
:return: the encoded size of the int according to Rissanen 1983's universal code for integers
"""
function universal_integer(z::Int64)
    NORMALIZATION_CONSTANT = 2.865064
    if z <= 0
        return 0.
    end
    c = log2(NORMALIZATION_CONSTANT)
    logstar_z = log2_star(z)
    return logstar_z + c
end

end