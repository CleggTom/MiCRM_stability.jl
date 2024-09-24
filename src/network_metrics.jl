#network structure functions
function get_A(com)
    N,M = com.N, com.M
    A = zeros(N+M, N+M)
    
    for i = 1:N
        for j = 1:M
            A[N + j, i] = com.U[i,j]
            A[i, N + j] = com.L[i,j]
        end
    end
    
    return(A')
end

#get trophic levels
function get_TL(A)
    v = sum(A, dims = 2)[:]
    v[v .== 0] .= 1
    Λ = diagm(v) - A
    
    if rank(Λ) >= length(A[1,:])
        s = Λ \ v
        return s
    else
        return Inf
    end
end

function get_TC(A)
    s = get_TL(A)

    if isa(s, Vector)
        S = length(s)
        x = [s[i] - s[j] for i = 1:S, j = 1:S] .* A
        v = (sum(x .* x) / sum(A)) - 1
        return abs(v) < 1e-6 ? 0 : sqrt(v) 
    else
        return Inf
    end
end