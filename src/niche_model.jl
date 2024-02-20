d_euclid(x,y) = sqrt( sum( (x .- y).^2 ) )
d_manhat(x,y) = sum( abs.(x .- y))

"""
    generate_L_mat(c::Vector{Vector{Float64}},fd::Function,mp::Float64)

Takes a vector of niche values and generates a leakage matrix with distance function `fd`.
"""
function generate_L_mat(M::Int64, ca::Vector{Vector{Float64}},fd::Function,p::Float64)
    #calculate distance matrix
    L = zeros(M,M)
    d = 0.2
    for i = 1:M
        for j = 1:M
            if i == j
                # dist[i,j] = Inf
            else
                if sum(ca[i]) < sum(ca[j])
                    L[i,j] = fd(ca[i],ca[j]) < p ? rand() : 0
                end
            end
        end
       
        if sum(L[i,:]) > 0
            L[i,:] .= L[i,:] ./ sum(L[i,:])
        end
    end

    # #get probablilty
    # d_filter = filter(x -> !isinf(x), dist[:])
    # l = -log(mp) / mean(d_filter)
    # p = exp.(-l * dist)

    # #return prob
    # return(rand(M,M) .< p)
    return(L)
end

"""
    generate_U_mat(N::Int64, M::Int64, ci::Vector{Vector{Float64}},ca::Vector{Vector{Float64}},fd::Function)

Calculate uptake matrix from the two niche vectors `ci` and `ca`.
"""
function generate_U_mat(N::Int64, M::Int64, ci::Vector{Vector{Float64}},ca::Vector{Vector{Float64}},fd::Function, β::Float64)
    #calculate distance matrix
    U = fill(0.0,N,M)
    for i = 1:N
        for a = 1:M
            U[i,a] = fd(ci[i],ca[a]) < r[i] ? rand() : 0
        end

    end

    return(U)

end