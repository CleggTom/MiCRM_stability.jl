
"""
    random_community(N,M,Cu,Cd)

Randomly generated communtiy with uptake and stoichometry matricies createcd with random graph structure and connectance 'Cu' and 'Cd'
"""
function random_community(N,M,Cu,Cd)
    U = Float64.(rand(N,M) .< Cu)
    # U .= U .* rand(N,M)
    U[U .< 1e-5] .= 0.0
    #normalise consumer uptakes to 1
    [U[i,:] .= U[i,:] ./ sum(U[i,:]) for i = 1:N]
    U[isnan.(U)] .= 0.0

    D =  Float64.(rand(M,M) .< Cd) 
    # D .= D .* rand(M,M)
    D[D .< 1e-5] .= 0.0
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
    D[isnan.(D)] .= 0.0

    return(Community(N,M,U,D))
end


"""
    rewire_row!(row::BitVector, σ::Real)

Randomly shuffles a row with probabiltiy σ
"""
function rewire_row!(row::BitVector, σ::Real)
    σ <= 0.0 && return
    for idx in findall(row)
        if rand() < σ
            avail = findall(.!row)
            if !isempty(avail)
                row[idx] = false
                row[rand(avail)] = true
            end
        end
    end
end

"""
    interpolated_niche_model(N, M, C, d, p)

Interpolates between a structured 1D niche model (p = 1.0) 
and a randomized consumer-resource network (p = 0.0).
"""
function niche_model(N::Int, M::Int, C::Real, d::Real, σ::Real; return_niche = false)

    ni = sort(rand(N)) # Consumers
    na = sort(rand(M)) # Resources

    #distributions for connectivity
    dist_C = Beta(1.0, (1 / C) - 1)
    dist_d = Beta(1.0, (1 / d) - 1)

    # 1. Consumer uptake matrix U (N x M)
    U = zeros(Float64, N, M)
    for i in 1:N
        row = falses(M)
        while !any(row)
            wi = rand(dist_C) / 2
            row .= (na .> (ni[i] - wi)) .&& (na .< (ni[i] + wi))
        end
        rewire_row!(row, σ)
        if sum(row) > 0
            U[i, :] .= row ./ sum(row)
        end
    end

    # 2. Resource leakage matrix L (M x M)
    L = zeros(Float64, M, M)
    for α in 1:(M - 1)
        row = falses(M)
        while !any(row)
            #draw position
            ca = rand(dist_d)
            pa = na[α] + (1.0 - na[α]) * ca
            #draw width
            wa = rand(dist_C) / 2
            
            # Niche leakage strictly downstream: na[α] < target
            row .= (na .>  max(na[α], pa - wa)) .&& (na .< (pa + wa))

        end
        rewire_row!(row, σ)
        row[α] = 0.0
        if sum(row) > 0
            L[α, :] .= row ./ sum(row)
        end
    end
    if !return_niche
        return MiCRM_stability.Community(N, M, U, L)
    else
        return MiCRM_stability.Community(N, M, U, L), (ni,na)
    end 
end