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

function modular_matrix(N, M, K, Q, total_links)
    mat = zeros(N,M)
    
    # 1. Define the "Inside" and "Outside" budgets
    internal_budget = round(Int, total_links * Q)
    external_budget = total_links - internal_budget
    
    # 2. Identify which coordinates are "Internal"
    rows_per_mod = div(N, K)
    cols_per_mod = div(M, K)
    
    internal_coords = Vector{CartesianIndex{2}}()
    external_coords = Vector{CartesianIndex{2}}()
    
    for r in 1:N, c in 1:M
        # Determine which module this row and column belong to
        mod_r = min(K, div(r-1, rows_per_mod) + 1)
        mod_c = min(K, div(c-1, cols_per_mod) + 1)
        
        if mod_r == mod_c
            push!(internal_coords, CartesianIndex(r, c))
        end
        
        push!(external_coords, CartesianIndex(r, c))
    end

    # 3. Sample from the budgets
    # We use 'sample' from StatsBase or just shuffle
    
    for idx in shuffle(internal_coords)[1:min(length(internal_coords), internal_budget)]
        mat[idx] = 1.0
    end
    
    for idx in shuffle(external_coords)[1:min(length(external_coords), external_budget)]
        mat[idx] = 1.0
    end

    return mat
end

function modular_community(N,M,Ku,Kd,Qu,Qd,Lu,Ld)
    U = MiCRM_stability.modular_matrix(N,M,Ku,Qu,Lu)
    [U[i,:] .= U[i,:] ./ sum(U[i,:]) for i = 1:N]
    U[isnan.(U)] .= 0.0
    
    D = MiCRM_stability.modular_matrix(M,M,Kd,Qd,Ld)
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
    D[isnan.(D)] .= 0.0

    return(MiCRM_stability.Community(N,M,U,D))
end

function hierarchical_modular_matrix(N, M, K, Q, total_links, off)
    mat = zeros(N, M)
    
    # 1. Define budgets
    internal_budget = round(Int, total_links * Q)
    external_budget = total_links - internal_budget
    
    rows_per_mod = div(N, K)
    cols_per_mod = div(M, K)
    
    # 2. Identify coordinates
    # "Internal" now means mod_c <= mod_r (Lower triangular block structure)
    internal_coords = Vector{CartesianIndex{2}}()
    external_coords = Vector{CartesianIndex{2}}()
    
    for r in 1:N, c in 1:M
        mod_r = min(K, div(r-1, rows_per_mod) + 1)
        mod_c = min(K, div(c-1, cols_per_mod) + 1)
        
        # This is the change: links are "Internal" if they point to 
        # a module at or below the current one (Trophic Flow)
        if any(mod_c .== (mod_r .+ off)) 
            push!(internal_coords, CartesianIndex(r, c))
        else
            push!(external_coords, CartesianIndex(r, c))
        end
    end

    # 3. Sample from the budgets
    # Internal: Filling the lower-triangular blocks
    for idx in shuffle(internal_coords)[1:min(length(internal_coords), internal_budget)]
        mat[idx] = 1.0 # Or use your specific distribution for weights
    end
    
    # External: Filling the upper-triangular blocks (the "Long loops")
    for idx in shuffle(external_coords)[1:min(length(external_coords), external_budget)]
        mat[idx] = 1.0
    end

    return mat
end

function hierarchical_community(N,M,Ku,Kd,Qu,Qd,Lu,Ld)
    U = MiCRM_stability.modular_matrix(N,M,Ku,Qu,Lu)
    [U[i,:] .= U[i,:] ./ sum(U[i,:]) for i = 1:N]
    U[isnan.(U)] .= 0.0
    
    D = MiCRM_stability.hierarchical_modular_matrix(M,M,M,Qd,Ld,1:Kd)
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
    D[isnan.(D)] .= 0.0

    return(MiCRM_stability.Community(N,M,U,D))
end