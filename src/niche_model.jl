d_euclid(x,y) = sqrt( sum( (x .- y).^2 ) )
d_manhat(x,y) = sum( abs.(x .- y))

function rand_community(N,M,C,Nin, Nout)
    U = rand(N,M) .* (rand(N,M) .< C)
    L = rand(N,M) .* (rand(N,M) .< C)

    U[:,rand(1:M, Nout)] .= 0
    L[:,rand(1:M, Nin)] .= 0

    return( Community(N,M,U,L))
end

"""
    rand_com(N,M,C, ci = rand(N), ca = rand(M))

Niche model
"""
function niche_community(N,M,C, ci = rand(N), ca = rand(M))
    ci = ci[sortperm(sum.(ci))]    
    ca = ca[sortperm(sum.(ca))]
    d = ci .- ca'
    
    fu(x, dist) = dist > x >= 0 ? rand() : 0.0
    fl(x, dist) = -dist < x < 0 ? rand() : 0.0
    
    U = zeros(N,M)
    L = zeros(N,M)

    b = (1/C - 1)
    dβ = Beta(1, b)

    for i = 1:N
        # ind = findmin(abs.(d[i,:]))[2]
        x = d[i,:] #.- d[i,ind]
        U[i,:] .= fu.(x, rand(dβ))
        L[i,:] .= fl.(x, rand(dβ))
    end

    #end products    
    U[:,end] .= 0.0
    L[:,1] .= 0.0

    return(MiCRM_stability.Community(N,M,U,L))
end


function modular_community(N,M,Nclass,pgen,Cg,Co,Ci,Csp,Csw,Css)
    res_class = repeat(collect(1:Nclass), inner = M ÷ Nclass)
    rclass = Dict(i => findall(res_class .== i) for i = 1:Nclass)

    #Consumers
    U = zeros(N,M)
    con_class = zeros(N)
    for i = 1:N
        if rand() < pgen
            con_class[i] = 0
            U[i,:] .= rand(M) .< Cg
        else
            con_class[i] = rand(1:(Nclass-1))
            U[i,:] .= rand(M) .< Co
            U[i, rclass[con_class[i]],:] .= rand(length(rclass[con_class[i]])) .< Ci
        end
    end
    U[:, rclass[Nclass]] .= 0
    U = U[sortperm(con_class),:]

    D = zeros(M,M)
    for a = 1:M
        #all leakage
        D[a, :] .= rand(M) .< Csp
        #leakage to waste
        D[a, rclass[Nclass]] .= rand(length(rclass[Nclass])) .< Csw
        #leakage to same group
        D[a, rclass[res_class[a]]] .= rand(length(rclass[res_class[a]])) .< Css
    end

    D[:, rclass[1]] .= 0

    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]

    return(Community(N,M,rand(N,M) .* U,rand(N,M) .* (U * D)))

end

function nested_community(N,M,κ,v)
    U = zeros(N,M)
    L = zeros(N,M)

    p = (1.0:M) .^ κ
    p = p ./ sum(p)
    
    kU = sample(1:M, Weights(p), N) |> sort |> reverse
    kL = sample(1:M, Weights(p), N) |> sort |> reverse
    
    for i = 1:N
        Z = kU[i]
        nj = (sum(U,dims=1)[:] .+ 1) .^ v
        p = nj ./ sum(nj)
        U[i,:] .= rand(Multinomial(Z,p)) .> 0

        Z = kL[i]
        nj = (sum(L,dims=1)[:] .+ 1) .^ v
        p = nj ./ sum(nj)
        L[i,:] .= rand(Multinomial(Z,p)) .> 0
    end
    ind = sortperm(sum(U, dims = 1)[:])
    return(Community(N,M,rand(N,M) .* U[:,ind], rand(N,M) .* L[:,ind]))
end

function empirical_communty(N, C, D)
    ind = rand(1:size(C)[1], N)
    return(Community(N,size(C)[2],C[ind,:],D[ind,:]))
end
