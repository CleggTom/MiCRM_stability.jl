function random_community(N,M,Cu,Cd)
    U = Float64.(rand(N,M) .< Cu) 
    #evey consumer has a link
    for i = 1:N
        if sum(U[i,:]) == 0
            U[i, rand(1:M)] = 1
        end
    end
    D =  Float64.(rand(M,M) .< Cd) 
    
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
    D[isnan.(D)] .= 0.0

    return(Community(N,M,U,D))
end



"""
    nested_community(N,M,κ,v, Cd)

Creates a nested `Communtiy` object with N consumers and M resources. The leakage matrix is random. `κ` is the exponent controling the distribution of node degrees and `v` the strength of attraction to high degree nodes. in this way `v` controls the nestedness. 
"""
function nested_community(N,M, k, v, Cd)
    U = zeros(N,M)

    p = (1:M) .^ -k
    p = p ./ sum(p)

    #sample numbers of resources and sort
    kU = sample(1:M, Weights(p), N, replace = true) |> sort |> reverse
    #loop over consumers
    for i = 1:N
        #number of links
        Z = kU[i]
        #get prob based on other consumption
        #v - controls the strength of attraction
        nj = (sum(U,dims=1)[:] .+ 1) .^ v
        p = nj ./ sum(nj)
        #sample multinomial
        U[i,sample(1:M, Weights(p), Z, replace = false)] .= 1.0     
    end
    #order by number of links
    ind = sortperm(sum(U .> 0, dims = 1)[:], rev = true)

    # U[:,rand(1:M,rand(1:M))] .= 0

    #random uptake matrix
    D =  Float64.(rand(M,M) .< Cd) 
    # D[:,rand(1:M,rand(1:M))] .= 0 #some resources not leaked
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
    D[isnan.(D)] .= 0.0

    return(MiCRM_stability.Community(N,M,U[:,ind],D))
end

"""
    modular_community(N,M,C,Nclass,r)

Creates a modular communtiy 
"""
function modular_community(N,M,C,Nclass,r)
    Mc = M ÷ Nclass
    Ci = (C*M) / (Mc + r*(M - Mc))
    Co = r*Ci

    res_class = repeat(collect(1:Nclass), inner = M ÷ Nclass)
    rclass = Dict(i => findall(res_class .== i) for i = 1:Nclass)

    #Consumers
    U = zeros(N,M)
    con_class = zeros(N)
    for i = 1:N
            con_class[i] = rand(1:(Nclass))
            U[i,:] .= rand(M) .< Co
            U[i, rclass[con_class[i]],:] .= rand(length(rclass[con_class[i]])) .< Ci
        # end
    end

    U = U[sortperm(con_class),:]

    D = zeros(M,M)
    for a = 1:M
        #all leakage
        D[a, :] .= rand(M) .< Co
        #leakage to waste
        #leakage to same group
        D[a, rclass[res_class[a]]] .= rand(length(rclass[res_class[a]])) .< Ci
    end

    D = D[sortperm(res_class),:]
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]

    return(Community(N,M, U, D))

end


# ##Trophic coherence
function coherence_D(M,λd,h,T)
    #create skeleton
    TL = fill(-Inf,M)
    TL[1] = 1
    
    #add iteratively to D-matrix
    D = zeros(M,M)
    for i = 2:M
        ind = 1:(i-1)
        #choose base
        pb = (TL[ind] / maximum(TL[ind])) .^ h
        b = sample(ind, Weights(pb[ind]))
        D[b,i] = 1
        TL[i] = TL[b] + 1
    end
    
    #sort by TL
    TLind = sortperm(TL)
    D = D[TLind, TLind]
    TL = TL[TLind]
    
    #add extra links
    k = rand(Poisson(λd), M)
    x = TL .- TL'
    pl = exp.(-abs.(x .- 1)/ T)
    pl[findall(D .== 1)] .= 0
    pl[diagind(pl)] .= 0
    pl .+= eps() 
    
    for i = 2:M
        D[sample(1:M, Weights(pl[i,:]), min(M,k[i]), replace = false), i] .= 1
    end
    
    return(D)
end

function coherent_community(N,M, λd, h, T)
    U = zeros(N,M)
    for i = 1:N
        U[i,rand(findall(sum(U,dims=1)[:] .== 0) )] = 1 
    end
    
    D = coherence_D(M,λd,h,T)
    [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
    D[isnan.(D)] .= 0.0

    return(Community(N,M,U,D))
end


# function overlap_community(N,M,C,o)
#     #generate uptake
#     covar = zeros(N,N)
#     covar .= o
#     covar[diagind(covar)] .= 1.0

#     means = zeros(N)
#     means .= quantile(Normal(), C)

#     mv = MultivariateNormal(means,covar)

#     U = rand(mv,M) .> 0
#     U[:,rand(1:M,rand(1:M))] .= 0 #some resources not consumed
    
#     D =  Float64.(rand(M,M) .< C) 
#     # D[:,rand(1:M,rand(1:M))] .= 0 #some resources not leaked
#     [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]
#     D[isnan.(D)] .= 0.0

#     return(Community(N,M,U,D))
# end

# function comp_community(N,M,C,a)
#     U = Float64.(rand(N,M) .< C)

#     for i = 1:N
#         s = sum(U[i,:])
#         sp_id = rand(findall(U[i,:] .!= 0)) 
#         #set all
#         U[i,:] .*= (1-a)*s / (s-1)
#         #set specialism
#         U[i, sp_id] = a*s
#     end

#     D =  Float64.(rand(M,M) .< C)
#     # D[:,1] .= 0
#     [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]

#     return(Community(N,M,U,D))
# end


# # function niche_model(N,M,D,Cu,Cl)
# #     ci = [rand(D) for i = 1:N]
# #     ci = ci[sortperm(sum.(ci))]

# #     ca = [rand(D) for i = 1:M]
# #     ca = ca[sortperm(sum.(ca))]

# #     β = (1/Cu) - 1

# #     U = MiCRM_stability.generate_U_mat(N,M,ci,ca, MiCRM_stability.d_manhat, β)
# #     L = MiCRM_stability.generate_L_mat(M,ca, MiCRM_stability.d_manhat, Cl)

# #     return(Community(N,M,U,L))
# # end

# # """
# #     rand_model(N,M, C = 0.3)

# # generate a uniform random community with a given connectance value
# # """
# # function rand_model(N,M, μ, σ, Cu, Cl)
# #     dMn = MultivariateNormal(μ,σ)
# #     x = rand(dMn, max(N,M), M)

# #     U = [cdf(Normal(μ[1],σ[1,1]),y[1]) for y = x[1:N,1:M]]
# #     L = [cdf(Normal(μ[2],σ[2,2]),y[2]) for y = x[1:M,1:M]]
    
# #     U = U .* (rand(N,M) .< Cu)
# #     L = L .* (rand(M,M) .< Cl)
# #     # U = 

# #     # for i = 1:M
# #     #     L[i,:] .= L[i,:] ./ sum(L[i,:])
# #     # end

# #     return Community(N,M,U,L)
# # end

# # function correlated_rand_model(N,M,σ, Cu, Cl)
# #     #correlation matrix
# #     V = zeros(2N, 2N)
    
# #     for i = 1:2N
# #         for j = 1:2N
# #             if i == j
# #                 V[i,j] = σ[1]
# #             elseif (i <= N) & (j <= N)
# #                 #uptake - uptake
# #                 V[i,j] = σ[2]
# #             elseif (i > N) & (j > N)
# #                 #leakage - leakage
# #                 V[i,j] = σ[3]
# #             else
# #                 #leakage - uptake
# #                 if (i == N + j) | (j == N + i)
# #                     #self uptake-leakage correlation
# #                     V[i,j] = σ[4]
# #                 else
# #                     V[i,j] = σ[5]
# #                 end
# #             end
# #         end
# #     end

# #     dMv = MultivariateNormal(V)
# #     x = rand(dMv, M)
# #     U = cdf.(Normal(), x[1:N,:]) .* (rand(N,M) .< Cu)
# #     K = cdf.(Normal(), x[N+1:end,:]) .* (rand(N,M) .< Cl)

# #     #input and output
# #     K[:,1] .= 0.0
# #     U[:,end] .= 0.0

# #     return(Community(N,M,U,K))
# # end

# # function check_community(c::Community)

# #     #every consumer has a resource
# #     cindx = sum(c.U, dims = 2)[:] .!= 0
# #     #every resource has a producer
# #     rindx = sum(c.L, dims = 1)[:] .!= 0

# #     U = c.U[cindx,rindx]
# #     K = c.L[cindx,rindx]
#     K[:,rand(1:sum(rindx))] .= 0.0

#     return(Community(sum(cindx),sum(rindx),U,K))


# end
# #distance
# du(a,b) = sqrt(sum((a .- b).^2))
# dl(a,b) = a < b ? sqrt(sum((a .- b).^2)) : Inf

# #swap
# function swap!(indx,i,j)
#     indx[i] = j
#     indx[j] = i
# end

# #discrete
# get_R(V,n) = V^(1/n) * (gamma(n/2 + 1)^(1/n)) / sqrt(n)
# f_discrete(d, R) = d < R

# #gausian
# f_gaussian(d) = exp(-((d)^2)/0.1 )


# function niche_model(N,M,C; nC = rand(N), nR = rand(M))
#     #sorting
#     nC = nC[[p for p = nC] |> sortperm]
#     nR = nR[[p for p = nR] |> sortperm]
    
#     eff_C = 1 - sqrt(1 - C)
#     RC = rand(Beta(1, (1/eff_C) - 1), N)
#     RR = rand(Beta(1, (1/eff_C) - 1), M)

#     U = [abs(nC - nR) <= RC[i] for (i,nC) = enumerate(nC), nR = nR]
#     L = [0 < nR2 - nR1 < 2RR[i]  for (i,nR1) = enumerate(nR), nR2 = nR]

#     return Community(N,M,U,L)
# end



# function block_model(N,M,A::Matrix{Int64},C)
#     # randomise blocks
#     nC = rand(1:size(A)[1], N) |> sort
#     nR = rand(1:size(A)[1], M) |> sort 
#     # generate matricies
#     U = zeros(N,M)
#     L = zeros(M,M)
#     # allocate
#     U_links = rand(N,M) .< C
#     L_links = rand(M,M) .< C

#     for i = 1:N
#         for a = 1:M
#             if nC[i] == nR[a]
#                 U[i,a] = U_links[i,a]
#             end
#         end
#     end

#     for a = 1:M
#         for b = 1:M
#             if A[nR[a],nR[b]] == 1
#                 L[a,b] = L_links[a,b]
#             end
#         end
#     end

#     return Community(N,M,U,L)

# end

# function spandrel_model(N,M,Nsteps,C,μ)
#     #allocate array
#     U = zeros(N,M)
#     L = rand(M,M) .< C
#     #add initial links
#     U[1,:] .= rand(M) .< C
    
#     #loop
#     for i = 1:Nsteps
#         #choose species
#         from = rand(1:N)
#         to = rand(1:N)

#         U[to,:] .=  U[from,:] 

#         #change
#         U[from,:] += Int.(rand(M) .< (μ/2)) .- Int.(rand(M) .< (μ/2))
#         U[to,:] += Int.(rand(M) .< (μ/2)) .- Int.(rand(M) .< (μ/2))
#         U = max.(U, 0)
#         U = min.(U, 1)
#     end

#     return Community(N,M,U,L)
# end

# function SNM_model(N,M,C, evalf::Function, t)
#     U = rand(N, M) .< C
#     L = rand(M, M) .< C
#     e = evalf(U)
    
#     k = 0
#     v = 0
#     #eval loop
#     while abs(e - t) > 0.05
#         # println(abs(e - t))
#         if k % 100 == 0
#             if v == abs(e - t)
#                 break
#             end
#             v = abs(e - t)
#             k = 0
#         end
#         k += 1

#         # row = rand(1:N)

#         ind_t = findall(U[:])
#         ind_f = findall(.!U[:])

#         if isempty(ind_t) || isempty(ind_f)
            
#         else
#             i = rand(ind_t)
#             j = rand(ind_f)

#             U[i] = 0
#             U[j] = 1 

#             if abs(evalf(U) - t) > abs(e - t)
#                 U[i] = 1
#                 U[j] = 0
#             else
#                 e = evalf(U)  
#             end
            
#         end
#     end

#     return Community(N,M,U,L)
# end