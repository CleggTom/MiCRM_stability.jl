# d_euclid(x,y) = sqrt( sum( (x .- y).^2 ) )
# d_manhat(x,y) = sum( abs.(x .- y))

# function rand_community(N,M,C,Nin, Nout)
#     U = rand(N,M) .* (rand(N,M) .< C)
#     L = rand(N,M) .* (rand(N,M) .< C)

#     U[:,rand(1:M, Nout)] .= 0
#     L[:,rand(1:M, Nin)] .= 0

#     return( Community(N,M,U,L))
# end

# function specialist_community(N,M,p,o,k,Nin,Nout)
#     U = fill(0.01, N,M) .+ (o .* (rand(N,M) .< p))
#     L = fill(0.01, N,M)
    
#     #sample L based on Usums
#     Usum = sum(U,dims=1)[:]
#     probs = (Usum .- minimum(Usum)) ./ (maximum(Usum) - minimum(Usum))
#     probs = [x > 0 ? x .^ k : 0 for x = probs]

#     for i = 1:N
#         L[i,sample(1:M, Weights(probs), rand(Binomial(M, p)), replace = false)] .+= o
#     end

#     U[:,rand(1:M, Nout)] .= 0
#     L[:,rand(1:M, Nin)] .= 0

#     return( Community(N,M,U,L))
# end


# """
#     rand_com(N,M,C, ci = rand(N), ca = rand(M))

# Niche model
# """
# function niche_community(N,M,C, ci = rand(N), ca = rand(M))
#     ci = ci[sortperm(sum.(ci))]    
#     ca = ca[sortperm(sum.(ca))]
#     d = ci .- ca'
    
#     fu(x, dist) = dist > x >= 0 ? rand() : 0.0
#     fl(x, dist) = -dist < x < 0 ? rand() : 0.0
    
#     U = zeros(N,M)
#     L = zeros(N,M)

#     b = (1/C - 1)
#     dβ = Beta(1, b)

#     for i = 1:N
#         # ind = findmin(abs.(d[i,:]))[2]
#         x = d[i,:] #.- d[i,ind]
#         U[i,:] .= fu.(x, rand(dβ))
#         L[i,:] .= fl.(x, rand(dβ))
#     end

#     #end products    
#     U[:,end] .= 0.0
#     L[:,1] .= 0.0

#     return(MiCRM_stability.Community(N,M,U,L))
# end


# function modular_community(N,M,Nclass,pgen,Cg,Co,Ci,Csp,Csw,Css)
#     res_class = repeat(collect(1:Nclass), inner = M ÷ Nclass)
#     rclass = Dict(i => findall(res_class .== i) for i = 1:Nclass)

#     #Consumers
#     U = zeros(N,M)
#     con_class = zeros(N)
#     for i = 1:N
#         if rand() < pgen
#             con_class[i] = 0
#             U[i,:] .= rand(M) .< Cg
#         else
#             con_class[i] = rand(1:(Nclass-1))
#             U[i,:] .= rand(M) .< Co
#             U[i, rclass[con_class[i]],:] .= rand(length(rclass[con_class[i]])) .< Ci
#         end
#     end
#     U[:, rclass[Nclass]] .= 0
#     U = U[sortperm(con_class),:]

#     D = zeros(M,M)
#     for a = 1:M
#         #all leakage
#         D[a, :] .= rand(M) .< Csp
#         #leakage to waste
#         D[a, rclass[Nclass]] .= rand(length(rclass[Nclass])) .< Csw
#         #leakage to same group
#         D[a, rclass[res_class[a]]] .= rand(length(rclass[res_class[a]])) .< Css
#     end

#     D[:, rclass[1]] .= 0

#     # [D[i,:] .= D[i,:] ./ sum(D[i,:]) for i = 1:M]

#     return(Community(N,M, U, D))

# end

# function nested_community(N,M,κ,v, Nin, Nout)
#     U = zeros(N,M)
#     L = zeros(N,M)

#     p = (1.0:M) .^ κ
#     p = p ./ sum(p)
    
#     kU = sample(1:M, Weights(p), N) |> sort |> reverse
#     kL = sample(1:M, Weights(p), N) |> sort |> reverse
    
#     for i = 1:N
#         Z = kU[i]
#         nj = (sum(U,dims=1)[:] .+ 1) .^ v
#         p = nj ./ sum(nj)
#         U[i,:] .= rand(Multinomial(Z,p)) .> 0

#         Z = kL[i]
#         nj = (sum(L,dims=1)[:] .+ 1) .^ v
#         p = nj ./ sum(nj)
#         L[i,:] .= rand(Multinomial(Z,p)) .> 0

        
#     end
#     ind = sortperm(sum(U .> 0, dims = 1)[:], rev = true)

#     U[:,rand(1:M, Nout)] .= 0
#     L[:,rand(1:M, Nin)] .= 0
    

#     return(Community(N,M, U[:,ind] .* rand(N,M),  L[:,ind] .* rand(N,M)))
# end

# function empirical_communty(N, C, D)
#     ind = rand(1:size(C)[1], N)
#     return(Community(N,size(C)[2],C[ind,:],D[ind,:]))
# end

# function attachment_community(N, M, C, T, pa = (x,T) -> exp(-abs(x - 1) / T))
#     #setup matricies
#     U = zeros(N,M)
#     L = zeros(N,M)

#     #initial skeleton
#     #first resource not leaked
#     B = 10
#     U[1,1:B] .= 1
    
#     ind = shuffle(1:(N+M-B-1))
#     x = vcat(fill(:N,N-1),fill(:M, M-1))[ind]

#     #numbers of links
#     l_vec = fill(1, N+M-B-1)
    
#     sN,sM = ones(N), ones(M)
#     kN,kM = 1,B
#     for i = eachindex(x)
#         if x[i] == :N
#             kN += 1
#             #inital connection
#             v = rand(1:kM)
#             U[kN,v] = 1
#             sN[kN] = sM[v] + 1
    
#         #     #other links
#         #     p_vec = pa.(sN[kN] .- sM[1:kM], T)
#         #     l = min(Int.(floor.(kM * l_vec[i])), sum(p_vec .> 0)) #
#         #     U[kN, sample(1:kM, Weights(p_vec), l, replace = false)] .= 1
            
#         else
#             kM += 1
#             #inital connection
#             v = rand(1:kN)
#             L[v, kM] = 1
#             sM[kM] = sN[v] + 1
    
#         #     #other links
#         #     p_vec = pa.(sM[kM] .- sN[1:kN], T) 
#         #     l = min(Int.(floor.(kN * l_vec[i])), sum(p_vec .> 0)) #
            
#         #     L[sample(1:kN, Weights(p_vec), l, replace = false), kM] .= 1
#         end
#     end

#     #get distancesa
#     Udiff = sN .- sM'
#     Ldiff = sM .- sN'
#     #convert to probs
#     pU = (pa.(Udiff, T))
#     pL = (pa.(Ldiff', T)) 
#     #normalise to desired connectance
#     Nl = C * ((N * M) - sum(U + L))
#     pU .= Nl .* pU ./ sum(pU)
#     pL .= Nl .* pL ./ sum(pL)

#     U .= min.(1, U .+ (pU .> rand(N,M)))
#     L .= min.(1, L .+ (pL .> rand(N,M)))
    
#     # #add additional links
#     # Up_mat = pa.(sN .- sM', T) .* rand(N,M) .* C
#     # Lp_mat = pa.((sM .- sN')', T) .* rand(N,M) .* C

#     # #add to matricies
#     # U .+= Up_mat
#     # L .+= Lp_mat

#     # #add for all but first resource
#     # kN,kM = 1,B
#     # for i = eachindex(x)
#     #     if x[i] == :N
#     #         kN += 1
#     #         l = min(Int.(floor.(kM * l_vec[i])), sum(p_vec .> 0))
#     #         v = sample(1:M, Weights(p_mat[kN,:]), l_vec[i])
#     #         U[kN, v] .= 1
#     #     else
#     #         kM += 1
#     #         v = sample(1:N, Weights(p_mat[:,kM]), l_vec[i])
#     #         U[kN, v] .= 1
#     #     end
#     # end


#     return(MiCRM_stability.Community(N, M, U, L))
# end
