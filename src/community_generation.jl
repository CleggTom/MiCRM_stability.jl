#distance
du(a,b) = sqrt(sum((a .- b).^2))
dl(a,b) = a < b ? sqrt(sum((a .- b).^2)) : Inf

#discrete
get_R(V,n) = V^(1/n) * (gamma(n/2 + 1)^(1/n)) / sqrt(n)

f_discrete(d, R) = d < R


#gausian
f_gaussian(d) = exp(-((d)^2)/0.1 )

function niche_model(N,M,C, nC = rand(N), nR = rand(M))

    #sorting
    nC = nC[[p for p = nC] |> sortperm]
    nR = nR[[p for p = nR] |> sortperm]
    
    RC = rand(Beta(1, (1/C) - 1), N)
    RR = rand(Beta(1, (1/C) - 1), N)

    U = [abs(nC - nR) < RC[i] for (i,nC) = enumerate(nC), nR = nR]
    L = [0 < nR2 - nR1 < 2RR[i]  for (i,nR1) = enumerate(nR), nR2 = nR]

    return Community(N,M,U,L)
end

function rand_model(N,M;D=1, C = 0.3)
    U = rand(N,M) .< C
    L = rand(M,M) .< 2C

    return Community(N,M,U,L)
end