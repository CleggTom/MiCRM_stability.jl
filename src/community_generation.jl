#distance
du(a,b) = sqrt(sum((a .- b).^2))
dl(a,b) = a < b ? sqrt(sum((a .- b).^2)) : Inf

#gausian 
f_gaussian(d) = exp(-((d)^2)/0.1 )

function niche_model(N,M,D; 
        nC =[rand(D) for i = 1:N],
        nR = [rand(D) for i = 1:M],
        f = f_gaussian)

    #sorting
    nC = nC[[p[1] for p = nC] |> sortperm]
    nR = nR[[p[1] for p = nR] |> sortperm]

    U = [f(du(nC,nR)) for nC = nC, nR = nR]
    L = [f(dl(nR1,nR2)) for nR1 = nR, nR2 = nR]

    return Community(N,M,U,L)
end