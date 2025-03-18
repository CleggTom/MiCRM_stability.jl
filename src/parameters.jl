function get_structural_params(U::Matrix{Float64},D::Matrix{Float64}, Λ::Vector{Float64})
    N,M = size(U)

    #get K topology
    K = diagm(Λ) * U * D
    
    #branching terms
    χ = zeros(N,M)
    γ = zeros(N,M)
    η = zeros(N,M)

    for i = 1:N
        χ[i,:] .= sum(U[i,:]) > 0 ? U[i,:] ./ sum(U[i,:]) : 0.0 #resources to consumer uptake
    end
    for a = 1:M
        γ[:,a] .= sum(U[:,a]) > 0 ? U[:,a] ./ sum(U[:,a]) : 0.0 #consumers to resource loss
        η[:,a] .= sum(K[:,a]) > 0 ? K[:,a] ./ sum(K[:,a]) : 0.0 #consumers to resource gain
    end

    #resource inflow/outflow - get the total internal flux for each resource
    Utot = sum(U, dims = 1)[:]
    Ktot = sum(K, dims = 1)[:]

    #get minimal external flux possible
    EF = Utot .- Ktot
    ρ = zeros(M)
    σ = zeros(M)
    for a = 1:M
        if EF[a] > 0
            ρ[a] = abs(EF[a]) / (abs(EF[a]) + Ktot[a])
        elseif EF[a] < 0
            σ[a] = abs(EF[a]) / (abs(EF[a]) + Utot[a])
        else
        end
    end

    #sample flux balances - add random influx + outflow
    ρ .+= (1 .- ρ) .* rand(M)
    σ .+= (1 .- σ) .* rand(M)
    ρ̃ = 1 .- ρ
    σ̃ = 1 .- σ

    # #structural
    α = rand() .+ rand(N) ./ 10
    β = rand() .+ rand(M) ./ 10

    ϕ = Matrix{Vector{Float64}}(undef,N,M)
    pmat = χ * D
    #match indexing in SI
    for i = 1:N
            for a = 1:M
                ϕ[i,a] = zeros(M)
                if pmat[i,a] > 0
                    for b = 1:M
                        ϕ[i,a][b] = (χ[i,b] * D[b,a])/ pmat[i,a]
                    end     
                end
            end
    end 

    # #normalise
    # [ϕ[i,a] .= sum(ϕ[i,a]) > 0 ? ϕ[i,a]  / sum(ϕ[i,a]) : ϕ[i,a] for i = 1:N, a = 1:M]

    return MiCRM_stability.structural_params(α,β,Λ,ρ,ρ̃,σ,σ̃,γ,η,χ,ϕ);
end

get_structural_params(c::Community, Λ::Vector{Float64}, extF::Float64) = get_structural_params(c.U, c.D, Λ)



function get_exponential_parameters(N,M)
    σ = 0.1
    gx = rand(Uniform(0.75,1.0)) .+ zeros(N)
    gs = rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    mx = rand(Uniform(1.0,1.25)) .+ zeros(N)
    
    fy = ones(N,M)

    λy = zeros(N,M)

    iy = rand(Uniform(0.0,0.5)) .+ rand(Uniform(-σ,σ), M)
    oy = rand(Uniform(1,2)) .+ rand(Uniform(-σ,σ), M)

    return MiCRM_stability.exponential_params(gx,gs,mx,fy,λy,iy,oy)
end

get_exponential_parameters(c::Community) = get_exponential_parameters(c.N, c.M)

# function get_stable_exp_parameters(N,M)
#     σ = 0.1
#     gx = 0.99ones(N)
#     gs = 2ones(N)
#     gw = -ones(N)
#     mx = ones(N)
    
#     fy = ones(N,M)
#     hy = -ones(N,M)
    
#     λy = zeros(N,M)
#     ωy = zeros(N,M)

#     iy = zeros(M) 
#     oy = ones(M) 

#     return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,iy,oy)
# end

# function get_size_exp_parameters(N,M)
#     σ = 0.1
#     gx = ones(N) .- rand(Uniform(-1e-6,0))
#     gs = 2ones(N) .+ rand(Uniform(-σ,σ), N)
#     gw = -ones(N) .+ rand(Uniform(-σ,σ), N)
    
#     mx = ones(N) #.+ rand(Uniform(-σ,σ), N)
    
#     fy = ones(N,M)
#     hy = -ones(N,M)
    
#     λy = zeros(N,M)
#     ωy = zeros(N,M)


#     iy = zeros(M) #rand(Uniform(0,1)) .+ rand(Uniform(-σ,σ), M)
#     oy = ones(M) #rand(Uniform()) .+ rand(Uniform(-σ,σ), M)

#     return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,iy,oy)
# end



function get_parameters(c::Community, Λ::Vector{Float64},extF::Float64)
    Parameters(c.N,c.M, 
    get_structural_params(c,Λ,extF), 
    get_exponential_parameters(c))
end