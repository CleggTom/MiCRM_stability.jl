function get_structural_params(U,D,Λ)
    N,M = size(U)

    #get K matrix
    G = sum(U,dims=2)[:]
    #consumer branching terms
    χ = U ./ (norm.(eachrow(U), 1))
    γ = U ./ (norm.(eachcol(U), 1)')
    #by-product release term
    K = G .* Λ .* (χ * D)
    #resource branching
    η = K ./ (norm.(eachcol(K), 1)')
    ϕ = K ./ (norm.(eachrow(K), 1))
 
    γ[isnan.(γ)] .= 0.0 
    η[isnan.(η)] .= 0.0 
    χ[isnan.(χ)] .= 0.0 
    ϕ[isnan.(ϕ)] .= 0.0 

    #branching
    #get excess flux
    flux = sum(K,dims=1)[:] .- sum(U,dims=1)[:]
    
    #add outflow to all
    Inflow = -minimum(flux)
    Outflow = flux .+ Inflow

    #supply
    # ρ = zeros(M)
    ρ = Inflow ./ (Inflow .+ sum(K,dims=1)[:])
    ρ̃ = 1 .- ρ

    #outflow
    σ = Outflow ./ (Outflow .+ sum(U,dims=1)[:])
    σ̃ = 1 .- σ
    
    #structural
    α = rand(N) #abs.(rand(Uniform(0,2)) .+ rand(Normal(0.0,0.1), N))
    β = rand(M) #abs.(rand(Uniform(0,2)) .+ rand(Normal(0.0,0.1), M))

    ψ = zeros(N,M,M)
    for i = 1:N
        for a = 1:M
            s = 0 
            for c = 1:M
                s += χ[i,a] * D[c,a]
            end

            for b = 1:M
                if s != 0
                    ψ[i,b,a] = (χ[i,a] * D[b,a]) / s
                else
                    ψ[i,b,a] = 0.0 
                end
            end
        end
    end 

    return MiCRM_stability.structural_params(α,β,Λ,ρ,ρ̃,σ,σ̃,γ,η,χ,ϕ,ψ);
end

get_structural_params(c::Community) = get_structural_params(c.U, c.L, c.N, c.M)



function get_exponential_parameters(N,M)
    σ = 0.1
    gx = rand(Uniform(0.75,1.0)) .+ zeros(N)
    gs = rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    gw = -rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    mx = rand(Uniform(1.0,1.25)) .+ zeros(N)
    
    # fy = rand(Uniform(1,1)) .+ rand(Uniform(-σ,σ),N,M)
    # hy = rand(Uniform(1,1)) .+ rand(Uniform(-σ,σ),N,M)

    fy = ones(N,M)
    hy = ones(N,M)
    
    # λy = rand(Uniform(-1,1)) .+ rand(Uniform(-σ,σ),N,M)
    # ωy = rand(Uniform(-1,1)) .+ rand(Uniform(-σ,σ),N,M)
    λy = zeros(N,M)

    iy = rand(Uniform(0.0,0.5)) .+ rand(Uniform(-σ,σ), M)
    oy = rand(Uniform(1,2)) .+ rand(Uniform(-σ,σ), M)

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,iy,oy)
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



function get_parameters(c::Community)
    Parameters(c.N,c.M, 
    get_structural_params(c), 
    get_exponential_parameters(c), 
    util_params(c.N,c.M,
        get_structural_params(c),
        get_exponential_parameters(c)))
end