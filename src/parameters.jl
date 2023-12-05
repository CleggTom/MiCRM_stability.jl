function get_structural_params(U,K,N,M)

    #structural
    α = 1 .+ rand(N) #abs.(rand(Uniform(0,2)) .+ rand(Normal(0.0,0.1), N))
    β = 1 .+ rand(M) #abs.(rand(Uniform(0,2)) .+ rand(Normal(0.0,0.1), M))
    
    #branching
    ρ = rand(M)
    # ρ[sum(K,dims=1)[:] .== 0] .= 1.0
    ρ̃ = 1 .- ρ
    
    σ = rand(M)
    # σ[sum(U,dims=1)[:] .== 0] .= 1.0
    σ̃ = 1 .- σ
    
    γ = U ./ (norm.(eachcol(U), 1)')
    η = K ./ (norm.(eachcol(K), 1)')
    χ = U ./ (norm.(eachrow(U), 1))
    ϕ = K ./ (norm.(eachrow(K), 1))
    
    γ[isnan.(γ)] .= 0.0 
    η[isnan.(η)] .= 0.0 
    χ[isnan.(χ)] .= 0.0 
    ϕ[isnan.(ϕ)] .= 0.0 

    return MiCRM_stability.structural_params(α,β,ρ,ρ̃,σ,σ̃,γ,η,χ,ϕ);
end

get_structural_params(c::Community) = get_structural_params(c.U, c.U * c.L, c.N, c.M)

function get_exponential_parameters(N,M)
    σ = 0.25
    gx = rand(Uniform(0,2)) .+ rand(Uniform(-σ,σ), N)
    gs = rand(Uniform(0,2)) .+ rand(Uniform(-σ,σ), N)
    gw = -rand(Uniform(0,2)) .+ rand(Uniform(-σ,σ), N)
    
    mx = rand(Uniform(1,2)) .+ rand(Uniform(-σ,σ),N)
    
    fy = ones(N,M) .+ rand(Uniform(-σ,σ), N, M)
    hy = ones(N,M) .+ rand(Uniform(-σ,σ), N, M)
    
    λy = zeros(N,M) .+ rand(Uniform(-σ,σ), N, M)
    ωy = zeros(N,M) .+ rand(Uniform(-σ,σ), N, M)
    
    iy = rand(Uniform(0,0.5)) .+ rand(Uniform(-σ,σ), M)
    oy = rand(Uniform(1.5,2.5)) .+ rand(Uniform(-σ,σ), M)

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,iy,oy)
end

get_exponential_parameters(c::Community) = get_exponential_parameters(c.N, c.M)

function get_parameters(c::Community)
    Parameters(c.N,c.M, get_structural_params(c), get_exponential_parameters(c))
end