function get_structural_params(U,K,N,M)

    #structural
    α = 1 .+ rand(N) #abs.(rand(Uniform(0,2)) .+ rand(Normal(0.0,0.1), N))
    β = 1 .+ rand(M) #abs.(rand(Uniform(0,2)) .+ rand(Normal(0.0,0.1), M))
    
    #branching
    #supply
    ρ = zeros(M)
    # ρ = 0.1 * rand(M)
    ρ[sum(K,dims=1)[:] .== 0] .= 1.0
    ρ̃ = 1 .- ρ

    #outflow
    σ = zeros(M)
    # σ = 0.1 * rand(M)
    σ[sum(U,dims=1)[:] .== 0] .= 1.0
    # σ[end] = 1
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
    σ = 0.1
    gx = rand(Uniform(0.75,1.0)) .+ zeros(N)
    gs = rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    gw = -rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    mx = rand(Uniform(1.0,1.25)) .+ zeros(N)
    
    # fy = rand(Uniform(1,1)) .+ rand(Uniform(-σ,σ),N,M)
    # hy = rand(Uniform(1,1)) .+ rand(Uniform(-σ,σ),N,M)

    fy = ones(N,M)
    hy = -ones(N,M)
    
    # λy = rand(Uniform(-1,1)) .+ rand(Uniform(-σ,σ),N,M)
    # ωy = rand(Uniform(-1,1)) .+ rand(Uniform(-σ,σ),N,M)
    λy = zeros(N,M)
    ωy = zeros(N,M)


    iy = rand(Uniform(0.0,0.5)) .+ rand(Uniform(-σ,σ), M)
    oy = rand(Uniform(1,2)) .+ rand(Uniform(-σ,σ), M)

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,iy,oy)
end

get_exponential_parameters(c::Community) = get_exponential_parameters(c.N, c.M)

function get_stable_exp_parameters(N,M)
    σ = 0.1
    gx = 0.99ones(N)
    gs = 2ones(N)
    gw = -ones(N)
    mx = ones(N)
    
    fy = ones(N,M)
    hy = -ones(N,M)
    
    λy = zeros(N,M)
    ωy = zeros(N,M)

    iy = zeros(M) 
    oy = ones(M) 

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,iy,oy)
end

function get_size_exp_parameters(N,M)
    σ = 0.1
    gx = ones(N) .- rand(Uniform(-1e-6,0))
    gs = 2ones(N) .+ rand(Uniform(-σ,σ), N)
    gw = -ones(N) .+ rand(Uniform(-σ,σ), N)
    
    mx = ones(N) #.+ rand(Uniform(-σ,σ), N)
    
    fy = ones(N,M)
    hy = -ones(N,M)
    
    λy = zeros(N,M)
    ωy = zeros(N,M)


    iy = zeros(M) #rand(Uniform(0,1)) .+ rand(Uniform(-σ,σ), M)
    oy = ones(M) #rand(Uniform()) .+ rand(Uniform(-σ,σ), M)

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,iy,oy)
end

function get_parameters(c::Community)
    Parameters(c.N,c.M, get_structural_params(c), get_exponential_parameters(c))
end