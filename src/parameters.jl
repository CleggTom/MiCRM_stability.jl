function get_structural_params(U,K)
    N,M = size(U)
    #structural
    α = ones(N)
    β = ones(M)
    
    #branching
    ρ = zeros(N)
    ρ[sum(K,dims=1)[:] .== 0] .= 1.0
    ρ̃ = 1 .- ρ
    
    σ = zeros(N)
    σ[sum(U,dims=1)[:] .== 0] .= 1.0
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

get_structural_params(c::Community) = get_structural_params(c.U, c.U * c.L)

function get_exponential_parameters(N,M)
    gx = (0.5) .+ (rand(N) / 2)
    gs = (0.5) .+ (rand(N))
    gw = -(0.5) .+ (rand(N))
    
    mx = ones(N)
    
    fy = ones(N,M) .+ rand(Normal(0,0.1), (N,M))
    hy = ones(N,M)
    
    λy = zeros(N,M)
    ωy = zeros(N,M)
    θy = zeros(N,M,M)
    
    iy = (0.5) .+ (rand(M) / 2)
    oy = (0.5) .+ (rand(M) / 2)

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,θy,iy,oy)
end

get_structural_params(c::Community) = get_structural_params(c.N, c.M)