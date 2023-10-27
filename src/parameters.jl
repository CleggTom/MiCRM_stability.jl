function get_structural_params(U,K,N,M)

    #structural
    α = abs.(rand(Uniform(0,1)) .+ rand(Normal(0.0,0.1), N))
    β = abs.(rand(Uniform(0,1)) .+ rand(Normal(0.0,0.1), N))
    
    #branching
    ρ = zeros(M)
    ρ[sum(K,dims=1)[:] .== 0] .= 1.0
    ρ̃ = 1 .- ρ
    
    σ = zeros(M)
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

get_structural_params(c::Community) = get_structural_params(c.U, c.U * c.L, c.N, c.M)

function get_exponential_parameters(N,M)
    σ = 0.5
    gx = rand(Uniform(0,1)) .+ rand(Normal(0.0,σ), N)
    gs = rand(Uniform(0,5)) .+ rand(Normal(0.0,σ), N)
    gw = -rand(Uniform(0,5)) .+ rand(Normal(0.0,σ), N)
    
    mx = rand(Uniform(1,2)) .+ rand(Normal(0.0,σ), N)
    
    fy = ones(N,M)
    hy = ones(N,M)
    
    λy = rand(Uniform(-0.5,0.5), N,M)
    ωy = rand(Uniform(-0.5,0.5), N,M)
    θy = rand(Uniform(-0.5,0.5), N,M,M)
    
    iy = rand(Uniform(0,2)) .+ rand(Normal(0.0,σ), M)
    oy = rand(Uniform(1,2)) .+ rand(Normal(0.0,σ), M)

    return MiCRM_stability.exponential_params(gx,gs,gw,mx,fy,hy,λy,ωy,θy,iy,oy)
end

get_exponential_parameters(c::Community) = get_exponential_parameters(c.N, c.M)

function get_parameters(c::Community)
    Parameters(c.N,c.M, get_structural_params(c), get_exponential_parameters(c))
end