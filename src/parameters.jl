function get_structural_params(N::Int, M::Int, χ::Matrix{Float64},D::Matrix{Float64}, Λ::Vector{Float64}, I::Vector{Float64})
    
    #sample G - random for now
    G = rand(N)

    #calculate M
    # Mor = G .* (1 .- Λ)

    # #get total resource uptake and leakage
    Utot = χ' * G
    Ktot = (χ * D)' * (G .* Λ)

    #outflow
    O = I .+ Ktot .- Utot
    mask = O .< 0
    # I = copy(I)  # avoid mutating input
    I[mask] .-= O[mask]
    O = I .+ Ktot .- Utot

    U = G .* χ
    K = (G .* Λ) .* (χ * D)

    ρ = I ./ (I .+ Ktot)
    σ = O ./ (O .+ Utot)
    
    #branching terms
    γ = zeros(N,M)
    η = zeros(N,M)

    for a = 1:M
        Ucolsum = sum(U[:,a])
        γ[:,a] .= Ucolsum > 0 ? U[:,a] ./ Ucolsum : 0.0
        Kcolsum = sum(K[:,a])
        η[:,a] .= Kcolsum > 0 ? K[:,a] ./ Kcolsum : 0.0
    end

    #sample flux balances - add random influx + outflow
    ρ̃ = 1 .- ρ
    σ̃ = 1 .- σ

    # #structural
    α = rand(Uniform(0.1,1.0)) .+ rand(Uniform(-0.1,0.1),N) 
    β = rand(Uniform(0.1,1.0)) .+ rand(Uniform(-0.1,0.1),M) 

    ϕ_norm_inv = 1 ./ (χ * D)  # precomputed
    ϕ_norm_inv[isinf.(ϕ_norm_inv)] .= 0.0

    return MiCRM_stability.structural_params(α,β,Λ,ρ,ρ̃,σ,σ̃,γ,η,χ, D, ϕ_norm_inv);
end

get_structural_params(c::Community, Λ::Vector{Float64}, extF::Float64, I::Float64) = get_structural_params(c.N, c.M, c.U, c.D, Λ, I)



function get_exponential_parameters(N,M)
    σ = 0.1
    gx = rand(Uniform(0.75,1.0)) .+ zeros(N)
    gs = rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    mx = rand(Uniform(1.0,1.25)) .+ zeros(N)
    
    fy = ones(N,M)
    λy = zeros(N)

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



function get_parameters(c::Community, Λ::Vector{Float64},I::Vector{Float64}, extF::Float64)
    Parameters(c.N,c.M, 
    get_structural_params(c,Λ,extF,I), 
    get_exponential_parameters(c))
end