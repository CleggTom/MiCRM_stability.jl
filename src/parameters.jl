"""
    get_structural_params(N::Int, M::Int, χ::Matrix{Float64},D::Matrix{Float64},
     Λ::Vector{Float64}, G::Vector{Float64}, α::Vector{Float64}, β::Vector{Float64},
     C_pert::Vector{Float64})

Generate structural parameters from a nomalised uptake (χ) and stoichometric (D) matrix. Also requires input parameters for consumer leaakage (Λ), growth (G) and the scale parameters (α, β).
"""
function get_structural_params(N::Int, M::Int, χ::Matrix{Float64},D::Matrix{Float64},
     Λ::Vector{Float64}, G::Vector{Float64}, α::Vector{Float64}, β::Vector{Float64},
     C_pert::Vector{Float64})
    
    #get total resource uptake and leakage
    Utot = χ' * G
    Ktot = (χ * D)' * (G .* Λ)

    #net demand in community (to be satisfied by ext env)
    Cnet = Utot - Ktot 
 
    #allocate inflow and outflow
    I = zeros(M)
    I .= C_pert
    O = zeros(M)
    O .= C_pert

    #find minimal environmental fluxes
    for a = 1:M
        if Cnet[a] > 0
            I[a] += Cnet[a]
        else
            O[a] -= Cnet[a] 
        end
    end

    ρ = I ./ (I .+ Ktot)
    σ = O ./ (O .+ Utot)
    
    
    #branching terms
    γ = zeros(N,M)
    η = zeros(N,M)

    U = G .* χ
    K = (G .* Λ) .* (χ * D)

    for a = 1:M
        Ucolsum = sum(U[:,a])
        γ[:,a] .= Ucolsum > 0 ? U[:,a] ./ Ucolsum : 0.0
        Kcolsum = sum(K[:,a])
        η[:,a] .= Kcolsum > 0 ? K[:,a] ./ Kcolsum : 0.0
    end

    #sample flux balances - add random influx + outflow
    ρ̃ = 1 .- ρ
    σ̃ = 1 .- σ

    ϕ_norm_inv = 1 ./ (χ * D)  # precomputed
    ϕ_norm_inv[isinf.(ϕ_norm_inv)] .= 0.0

    return MiCRM_stability.structural_params(α, β, Λ, ρ, ρ̃, σ, σ̃, γ, η, χ, D, ϕ_norm_inv);
end
