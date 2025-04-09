struct structural_params
    # Scale params
    α::Vector{Float64}  # Scale param for consumers
    β::Vector{Float64}  # Scale param for resources
    Λ::Vector{Float64}  # Efficiency param
    Ω::Vector{Float64}  # Precomputed (1 - Λ)

    # Branching params
    ρ::Vector{Float64}  # Proportion of resource gain from inflow
    ρ̃::Vector{Float64}  # Proportion of resource gain from leakage

    σ::Vector{Float64}  # Proportion of resource loss from outflow
    σ̃::Vector{Float64}  # Proportion of resource loss from uptake

    γ::Matrix{Float64}  # Proportion of uptake loss from consumer i
    η::Matrix{Float64}  # Proportion of leakage gain from consumer i

    χ::Matrix{Float64}  # Proportion of substrate gain from resource a
    D::Matrix{Float64}  # Stoichometric matrix

    ϕ_norm_inv::Array{Float64,2} # Proportion of substrate production to resource a

    function structural_params(α, β, Λ, ρ, ρ̃, σ, σ̃, γ, η, χ, D, ϕ_norm_inv)
        Ω = 1.0 .- Λ  # Precompute (1 - Λ)
        new(α, β, Λ, Ω, ρ, ρ̃, σ, σ̃, γ, η, χ, D, ϕ_norm_inv)
    end
end

struct exponential_params
    # Consumer growth
    gx::Vector{Float64}  # Growth to consumer
    gs::Vector{Float64}  # Growth to substrates
    
    # Consumer mortality
    mx::Vector{Float64}  # Mortality to consumer
    
    # Consumer uptake + leakage stoichiometry
    fy::Matrix{Float64}  # Substrate proportion to resource

    # Leakage structure
    λg::Vector{Float64}  # Total leakage to uptake

    # External resource dynamics
    iy::Vector{Float64}
    oy::Vector{Float64}

    function exponential_params(gx, gs, mx, fy, λy, iy, oy)
        new(gx, gs, mx, fy, λy, iy, oy)
    end
end


struct Parameters
    N::Int64 #number of consumers
    M::Int64 #number of resources
    s::structural_params
    e::exponential_params
end

#struct for storing communtiy uptake and leakage matricies
struct Community
    N::Int64
    M::Int64
    U::Matrix{Float64}
    D::Matrix{Float64}
end