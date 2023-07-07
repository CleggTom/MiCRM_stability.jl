struct Parameters
    N::Int64 #number of consumers
    M::Int64 #number of resources

    #scale params
    a_x::Vector{Float64} #scale param for consumers
    a_y::Vector{Float64} #scale param for resources

    #Proportions
    β::Matrix{Float64} #prop of growth of i from resource a
    T::Matrix{Float64} #total leakage prop of a by i
    ψ::Array{Float64, 3} #prop of leakage of a->b by i

    ρ::Vector{Float64} #prop of resource gain from inflow
    ρ̃::Vector{Float64}  #prop of resource gain from leakage
    ϕ::Array{Float64, 3}   #prop of leakage of b->a by i

    σ::Vector{Float64}  #prop of resource loss from outlfow
    σ̃::Vector{Float64}  #prop of resource loss from uptake
    γ::Matrix{Float64}  #prop of resource loss from uptake by i

    #exponential
    μ_x::Matrix{Float64} #elasticity of uptake of a by i to change in consumer mass
    μ_y::Matrix{Float64} #elasticity of uptake of a by i to change in resource mass
    ζ_x::Array{Float64, 3} #elasticity of leakage to consumer
    ζ_S::Array{Float64, 3}
    ζ_P::Array{Float64, 3} #elasticity of inflow to resources
    ω::Vector{Float64} #elasticity of consumer mortality
    ι::Vector{Float64} #elasticity of resource inflow
    ϱ::Vector{Float64} #elasticity of resource outlfow
end