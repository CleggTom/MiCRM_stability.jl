struct structural_params
    #scale params
    α::Vector{Float64} #scale param for consumers
    β::Vector{Float64} #scale param for resources
    Λ::Vector{Float64} #efficency param

    #branching params
    ρ::Vector{Float64} #prop of resource gain from inflow
    ρ̃::Vector{Float64} #prop of resource gain from leakage
   
    σ::Vector{Float64} #prop of resource loss from outlfow
    σ̃::Vector{Float64} #prop of resource loss from uptake
   
    γ::Matrix{Float64} #prop of uptake loss from consumer i
    η::Matrix{Float64} #prop of leakage gain from consumer i

    χ::Matrix{Float64} #prop of substrate gain from resource a
    ϕ::Vector{Matrix{Float64}} #prop of waste product of b from uptake of a by i
end

struct exponential_params
    #consumer growth
    gx::Vector{Float64} #growth to consumer
    gs::Vector{Float64} #growth to substrates
    
    #consumer mortality
    mx::Vector{Float64} #mortality to consumer
    
    #consumer uptake + leakage stochiometry
    fy::Matrix{Float64} #substrate proportion to resource

    #Leakage structure
    λy::Matrix{Float64} #total leakage to resources

    #external resource dynamics
    iy::Vector{Float64}
    oy::Vector{Float64}
end

struct Parameters
    N::Int64 #number of consumers
    M::Int64 #number of resources
    s::structural_params
    e::exponential_params
    tmp::Vector{Matrix{Float64}}
end

#struct for storing communtiy uptake and leakage matricies
struct Community
    N::Int64
    M::Int64
    U::Matrix{Float64}
    D::Matrix{Float64}
end