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
    ϕ::Matrix{Float64} #prop of waste product to resource a
    ψ::Array{Float64} #prop of product of a released as b
end

struct exponential_params
    #consumer growth
    gx::Vector{Float64} #growth to consumer
    gs::Vector{Float64} #growth to substrates
    gw::Vector{Float64} #growth to waste products
    
    #consumer mortality
    mx::Vector{Float64} #mortality to consumer
    
    #consumer uptake + leakage stochiometry
    fy::Matrix{Float64} #substrate proportion to resource
    hy::Matrix{Float64} #waste proportion to resource

    #Leakage structure
    λy::Matrix{Float64} #total leakage to resources

    #external resource dynamics
    iy::Vector{Float64}
    oy::Vector{Float64}
end

#stores precalculated values
struct util_params
    py::Array{Float64} #y - derivatives of p_ia
    wy::Matrix{Float64} #derivatives of waste product
end

#fucntion to precalculate
function util_params(N::Int64, M::Int64, s::structural_params, e::exponential_params)
    py = similar(s.ψ)
    wy = similar(s.ϕ)

    for i = 1:N
        for a = 1:M
            for b = 1:M
                #d p_ia / d y_b
                py[i,a,b] = e.fy[i,b] * (s.ψ[i,b,a] - s.χ[i,b])
            end
        end
    end

    for i = 1:N
        for a = 1:M
            #dw_i / d_ya
            wy[i,a] = s.ϕ[i,a] * e.hy[i,a]
            for b = 1:M
             #dw_i / d_ya
             wy[i,a] +=  s.ϕ[i,b] * py[i,b,a]
            end
        end
    end

    return(util_params(py,wy))
end


struct Parameters
    N::Int64 #number of consumers
    M::Int64 #number of resources
    s::structural_params
    e::exponential_params
    u::util_params
end

#struct for storing communtiy uptake and leakage matricies
struct Community
    N::Int64
    M::Int64
    U::Matrix{Float64}
    D::Matrix{Float64}
end