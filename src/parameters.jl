function get_structural_params(U::Matrix{Float64},D::Matrix{Float64}, Λ::Vector{Float64}, extF::Float64)
    N,M = size(U)

    #get K topology
    K = diagm(Λ) * U * D

    #uptake branching terms
    χ = U ./ (norm.(eachrow(U), 1)) #resources to consumer uptake
    γ = U ./ (norm.(eachcol(U), 1)') #consumers to resource loss
    #release branching
    η = K ./ (norm.(eachcol(K), 1)') #consumers to resource gain

    γ[isnan.(γ)] .= 0.0 
    χ[isnan.(χ)] .= 0.0 
    η[isnan.(η)] .= 0.0 

    #get inflow outflow
    #calculate excess flux
    EF = sum(K,dims=1)[:] .- sum(U,dims=1)[:]
   

    external = hcat([EF[a] < 0 ? [abs(EF[a])*(1 + extF), abs(EF[a])*extF] : [abs(EF[a])*extF, abs(EF[a])*(1 + extF)] for a = 1:M]...)
    
    ρ = external[1,:] ./ (external[1,:] .+ sum(K,dims=1)[:])
    ρ̃ = 1 .- ρ
   
    σ = external[2,:] ./ (external[2,:] .+ sum(U,dims=1)[:])
    σ̃ = 1 .- σ
   
   
    ρ[isnan.(ρ)] .= 0.0
    ρ̃[isnan.(ρ̃)] .= 0.0
    σ[isnan.(σ)] .= 0.0
    σ̃[isnan.(σ̃)] .= 0.0

    #structural
    α = sum(U,dims=2)[:]
    β = (external[1,:] .+ sum(K,dims=1)[:])

    ϕ = [zeros(M,M) for i = 1:N]
    for i = 1:N
        for a = 1:M
            for b = 1:M
                ϕ[i][a,b] = (χ[i,b] * D[b,a]) 
            end
        end
    end 
    #normalise
    [ϕ[i][a,:] .= ϕ[i][a,:]  / sum(ϕ[i][a,:] ) for i = 1:N, a = 1:M]
    [x[isnan.(x)] .= 0.0 for x = ϕ]

    return MiCRM_stability.structural_params(α,β,Λ,ρ,ρ̃,σ,σ̃,γ,η,χ,ϕ);
end

get_structural_params(c::Community, Λ::Vector{Float64}, extF::Float64) = get_structural_params(c.U, c.D, Λ, extF)



function get_exponential_parameters(N,M)
    σ = 0.1
    gx = rand(Uniform(0.75,1.0)) .+ zeros(N)
    gs = rand(Uniform(0,3)) .+ rand(Uniform(-σ,σ), N)
    mx = rand(Uniform(1.0,1.25)) .+ zeros(N)
    
    fy = ones(N,M)

    λy = zeros(N,M)

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



function get_parameters(c::Community, Λ::Vector{Float64},extF::Float64)
    Parameters(c.N,c.M, 
    get_structural_params(c,Λ,extF), 
    get_exponential_parameters(c),
    [])
end