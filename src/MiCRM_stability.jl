module MiCRM_stability

using LinearAlgebra, Distributions, SpecialFunctions, StatsBase

# Write your package code here.
include("types.jl")
include("parameters.jl")
include("jacobian.jl")
include("niche_model.jl")
include("community_generation.jl")


end
