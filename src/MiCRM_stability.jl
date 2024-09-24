module MiCRM_stability

using LinearAlgebra, Distributions, SpecialFunctions, StatsBase, Random

# Write your package code here.
include("types.jl")
include("parameters.jl")
include("jacobian.jl")
include("niche_model.jl")
include("community_generation.jl")
include("network_metrics.jl")


end
