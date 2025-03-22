# function calculate_g_s(p::Parameters)
#     res = zeros(p.N,p.M)
#     for i = 1:p.N
#         for a = 1:p.M
#             res[i,a] = p.e.fy[i,a] * p.s.χ[i,a] * (p.e.gs[i] - 1)
#         end
#     end
#     return(res)
# end

# Optimized Jacobian function with reduced setindex! operations
function jacobian!(p::Parameters, J::Matrix{Float64})
    J .= 0
    # Loop over consumers
    @inbounds for i = 1:p.N
        α_i = p.s.α[i]
        Ω_i = p.s.Ω[i]
        gx_i = p.e.gx[i]
        mx_i = p.e.mx[i]
        gs_i = p.e.gs[i]
        Λ_i = p.s.Λ[i]
        
        #consumer intraspecific element
        @fastmath J[i, i] = α_i * Ω_i * (gx_i - mx_i)

        #consumer to resources
        for a = 1:p.M
            χ_ia = p.s.χ[i, a]
            fy_ia = p.e.fy[i, a]
            λy_ia = p.e.λy[i, a]
            
            # dxi/dya
            @fastmath J[i, a + p.N] = α_i * Ω_i * ( (gs_i * χ_ia * fy_ia) - Λ_i * λy_ia)
        end
    end

    @inbounds for a = 1:p.M
        β_a = p.s.β[a]
        ρ̃_a = p.s.ρ̃[a]
        σ̃_a = p.s.σ̃[a]
        ρ_a = p.s.ρ[a]
        σ_a = p.s.σ[a]
        iy_a = p.e.iy[a]
        oy_a = p.e.oy[a]
        
        #resource to consumer
        for i = 1:p.N
            η_ia = p.s.η[i, a]
            γ_ia = p.s.γ[i, a]
            gx_i = p.e.gx[i]
            
            # dya/dxi
            @fastmath J[a + p.N, i] += β_a * gx_i * (ρ̃_a * η_ia - σ̃_a * γ_ia)
        end

        # Uptake and leakage
        @inbounds @simd for b = 1:p.M
            idx = a + p.N  # Precompute row index
            
            # External supply and outflow
            if a == b
                @fastmath J[idx, idx] += ρ_a * iy_a - σ_a * oy_a
            end
            
            # Internal dynamics
            for i = 1:p.N
                γ_ia = p.s.γ[i, a]
                η_ia = p.s.η[i, a]
                gs_i = p.e.gs[i]
                
                fy_ib = p.e.fy[i, b]
                χ_ib = p.s.χ[i, b]
                λy_ib = p.e.λy[i, b]
                
                uptake_term = σ̃_a * γ_ia * (-fy_ib * χ_ib + gs_i * χ_ib * fy_ib)
                # Use a temporary buffer to accumulate results
                temp_J = J[idx, b + p.N] - uptake_term
                
                # Direct loss term - diagonal only
                if a == b
                    temp_J -= σ̃_a * γ_ia * fy_ib
                end
                
                # Leakage
                leakage_direct = ρ̃_a * η_ia * p.s.ϕ[i, b, a] * fy_ib
                temp_J += leakage_direct
                
                # Indirect terms
                @inbounds for c = 1:p.M
                    ϕ_ic_a = p.s.ϕ[i, c, a]
                    indirect_leakage = ρ̃_a * η_ia * ϕ_ic_a * (-fy_ib * χ_ib + gs_i * fy_ib * χ_ib + λy_ib)
                    temp_J += indirect_leakage
                end
                
                # Single assignment to J to minimize setindex! overhead
                J[idx, b + p.N] = temp_J
            end
        end
    end
end







##Metrics
function gc!(G, J, x0)
    # G = similar(J)
    for i = eachindex(x0)
        x0[i] = 1.0
        G[i,:] = -J \ x0
        x0[i] = 0.0
    end
end


"""
    se(J, i)

calculate sensitvtiy of the ith element from a jacobian
"""
function se(J) 
    L = size(J)[1]
    vals = eigvals(J)
    v = eigvecs(J)
    s = zeros(Complex, L)
    for i = 1:L
        for n = eachindex(vals)
            s[i] += norm(v[i,n]) / vals[n]
        end
    end
    return(log.(-s)) 
end

"""
    im(J, i)

calculate the impact of perturbing the ith element in a system
"""
function im(J) 
    L = size(J)[1]
    Jinv = inv(J)
    vals = eigvals(Jinv)
    w = eigvecs(Jinv)
    s = zeros(Complex, L)
    for i = 1:L
        for n = eachindex(vals)
            s[i] += norm(w[i,n]) / vals[n]
        end
    end
    return(log.(-s)) 
end