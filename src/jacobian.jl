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
function jacobian!(p::Parameters, J::AbstractMatrix{Float64})
    J .= 0
    # Loop over consumers
    @inbounds for i = 1:p.N
        α_i = p.s.α[i]
        Ω_i = p.s.Ω[i]
        gx_i = p.e.gx[i]
        mx_i = p.e.mx[i]
        gs_i = p.e.gs[i]
        λg_i = p.e.λg[i]
        Λ_i = p.s.Λ[i]
        
        #consumer intraspecific element
        J[i, i] = α_i * (gx_i * (1 - Λ_i * (1 + λg_i)) - Ω_i * mx_i)

        #consumer to resources
        for a = 1:p.M
            χ_ia = p.s.χ[i, a]
            fy_ia = p.e.fy[i, a]            
            # dxi/dya
            J[i, a + p.N] = α_i * gs_i * χ_ia * fy_ia *(1 - Λ_i * (1 + λg_i))
        end
    end

    #### Resource dynamics (C + D blocks) ####
    @inbounds for a = 1:p.M
        idx_a = a + p.N
        β_a = p.s.β[a]
        ρ̃_a = p.s.ρ̃[a]
        σ̃_a = p.s.σ̃[a]
        ρ_a = p.s.ρ[a]
        σ_a = p.s.σ[a]
        iy_a = p.e.iy[a]
        oy_a = p.e.oy[a]


        #C-block
        @inbounds for i = 1:p.N
            gx_i = p.e.gx[i]
            η_ia = p.s.η[i, a]
            γ_ia = p.s.γ[i, a]
            λg_i = p.e.λg[i]

            # C block (resource response to consumers)
            J[idx_a, i] += β_a * gx_i * (ρ̃_a * η_ia * (1 + λg_i) - σ̃_a * γ_ia)
        end

        #D-block
        # Inflow/outflow ONCE
        J[idx_a, idx_a] += ρ_a * iy_a - σ_a * oy_a
    
        @inbounds for b = 1:p.M
            idx_b = b + p.N
            D_ba = p.s.D[b, a]
            Tmp = 0.0
            @inbounds  @simd for i = 1:p.N
                gs_i = p.e.gs[i]
                gs_shift = (gs_i - 1)
                λg_i = p.e.λg[i]

                fy_ib = p.e.fy[i, b]
                χ_ib = p.s.χ[i, b]
                

                η_ia = p.s.η[i, a]
                γ_ia = p.s.γ[i, a]
                ϕ_norm_inv_ia = p.s.ϕ_norm_inv[i, a]


                ϕ_iba = (χ_ib * D_ba) * ϕ_norm_inv_ia
 
                du = gs_shift * χ_ib * fy_ib 
                dk = ϕ_iba * fy_ib + χ_ib * fy_ib * (gs_i * (1 + λg_i) - 1)
                
                uptake_term = σ̃_a * γ_ia * (du + (a==b)*fy_ib)
                leakage_term = ρ̃_a * η_ia * dk

                Tmp += leakage_term - uptake_term
            end
             J[idx_a, idx_b] += Tmp
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