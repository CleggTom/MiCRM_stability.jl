function calculate_g_s(p::Parameters)
    res = zeros(p.N,p.M)
    for i = 1:p.N
        for a = 1:p.M
            res[i,a] = p.e.fy[i,a] * p.s.χ[i,a] * (p.e.gs[i] - 1)
        end
    end
    return(res)
end

#jacobian
function jacobian!(p::Parameters, J)
    J .= 0
    #loop over consumers
    for i = 1:p.N
        for j = 1:p.N
            #dxi/dxi
            if i == j
                J[i,i] = p.s.α[i] * (1 - p.s.Λ[i]) * (p.e.gx[i] - p.e.mx[i])
            end
        end
        
        for a = 1:p.M
            #dxi/dya
            J[i, a + p.N] = p.s.α[i] * ( 
                (1 - p.s.Λ[i]) * ((p.e.gs[i] * p.s.χ[i,a] * p.e.fy[i,a])) +
                (p.s.Λ[i]) * (p.e.λy[i,a])
            )
        end

    end

    for a = 1:p.M
        for i = 1:p.N
            #dya/dxi
            J[a + p.N, i] = p.s.β[a] * p.e.gx[i] * (
                p.s.ρ̃[a] * p.s.η[i,a] - p.s.σ̃[a] * p.s.γ[i,a]
            )
        end

        for b = 1:p.M
            #dya/dya
            if a == b
                #external supply and outflow
                J[a + p.N, a + p.N] = p.s.ρ[a] * p.e.iy[a] - p.s.σ[a] * p.e.oy[a]
            end

            for i = 1:p.N
                if a == b
                    #direct uptake effects
                    J[a + p.N, a + p.N] -= p.s.σ̃[a] * p.s.γ[i,a] * p.e.fy[i,a]
                end
        #         #uptake effect
                J[a + p.N, b + p.N] += p.s.σ̃[a] * p.s.γ[i,b] * p.tmp[1][i,b]
                
                #release direct effect
                J[a + p.N, b + p.N] += p.s.ρ̃[a] * p.s.η[i,b] * p.e.fy[i,b] * p.s.ϕ[i][b,a] 
                
                #release total effects
                tmp = 0.0
                for c = 1:p.M
                    tmp +=  p.s.ϕ[i][c,a] * (p.tmp[1][i,b] + p.e.λy[i,b])
                end
                J[a + p.N, b + p.N] += p.s.ρ̃[a] * p.s.η[i,a] * tmp
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
            s[i] += abs(v[i,n]) / vals[n]
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
    vals = eigvals(J)
    w = inv(eigvecs(J))
    s = zeros(Complex, L)
    for i = 1:L
        for n = eachindex(vals)
            s[i] += abs(w[i,n]) / vals[n]
        end
    end
    return(log.(-s)) 
end