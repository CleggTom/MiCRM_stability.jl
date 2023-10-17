#jacobian
function jacobian!(p::Parameters, J)
    #loop over consumers
    for i = 1:p.N
        for j = 1:p.N
            if i == j
                J[i,i] = p.e.gx[i] - p.e.mx[i]
            else
                J[i,j] = 0
            end
            
            J[i,j] *= p.s.α[i]
        end

        for a = 1:p.M
            J[i, a + p.M] = p.e.gs[i] * p.s.χ[i,a] * p.e.fy[i,a] + p.e.gw[i] * p.s.ϕ[i,a] * p.e.hy[i,a] + p.e.ωy[i,a]
            
            J[i, a + p.M] *= p.s.α[i]
        end

    end

    for a = 1:p.M
        for i = 1:p.N
            J[a + p.M, i] = p.e.gx[i] * (p.s.ρ̃[a] * p.s.η[i,a] - p.s.σ̃[a] * p.s.γ[i,a])
        end

        for b = 1:p.M
            if a == b
                J[a + p.M, a + p.M] = p.s.ρ[a] * p.e.iy[a] - p.s.σ[a] * p.e.oy[a]

                for i = 1:p.N
                    J[a + p.M, a + p.M] -= p.s.σ̃[a] * p.s.γ[i,a] * (p.e.fy[i,a] * (p.s.χ[i,a] * (p.e.gs[i] - 1) +1) + 
                    p.e.gw[i] * p.s.ϕ[i,a] * p.e.hy[i,a])

                    J[a + p.M, a + p.M] += p.s.ρ̃[a] * p.s.η[i,a] * (p.e.gs[i] * p.s.χ[i,a] * p.e.fy[i,a] + p.e.gw[i] * p.s.ϕ[i,a] * p.e.hy[i,a] + p.e.λy[i,a] + p.e.θy[i,a,a])
                end

            else 
                J[a + p.M, b + p.M] = 0
                
                for i = 1:p.N
                    J[a + p.M, b + p.M] -= p.s.σ̃[a] * p.s.γ[i,a] * (p.e.fy[i,b] * p.s.χ[i,b] * (p.e.gs[i] - 1) + p.e.gw[i] * p.s.ϕ[i,b] * p.e.hy[i,b])

                    J[a + p.M, b + p.M] += p.s.ρ̃[a] * p.s.η[i,a] * (p.e.gs[i] * p.s.χ[i,b] * p.e.fy[i,b] + p.e.gw[i] * p.s.ϕ[i,b] * p.e.hy[i,b] + p.e.λy[i,b] + p.e.θy[i,a,b])
                end

            end
        end
    end
end