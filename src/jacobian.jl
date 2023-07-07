#jacobian

function jacobian!(p::Parameters, J)
    for i = 1:p.N
        #loop over consumers
        for j = 1:p.N
            #diagonal
            if i == j
                J[i,j] = -p.ω[i]
                for a = 1:p.M
                    J[i,j] += p.β[i,a] * p.μ_x[i,a]
                    for b = 1:p.M
                        J[i,j] -= p.β[i,a] * p.T[i,a] * p.ψ[i,a,b] * (p.μ_x[i,a] + p.ζ_x[i,a,b])
                    end
                end
            else
                J[i,j] = 0
            end
            J[i,j] *= p.a_x[i]
        end

        #loop over resources
        for c = 1:p.M
            J[i,p.N + c] = p.β[i,c] * p.μ_x[i,c]

            for a = 1:p.M
                J[i,p.N + c] -= p.β[i,a] * p.T[i,a] * p.ψ[i,a,c] * p.ζ_P[i,a,c]
            end

            for b = 1:p.M
                J[i,p.N + c] -= p.β[i,c] * p.T[i,c] * p.ψ[i,c,b] * (p.μ_y[i,c] + p.ζ_S[i,c,b])
            end
            J[i,p.N + c] *= p.a_x[i]
        end
        
        
    end

    #loop over resources
    for a = 1:p.M
        #over consumers
        for j = 1:p.N
            J[a + p.N, j] = -p.σ̃[a] * p.γ[j,a] * p.μ_x[j,a]
            for b = 1:p.M
                J[a + p.N, j] -= p.ρ̃[a] * p.ϕ[j,b,a] * (p.μ_x[j,a] + p.ζ_x[j,b,a])
            end
            J[p.N + a, j] *= p.a_y[a]
        end

        #over resources
        for c = 1:p.M
            if a == c
                J[a + p.N, c + p.N] = p.ρ[a] * p.ι[a] - p.σ[a] * p.ϱ[a]
                    for i = 1:p.N
                        J[a + p.N, c + p.N] -= p.σ̃[a] * p.γ[i,a] * p.μ_y[i,a] 
                        for b = 1:p.M
                            J[a + p.N, c + p.N] += p.ρ̃[a] * p.ϕ[i,b,a] * p.ζ_P[i,b,a]
                        end
                    end
            else
                J[a + p.N, c + p.N] = 0.0
                for i = 1:p.N
                    J[a + p.N, c + p.N] += p.ρ̃[a] * p.ϕ[i,c,a] * (p.μ_y[i,c] + p.ζ_S[i,c,a])
                end
            end
            J[p.N + a, p.N + c] *= p.a_y[a]
        end
    end
end