module ElementaryPDESolutionsFixedPolynomialsExt

import ElementaryPDESolutions
import DynamicPolynomials: @polyvar
import FixedPolynomials

function __init__()
    @info "Loading ElementaryPDESolutions.jl FixedPolynomials extension"
end

struct PolyFastEvaluator{N,S,T}
    N::Int64
    PolSystem::FixedPolynomials.System{S}
    cfg::FixedPolynomials.JacobianConfig{T,S}
end

function PolyFastEvaluator(N, PolSystem::FixedPolynomials.System{S}, cfg::FixedPolynomials.JacobianConfig{T,S}) where {S,T}
    return PolyFastEvaluator{N,S,T}(N, PolSystem, cfg)
end

function ElementaryPDESolutions.fast_evaluate_with_jacobian!(vals::Array{S}, grad::Array{T}, x::Vector{Vector{S}},
                                PFE::PolyFastEvaluator{N,S,T}) where {N,S,T}
    @assert size(vals) == (length(PFE.PolSystem), length(x))
    @assert size(grad) == (length(PFE.PolSystem), N, length(x))
    for i in 1:length(x)
        FixedPolynomials.evaluate_and_jacobian!(view(vals, :, i), view(grad, :, :, i), PFE.PolSystem, x[i],
                                PFE.cfg)
    end
end

function ElementaryPDESolutions.fast_evaluate!(vals::Array{S}, x::Vector{Vector{S}},
                                PFE::PolyFastEvaluator{N,S,T}) where {N,S,T}
    @assert size(vals) == (length(PFE.PolSystem), length(x))
    for i in 1:length(x)
        FixedPolynomials.evaluate!(view(vals, :, i), PFE.PolSystem, x[i], PFE.cfg)
    end
end

function ElementaryPDESolutions.assemble_fastevaluator(
        Pols::Vector{ElementaryPDESolutions.Polynomial{N,T}}
    ) where {N, T}

    if N == 1
        @polyvar x
    elseif N == 2
        @polyvar x y
    elseif N == 3
        @polyvar x y z
    else
        error("assemble_fastevaluator only implemented for dimensions ≤ 3")
    end

    PolArray = Array{FixedPolynomials.Polynomial{Float64}}(undef, length(Pols))
    for (polind, poly) in enumerate(Pols)
        coeffs = poly.order2coeff
        exp_data = Matrix{Int64}(undef, N, length(coeffs))
        coeff_data = Vector{Float64}(undef, length(coeffs))
        for (i, pol) in enumerate(coeffs)
            exp_data[:, i] = [q for q in pol[1]]
            coeff_data[i] = pol[2]
        end
        if N == 1
            PolArray[polind] = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x])
        elseif N == 2
            PolArray[polind] = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x, :y])
        else
            PolArray[polind] = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x, :y, :z])
        end
    end
    PolSystem = FixedPolynomials.System(PolArray)
    pts = Vector{Vector{Float64}}(undef, 1)
    pts[1] = zeros(N) 
    cfg = FixedPolynomials.JacobianConfig(PolSystem, pts[1])
    return PolyFastEvaluator(N, PolSystem, cfg)
end

end
