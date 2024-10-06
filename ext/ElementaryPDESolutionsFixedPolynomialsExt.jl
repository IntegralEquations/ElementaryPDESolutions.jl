module ElementaryPDESolutionsFixedPolynomialsExt

using ElementaryPDESolutions: ElementaryPDESolutions
import DynamicPolynomials: @polyvar
using FixedPolynomials: FixedPolynomials

function __init__()
    @info "Loading ElementaryPDESolutions.jl FixedPolynomials extension"
end

# point evaluations at points x::AbstractVector{S} of vectors of polynomials F::Vector{Polynomial{T}}
struct VecPolyFastEvaluator{N,S,T}
    N::Int64
    PolSystem::FixedPolynomials.System{S}
    cfg::FixedPolynomials.JacobianConfig{T,S}
end

# point evaluations at points x::AbstractVector{S} of a polynomial F::Polynomial{T}
struct ScaPolyFastEvaluator{N,T}
    N::Int64
    Pol::FixedPolynomials.Polynomial{T}
    r::FixedPolynomials.GradientDiffResult{T,Vector{T}}
    cfg::FixedPolynomials.GradientConfig{T}
end

function VecPolyFastEvaluator(N, PolSystem::FixedPolynomials.System{S},
                              cfg::FixedPolynomials.JacobianConfig{T,S}) where {S,T}
    return VecPolyFastEvaluator{N,S,T}(N, PolSystem, cfg)
end

function ScaPolyFastEvaluator(N, Pol, r::FixedPolynomials.GradientDiffResult{T},
                              cfg::FixedPolynomials.GradientConfig{T}) where {T}
    return ScaPolyFastEvaluator{N,T}(N, Pol, r, cfg)
end

Base.length(VPFE::VecPolyFastEvaluator{N,S,T}) where {N,S,T} = length(VPFE.PolSystem)

function ElementaryPDESolutions.fast_evaluate_with_jacobian!(vals::AbstractArray{S},
                                                             grad::AbstractArray{S},
                                                             x::AbstractVector{S},
                                                             VPFE::VecPolyFastEvaluator{N,S,
                                                                                        T}) where {N,
                                                                                                   S,
                                                                                                   T}
    @assert length(vals) == length(VPFE.PolSystem)
    @assert size(grad) == (length(VPFE.PolSystem), N)
    return FixedPolynomials.evaluate_and_jacobian!(vals, grad, VPFE.PolSystem, x,
                                                   VPFE.cfg)
end

function ElementaryPDESolutions.fast_evaluate_with_gradient!(vals::AbstractArray{S},
                                                             grad::AbstractArray{S},
                                                             x::AbstractVector{S},
                                                             SPFE::ScaPolyFastEvaluator{N,
                                                                                        T}) where {N,
                                                                                                   S,
                                                                                                   T}
    @assert SPFE.N == length(x)
    @assert length(grad) == N
    FixedPolynomials.gradient!(SPFE.r, SPFE.Pol, x, SPFE.cfg)
    vals = FixedPolynomials.value(SPFE.r)
    return grad .= FixedPolynomials.gradient(SPFE.r)
end

function ElementaryPDESolutions.fast_evaluate!(vals::AbstractArray{S}, x::AbstractVector{S},
                                               VPFE::VecPolyFastEvaluator{N,S,T}) where {N,
                                                                                         S,
                                                                                         T}
    @assert VPFE.N == length(x)
    @assert length(vals) == length(VPFE.PolSystem)
    return FixedPolynomials.evaluate!(vals, VPFE.PolSystem, x, VPFE.cfg)
end

function ElementaryPDESolutions.fast_evaluate!(vals::AbstractArray{S}, x::AbstractVector{S},
                                               SPFE::ScaPolyFastEvaluator{N,T}) where {N,S,
                                                                                       T}
    @assert SPFE.N == length(x)
    return vals = FixedPolynomials.evaluate(SPFE.Pol, x, SPFE.cfg)
end

"""
    assemble_fastevaluator(Pols::Vector{ElementaryPDESolutions.Polynomial{N,T}}, ::Type{S})

Perform precomputations for fast evaluation of a `Vector` of `N`-variate
`Polynomial` of input type `S`, with coefficients of type T.  Return a
`VecPolyFastEvaluator` object for online use.
"""
function ElementaryPDESolutions.assemble_fastevaluator(Pols::Vector{ElementaryPDESolutions.Polynomial{N,
                                                                                                      T}},
                                                       ::Type{S}) where {N,S,T}
    if N == 1
        @polyvar x
    elseif N == 2
        @polyvar x y
    elseif N == 3
        @polyvar x y z
    else
        # implementing for dimensions > 3 not hard if desired, just need to
        # modify the Polynomial() call with new variable names
        error("assemble_fastevaluator only implemented for dimensions ≤ 3")
    end

    PolArray = Array{FixedPolynomials.Polynomial{T}}(undef, length(Pols))
    for (polind, poly) in enumerate(Pols)
        coeffs = poly.order2coeff
        exp_data = Matrix{Int64}(undef, N, length(coeffs))
        coeff_data = Vector{T}(undef, length(coeffs))
        for (i, pol) in enumerate(coeffs)
            exp_data[:, i] = [q for q in pol[1]]
            coeff_data[i] = pol[2]
        end
        if N == 1
            PolArray[polind] = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x])
        elseif N == 2
            PolArray[polind] = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x, :y])
        else
            PolArray[polind] = FixedPolynomials.Polynomial(exp_data, coeff_data,
                                                           [:x, :y, :z])
        end
    end
    PolSystem = FixedPolynomials.System(PolArray)
    cfg = FixedPolynomials.JacobianConfig(PolSystem, S)
    return VecPolyFastEvaluator(N, PolSystem, cfg)
end

"""
    assemble_fastevaluator(Pols::Vector{ElementaryPDESolutions.Polynomial{N,T}}, ::Type{S})

Perform precomputations for fast evaluation of a `Polynomial` on points of type `S`.
Return a `ScaPolyFastEvaluator` object for online use.
"""
function ElementaryPDESolutions.assemble_fastevaluator(Pol::ElementaryPDESolutions.Polynomial{N,
                                                                                              T},
                                                       ::Type{S}) where {N,S,T}
    if N == 1
        @polyvar x
    elseif N == 2
        @polyvar x y
    elseif N == 3
        @polyvar x y z
    else
        # implementing for dimensions > 3 not hard if desired, just need to
        # modify the Polynomial() call with new variable names
        error("assemble_fastevaluator only implemented for dimensions ≤ 3")
    end

    coeffs = Pol.order2coeff
    exp_data = Matrix{Int64}(undef, N, length(coeffs))
    coeff_data = Vector{Float64}(undef, length(coeffs))
    for (i, pol) in enumerate(coeffs)
        exp_data[:, i] = [q for q in pol[1]]
        coeff_data[i] = pol[2]
    end
    if N == 1
        Pol = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x])
    elseif N == 2
        Pol = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x, :y])
    else
        Pol = FixedPolynomials.Polynomial(exp_data, coeff_data, [:x, :y, :z])
    end
    cfg = FixedPolynomials.GradientConfig(Pol, S)
    r = FixedPolynomials.GradientDiffResult(cfg)
    return ScaPolyFastEvaluator(N, Pol, r, cfg)
end

end
