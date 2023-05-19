module PolynomialSolutions

using LinearAlgebra
using StaticArrays

"""
    struct Polynomial{N,T}

A polynomial in `N` variables with coefficients of type `T`.

The functor interface is implemented, so that `p(x)` evaluates the polynomial.
For performance reasons, `x` is expected to be an `Tuple` (or `SVector`).

# Examples

A polynomial with a single term can be created by passing a pair mapping the
order to the coefficient:

```jldoctest; output = false
julia> Polynomial((1,1)=>2)
2xy
```

When multiple terms are present, they must be passed as vector (or tuple) of pairs:

```jldoctest
julia> Polynomial([(1,1)=>2, (0,1)=>-1])
-y + 2xy
```

The spatial dimension is automatically inferred from the length of the order
tuple:

```jldoctest
julia> Polynomial((1,1,1)=>2)
2xyz
```

"""
struct Polynomial{N,T}
    order2coeff::Dict{NTuple{N,Int},T}
end

# empty constructor
Polynomial{N,T}() where {N,T} = Polynomial{N,T}(Dict{NTuple{N,Int},T}())

# construct a polynomial from a tuple of pairs
Polynomial(t::NTuple{<:Any,Pair{NTuple{N,Int},T}}) where {N,T} = Polynomial{N,T}(Dict(t))

# construct a polynomial from a vector of pairs
Polynomial(v::Vector{Pair{NTuple{N,Int},T}}) where {N,T} = Polynomial{N,T}(Dict(v))

# construct a polynomial from a single pair
Polynomial(p::Pair{NTuple{N,Int},T}) where {N,T} = Polynomial{N,T}(Dict(p))

# functor interface
(p::Polynomial{N})(x::SVector{N}) where {N} = sum(c* prod(x.^θ) for (θ,c) in p.order2coeff)
(p::Polynomial)(x::Tuple) = p(SVector(x))
(p::Polynomial{1})(x::Number) = p(SVector(x))
(p::Polynomial)(x) = error("Expected a tuple or a static vector for argument, got $(typeof(x))")

"""
    is_homogeneous(p::Polynomial)

Return `true` if `p` is homogeneous, i.e. if all the monomials in `p` have the
same degree.
"""
function is_homogeneous(p::Polynomial{N,T}) where {N,T}
    allequal(sum(θ) for θ in keys(p.order2coeff))
end

function Base.iszero(p::Polynomial)
    q = drop_zeros!(deepcopy(p))
    isempty(q.order2coeff)
end

"""
    drop_zeros!(q::Polynomial,tol=0,p=2)

Drop all coefficients in `q` for which the `norm(p,2) ≤ tol`.
"""
function drop_zeros!(q::Polynomial,tol=0,p=2)
    for (k,v) in q.order2coeff
        if norm(v,p) ≤ tol
            delete!(q.order2coeff, k)
        end
    end
    return q
end

"""
    multiply_by_r(p::Polynomial, k::Int = 2)

Multiply a polynomial `p` by the monomial `r^k`, where `r = |𝐱|` and `k` is an
even positive integer.
"""
function multiply_by_r(p::Polynomial{N,T}, k::Int) where {N,T}
    @assert iseven(k)
    k == 0 && return p
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        for d in 1:N
            θ′ = ntuple(i-> i==d ? θ[i] + 2 : θ[i], length(θ))
            order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c
        end
    end
    q = Polynomial(order2coeff)
    return multiply_by_r(q,k-2)
end

"""
    degree(p::Polynomial)

The largest degree of any monomial in `p`.
"""
function degree(p::Polynomial{N,T}) where {N,T}
    deg = 0
    for θ in keys(p.order2coeff)
        deg = max(deg,sum(θ))
    end
    return deg
end

function Base.:+(p1::Polynomial{N,S}, p2::Polynomial{N,T}) where {N,S,T}
    V   = promote_type(S,T)
    acc = Polynomial{N,V}()
    # re-build the p1 elements in the promoted datatype; this is a bit wasteful..
    for (order, coeff) in p1.order2coeff
        acc.order2coeff[order] = coeff
    end
    # loop over the elements of p2. If already in p1, add the coefficients,
    # otherwise add the pair
    for (order, coeff) in p2.order2coeff
        if haskey(acc.order2coeff, order)
            acc.order2coeff[order] += coeff
        else
            acc.order2coeff[order] = coeff
        end
    end
    return acc
end

function Base.:-(p::Polynomial)
    q = deepcopy(p)
    for (order, coeff) in q.order2coeff
        q.order2coeff[order] = - coeff
    end
    return q
end
Base.:-(p1::Polynomial{N,S}, p2::Polynomial{N,T}) where {N,S,T} = p1 + (-p2)

function Base.:(==)(p1::Polynomial{N}, p2::Polynomial{M}) where {N,M}
    N == M ? iszero(p1-p2) : false
end

# multiply a polynomial by a scalar
function Base.:*(c::S, p::Polynomial{N,T}) where {S,N,T}
    V = promote_type(S,T)
    acc = Polynomial{N,V}()
    for (order, coeff) in p.order2coeff
        acc.order2coeff[order] = c * coeff
    end
    return acc
end
Base.:*(p::Polynomial, c) = c * p

"""
    convert_coefs(p::Polynomial, T)

Return a version of `p` where the coefficients have been converted to type `T`
(is such a conversion is possible).
"""
function convert_coefs(p::Polynomial{N,S}, ::Type{T}) where {N,S,T}
    q = Polynomial{N,T}()
    for (order, coeff) in p.order2coeff
        q.order2coeff[order] = T(coeff)
    end
    return q
end

Base.convert(::Type{Polynomial{N,T}}, p::Polynomial{N,S}) where {N,T,S} = convert_coefs(p, T)

"""
    derivative(p::Polynomial, i::Int)

Differentiate `p` with respect to the `i`th variable.
"""
function derivative(p::Polynomial{N,T},d) where {N,T}
    @assert d ∈ 1:N
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        θ[d] < 1 && continue
        θ′ = ntuple(i-> i==d ? θ[d]-1 : θ[i], N)
        c′ = c*(θ[d])
        order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c′
    end
    return Polynomial{N,T}(order2coeff)
end

"""
    gradient(p::Polynomial)

Return an `N`-tuple of the derivatives of `p` with respect to each variable.
"""
function gradient(p::Polynomial{N,T}) where {N,T}
    ntuple(N) do d
        derivative(p,d)
    end |> SVector
end

function laplacian(p::Polynomial{N,T}) where {N,T}
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        for d in 1:N
            θ[d] < 2 && continue
            θ′ = ntuple(i-> i==d ? θ[d]-2 : θ[i], N)
            c′ = c*(θ[d])*(θ[d]-1)
            order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c′
        end
    end
    return Polynomial{N,T}(order2coeff)
end

function divergence(P::SVector{N,Polynomial{N,T}}) where {N,T}
    sum(derivative(P[i],i) for i in 1:N)
end

function curl(P::SVector{N, Polynomial{N,T}}) where {N,T}
    ∇P = gradient.(P)
    if N == 2
        curlP = SVector{3, Polynomial{N,T}}(Polynomial{N,T}(), Polynomial{N,T}(), ∇P[2][1] - ∇P[1][2])
    elseif N == 3
        curlP = SVector{3, Polynomial{N,T}}(∇P[3][2] - ∇P[2][3], ∇P[1][3] - ∇P[3][1], ∇P[2][1] - ∇P[1][2])
    else
        print("Curl not implemented for this value of N")
    end
    return curlP
end

# general show
function Base.show(io::IO, p::Polynomial{N,T}) where {N,T}
    order2coeff = sort(collect(p.order2coeff), by = x -> sum(x[1]))
    isempty(order2coeff) && return print(io, zero(T))
    for (order, coeff) in order2coeff
        # first term is special case
        order == first(order2coeff)[1] || print(io, " + ")
        # print the coefficient
        print(io, "(", coeff, ")")
        # finally print the monomials
        for (i, o) in enumerate(order)
            _print_variable(io, i, o)
        end
    end
end

# adapt show to reals
function Base.show(io::IO, p::Polynomial{N,T}) where {N,T<:Real}
    order2coeff = sort(collect(p.order2coeff), by = x -> sum(x[1]))
    isempty(order2coeff) && return print(io, "0")
    for (order, coeff) in order2coeff
        # first term is special case
        first_coeff = order == first(order2coeff)[1]
        if !first_coeff
            if coeff < 0
                print(io, " - ")
            else
                print(io, " + ")
            end
        else
            if abs(coeff) == 1 && coeff < 0
                print(io, "-")
            end
        end
        # print the coefficient if it is not one
        if sum(order) == 0
            print(io, coeff)
        elseif abs(coeff) != 1
            first_coeff ?  print(io, coeff) : print(io, abs(coeff))
        end
        # finally print the monomials
        for (i, o) in enumerate(order)
            _print_variable(io, i, o)
        end
    end
end

# adapt show to complex
function Base.show(io::IO, p::Polynomial{N,T}) where {N,T<:Complex}
    order2coeff = sort(collect(p.order2coeff), by = x -> sum(x[1]))
    isempty(order2coeff) && return print(io, "0")
    for (order, coeff) in order2coeff
        # first term is special case
        first_coeff = order == first(order2coeff)[1]
        if !first_coeff
            print(io, " + ")
        else
            if coeff == -1
                print(io, "-")
            end
        end
        # print the coefficient if it is not one
        if coeff.im == 0 && coeff != 1
            print(io, coeff.re)
        elseif coeff != 1
            print(io, "(", coeff, ")")
            #print(io, "(", coeff.re, " + ", coeff.im, "ı)")
        end
        # finally print the monomials
        for (i, o) in enumerate(order)
            _print_variable(io, i, o)
        end
    end
end

# verbose code for pretty printing of monomials using unicode
function _print_variable(io, i, p)
    if i == 1
        if p == 0
            print(io, "")
        elseif p == 1
            print(io, "x")
        elseif p == 2
            print(io, "x²")
        elseif p == 3
            print(io, "x³")
        elseif p == 4
            print(io, "x⁴")
        elseif p == 5
            print(io, "x⁵")
        elseif p == 6
            print(io, "x⁶")
        elseif p == 7
            print(io, "x⁷")
        elseif p == 8
            print(io, "x⁸")
        elseif p == 9
            print(io, "x⁹")
        else
            print(io, "x", "^", p)
        end
    elseif i == 2
        if p == 0
            print(io, "")
        elseif p == 1
            print(io, "y")
        elseif p == 2
            print(io, "y²")
        elseif p == 3
            print(io, "y³")
        elseif p == 4
            print(io, "y⁴")
        elseif p == 5
            print(io, "y⁵")
        elseif p == 6
            print(io, "y⁶")
        elseif p == 7
            print(io, "y⁷")
        elseif p == 8
            print(io, "y⁸")
        elseif p == 9
            print(io, "y⁹")
        else
            print(io, "y", "^", p)
        end
    elseif i == 3
        if p == 0
            print(io, "")
        elseif p == 1
            print(io, "z")
        elseif p == 2
            print(io, "z²")
        elseif p == 3
            print(io, "z³")
        elseif p == 4
            print(io, "z⁴")
        elseif p == 5
            print(io, "z⁵")
        elseif p == 6
            print(io, "z⁶")
        elseif p == 7
            print(io, "z⁷")
        elseif p == 8
            print(io, "z⁸")
        elseif p == 9
            print(io, "z⁹")
        else
            print(io, "z", "^", p)
        end
    else
        print(io, "x", i, "^", p)
    end
end

"""
    solve_helmholtz(Q::Polynomial;k=1)

Return the unique polynomial `P` satisfying `ΔP + k²P = Q`.
"""
function solve_helmholtz(Q::Polynomial,k²)
    n = degree(Q)
    m = floor(Int, n/2)
    P = Q
    ΔⁱQ = laplacian(Q)
    for i in 1:m
        P   = P + (-1/k²)^i*ΔⁱQ
        ΔⁱQ = laplacian(ΔⁱQ) # next laplacian
    end
    return 1/k²*P
end
solve_helmholtz(Q::Polynomial;k=1) = solve_helmholtz(Q,k^2)

"""
    solve_laplace(Q::Polynomial)

Return a polynomial `P` satisfying `ΔP = Q`. `Q` is required to be homogeneous.
"""
function solve_laplace(Q::Polynomial{N,T}) where {N,T}
    @assert is_homogeneous(Q) "source term `Q` must be a homogeneous polynomial"
    n = degree(Q)
    γ = (k,p) -> big(2*(k+1)*(2k+2p+N)) # γₖᵖ
    cₖ  = 1//γ(0,n) # c₀
    P   = cₖ*multiply_by_r(deepcopy(Q),2)
    ΔᵏQ = deepcopy(Q)
    m = floor(Int, n/2)
    for k in 1:m
        cₖ  = -cₖ/(γ(k,n-2k))
        ΔᵏQ = laplacian(ΔᵏQ)
        ΔP = cₖ*(multiply_by_r(ΔᵏQ,2k+2))
        P   = P + ΔP
    end
    return P
end

"""
    solve_bilaplace(Q::Polynomial)

Compute a polynomial solution to `Δ²P = Q`. `Q` is required to be homogeneous.
"""
function solve_bilaplace(Q::Polynomial{N}) where {N}
    P′ = solve_laplace(Q)
    P  = solve_laplace(P′)
    return P
end

"""
    solve_stokes(Q::SVector{N,Polynomial{N,T}};μ=1)

Compute a vector of polynomials `U` and a polynomial `P` satisfying `μΔU - ∇P =
Q` with `∇ ⋅ U = 0`.
```
"""
function solve_stokes(Q::SVector{N,Polynomial{N,T}};μ=1//1) where {N,T}
    # u = Δg - ∇ (∇ ⋅ g), p = -μ Δ (∇ ⋅ g), where g solves μΔΔg = Q
    g = 1/μ .* map(q->solve_bilaplace(q),Q)
    h = -divergence(g)
    u = laplacian.(g) .+ gradient(h)
    p = μ*laplacian(h)
    return u,p
end
solve_stokes(Q::NTuple;kwargs...) = solve_stokes(SVector(Q);kwargs...)

"""
    solve_elastodynamics(Q::SVector{N,Polynomial{N,T}};ρ=1,μ=1,ν=1/4,ω=1)

Compute a vector of polynomials `U` satisfying `-μ/(1-2ν) ∇(div U) - μ ΔU - μ
k₂² U = Q`.
"""
function solve_elastodynamics(Q::SVector{N,Polynomial{N,T}};ρ=1//1,μ=1//1,ν=1//4,ω=1//1) where {N,T}
    k₁² = ω^2/(2*μ*(1-ν)/(ρ*(1 - 2ν)))
    k₂² = ω^2*ρ/μ
    g = -1/(2*μ*(1-ν))*map(q->solve_helmholtz(solve_helmholtz(q,k₁²),k₂²),Q)
    u = 2*(1-ν)*(laplacian.(g) + k₁²*g) - gradient(divergence(g))
    return u
end

"""
    solve_elastostatic(Q::SVector{N,Polynomial{N,T}};μ=1,ν=1)

Compute a vector of polynomials `U` satisfying `μ/(1-2ν) ∇(div U) + μΔU = Q`.
"""
function solve_elastostatic(Q::SVector{N, Polynomial{N, T}};μ=1,ν=0) where {N,T}
    g = 1/(2 * μ * (1 - ν)) .* map(q->solve_bilaplace(q), Q)
    u = 2(1 - ν) * laplacian.(g) - gradient(divergence(g))
    return u
end

"""
    solve_maxwell(J::SVector{N,Polynomial{N,T}}, ρ::Polynomial{N, T};ϵ=1,μ=1,ω=1)

Compute a pair of vectors of polynomials `E` and `H` satisfying the Maxwell
system. Also returns the polynomial vector potential `A` and scalar potential
`φ`.
"""
function solve_maxwell(J::SVector{N, Polynomial{N, T}},ρ::Polynomial{N, T};ϵ=1,μ=1,ω=1) where {N,T}
    #@assert divergence(J) - im*ω*ρ == 0
    k² = ω^2 * ϵ * μ
    A = -μ * map(j->solve_helmholtz(j,k²), J)
    φ = -1/ϵ * solve_helmholtz(ρ,k²)
    E = im * ω * A - gradient(φ)
    H = 1/μ * curl(A)
    #return E, H
    return E, H, A, φ
end
function solve_maxwell(J::SVector{N, Polynomial{N, T}},ρ::Polynomial{N, S};kwargs...) where {N,T,S}
    TS = promote_type(T,S)
    solve_maxwell(map(j->convert(Polynomial{N,TS},j),J),convert(Polynomial{N,TS},ρ);kwargs...)
end

export
    Polynomial,
    convert_coefs,
    solve_helmholtz,
    solve_laplace,
    solve_bilaplace,
    solve_stokes,
    solve_elastostatic,
    solve_elastodynamics,
    solve_maxwell

end # module (Polynomials)
