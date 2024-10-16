module ElementaryPDESolutions
using LinearAlgebra

"""
    struct Polynomial{N,T}

A polynomial in `N` variables with coefficients of type `T`.

The functor interface is implemented, so that `p(x)` evaluates the polynomial.
For performance reasons, `x` is expected to be a `Tuple`.

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
function (p::Polynomial{N,T})(x) where {N,T}
    @assert length(x) == N "Expected input of length $N, got $(length(x))"
    return sum(c * prod(x .^ θ) for (θ, c) in p.order2coeff; init=zero(T))
end
(p::Polynomial)(x...) = p(x) # so that e.g. p(x,y) works

"""
    is_homogeneous(p::Polynomial)

Return `true` if `p` is homogeneous, i.e. if all the monomials in `p` have the
same degree.
"""
function is_homogeneous(p::Polynomial{N,T}) where {N,T}
    return allequal(sum(θ) for θ in keys(p.order2coeff))
end

function Base.iszero(p::Polynomial)
    q = drop_zeros!(deepcopy(p))
    return isempty(q.order2coeff)
end

"""
    drop_zeros!(q::Polynomial,tol=0,p=2)

Drop all coefficients in `q` for which the `abs(p) ≤ tol`.
"""
function drop_zeros!(q::Polynomial, tol=0)
    for (k, v) in q.order2coeff
        if abs(v) ≤ tol
            delete!(q.order2coeff, k)
        end
    end
    return q
end

"""
    multiply_by_r(p::Polynomial, k::Int = 2)

Multiply a polynomial `p` by the polynomial `r^k`, where `r = |𝐱|` and `k` is an
even positive integer.
"""
function multiply_by_r(p::Polynomial{N,T}, k::Int) where {N,T}
    @assert iseven(k)
    k == 0 && return p
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        for d in 1:N
            θ′ = ntuple(i -> i == d ? θ[i] + 2 : θ[i], length(θ))
            order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c
        end
    end
    q = Polynomial(order2coeff)
    return multiply_by_r(q, k - 2)
end

"""
    multiply_by_anisotropic_anisotropic_r(A::AbstractMatrix{T}, p::Polynomial, k::Int = 2)

Multiply a polynomial `p` by the polynomial `r_A^k`, where `r_A = |r^T A^{-1} r|`,
r = (x_1, x_2, ..., x_n), and `k` is an even positive integer.
"""
function multiply_by_anisotropic_r(A::AbstractMatrix{T}, p::Polynomial{N,T},
                                   k::Int) where {N,T}
    @assert LinearAlgebra.checksquare(A) == N
    @assert iseven(k)
    # This slows us down, but prevents a degradation relative to Base.LinearAlgebra
    # when `using StaticArrays`
    # https://github.com/JuliaArrays/StaticArrays.jl/issues/434
    @assert det(A) ≠ zero(T) "anisotropic tensor must be invertible"
    k == 0 && return p
    order2coeff = empty(p.order2coeff)
    invA = inv(A)
    for (θ, c) in p.order2coeff
        for i in 1:N
            for j in 1:N
                θ′ = ntuple(l -> θ[l] + Int(l == j) + Int(l == i), length(θ))
                order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c * invA[i, j]
            end
        end
    end
    q = Polynomial(order2coeff)
    return multiply_by_anisotropic_r(A, q, k - 2)
end

"""
    multiply_by_anisotropic_β_r(β::AbstractVector{T}, p::Polynomial, k::Int)

Multiply a polynomial `p` by the polynomial (β ⋅ 𝐫)ᵏ, 𝐫 = (x_1, x_2, ..., x_n),
and `k` is a non-negative integer.
"""
function multiply_by_anisotropic_β_r(β::AbstractVector, p::Polynomial{N,T},
                                     k::Int) where {N,T}
    @assert length(β) == N
    @assert k ≥ 0
    k == 0 && return p
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        for i in 1:N
            θ′ = ntuple(l -> θ[l] + Int(l == i), length(θ))
            order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c * β[i]
        end
    end
    return multiply_by_anisotropic_β_r(β, Polynomial(order2coeff), k - 1)
end

"""
    degree(p::Polynomial)

The largest degree of any monomial in `p`.
"""
function degree(p::Polynomial{N,T}) where {N,T}
    deg = 0
    for θ in keys(p.order2coeff)
        deg = max(deg, sum(θ))
    end
    return deg
end

function Base.:+(p1::Polynomial{N,S}, p2::Polynomial{N,T}) where {N,S,T}
    V = promote_type(S, T)
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
        q.order2coeff[order] = -coeff
    end
    return q
end
Base.:-(p1::Polynomial{N,S}, p2::Polynomial{N,T}) where {N,S,T} = p1 + (-p2)

function Base.:(==)(p1::Polynomial{N}, p2::Polynomial{M}) where {N,M}
    return N == M ? iszero(p1 - p2) : false
end

# multiply a polynomial by a scalar
function Base.:*(c::Number, p::Polynomial{N,T}) where {N,T}
    S = typeof(c)
    V = promote_type(S, T)
    acc = Polynomial{N,V}()
    for (order, coeff) in p.order2coeff
        acc.order2coeff[order] = c * coeff
    end
    return acc
end
Base.:*(p::Polynomial, c::Number) = c * p

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

function Base.convert(::Type{Polynomial{N,T}}, p::Polynomial{N,S}) where {N,T,S}
    return convert_coefs(p, T)
end

"""
    derivative(p::Polynomial, i::Int)

Differentiate `p` with respect to the `i`th variable.
"""
function derivative(p::Polynomial{N,T}, d) where {N,T}
    @assert d ∈ 1:N
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        θ[d] < 1 && continue
        θ′ = ntuple(i -> i == d ? θ[d] - 1 : θ[i], N)
        c′ = c * (θ[d])
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
        return derivative(p, d)
    end
end

function laplacian(p::Polynomial{N,T}) where {N,T}
    order2coeff = empty(p.order2coeff)
    for (θ, c) in p.order2coeff
        for d in 1:N
            θ[d] < 2 && continue
            θ′ = ntuple(i -> i == d ? θ[d] - 2 : θ[i], N)
            c′ = c * (θ[d]) * (θ[d] - 1)
            order2coeff[θ′] = get(order2coeff, θ′, zero(T)) + c′
        end
    end
    return Polynomial{N,T}(order2coeff)
end

"""
    anisotropic_laplacian(A::AbstractMatrix, P::Polynomial)

Evaluate the anisotropic Laplacian `∇ ⋅ (A ∇P)`.
"""
function anisotropic_laplacian(A::AbstractMatrix, p::Polynomial{N}) where {N}
    @assert LinearAlgebra.checksquare(A) == N
    ∇p = gradient(p)
    Δp = sum(derivative(sum(A[i, j] * ∇p[j] for j in 1:N), i) for i in 1:N)
    return Δp
end

function divergence(P::NTuple{N,Polynomial{N,T}}) where {N,T}
    return sum(derivative(P[i], i) for i in 1:N)
end

function curl(P::NTuple{N,Polynomial{N,T}}) where {N,T}
    ∇P = gradient.(P)
    if N == 2
        curlP = (Polynomial{N,T}(), Polynomial{N,T}(), ∇P[2][1] - ∇P[1][2])
    elseif N == 3
        curlP = (∇P[3][2] - ∇P[2][3], ∇P[1][3] - ∇P[3][1], ∇P[2][1] - ∇P[1][2])
    else
        print("Curl not implemented for N = $N")
    end
    return curlP
end

# general show
function Base.show(io::IO, p::Polynomial{N,T}) where {N,T}
    order2coeff = sort(collect(p.order2coeff); by=x -> sum(x[1]))
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
    order2coeff = sort(collect(p.order2coeff); by=x -> sum(x[1]))
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
            first_coeff ? print(io, abs(coeff)) : print(io, coeff)
        elseif abs(coeff) != 1
            first_coeff ? print(io, coeff) : print(io, abs(coeff))
        end
        # finally print the monomials
        for (i, o) in enumerate(order)
            _print_variable(io, i, o)
        end
    end
end

# adapt show to complex
function Base.show(io::IO, p::Polynomial{N,T}) where {N,T<:Complex}
    order2coeff = sort(collect(p.order2coeff); by=x -> sum(x[1]))
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
        print(io, "x_", i, "^", p)
    end
end

"""
    solve_helmholtz(Q::Polynomial;k=1)

Return the unique polynomial `P` satisfying `ΔP + k²P = Q`.

# Examples

```jldoctest
julia> Q = Polynomial((1,2)=>1)
xy²

julia> P = solve_helmholtz(Q, k=1)
-2.0x + xy²
```
"""
function solve_helmholtz(Q::Polynomial, k²)
    n = degree(Q)
    m = floor(Int, n / 2)
    P = Q
    ΔⁱQ = laplacian(Q)
    for i in 1:m
        P = P + (-1 / k²)^i * ΔⁱQ
        ΔⁱQ = laplacian(ΔⁱQ) # next laplacian
    end
    return 1 / k² * P
end
solve_helmholtz(Q::Polynomial; k=1) = solve_helmholtz(Q, k^2)

"""
    solve_laplace(Q::Polynomial)

Return a polynomial `P` satisfying `ΔP = Q`. `Q` is required to be homogeneous.

# Examples

```jldoctest
julia> Q = Polynomial((1,0)=>1.0)
x

julia> P = solve_laplace(Q)
0.125xy² + 0.125x³
```
"""
function solve_laplace(Q::Polynomial{N,T}) where {N,T}
    @assert is_homogeneous(Q) "source term `Q` must be a homogeneous polynomial"
    n = degree(Q)
    γ = (k, p) -> 2 * (k + 1) * (2k + 2p + N) # γₖᵖ
    cₖ = big(1) // γ(0, n) # c₀
    P = cₖ * multiply_by_r(deepcopy(Q), 2)
    ΔᵏQ = deepcopy(Q)
    m = floor(Int, n / 2)
    for k in 1:m
        cₖ = -cₖ / (γ(k, n - 2k))
        ΔᵏQ = laplacian(ΔᵏQ)
        ΔP = cₖ * (multiply_by_r(ΔᵏQ, 2k + 2))
        P = P + ΔP
    end
    return P
end

"""
    solve_anisotropic_laplace(A::AbstractMatrix{T}, Q::Polynomial)

Return a polynomial `P` satisfying the anisotropic Laplace equation `∇ ⋅ (A ∇P) = Q`, `A` an invertible
matrix. `Q` is required to be homogeneous. Inverse is [`anisotropic_laplacian`](@ref).

# Examples

```jldoctest
using StaticArrays
A = SMatrix{2,2,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, 3 // 1)
Q = Polynomial([(1, 1) => 2 // 1])
P = solve_anisotropic_laplace(A, Q)

# output

-3//400x⁴ + 11//100x³y + 11//150xy³ - 2//25x²y² - 1//300y⁴
```
"""
function solve_anisotropic_laplace(A::AbstractMatrix{T}, Q::Polynomial{N,T}) where {N,T}
    @assert LinearAlgebra.checksquare(A) == N
    @assert A == transpose(A) "anisotropic tensor must be symmetric"
    @assert is_homogeneous(Q) "source term `Q` must be a homogeneous polynomial"

    n = degree(Q)
    γ = (k, p) -> 2 * (k + 1) * (2k + 2p + N) # γₖᵖ
    cₖ = big(1) // γ(0, n) # c₀
    P = cₖ * multiply_by_anisotropic_r(A, deepcopy(Q), 2)
    ΔᵏQ = deepcopy(Q)
    m = floor(Int, n / 2)
    for k in 1:m
        cₖ = -cₖ / (γ(k, n - 2k))
        ΔᵏQ = anisotropic_laplacian(A, ΔᵏQ)
        ΔP = cₖ * (multiply_by_anisotropic_r(A, ΔᵏQ, 2k + 2))
        P = P + ΔP
    end
    return P
end

"""
    solve_anisotropic_advect_diffuse(A::SMatrix{N, N}, β::AbstractVector{T}, Q::Polynomial)

Return a polynomial `P` satisfying the anisotropic advection-diffusion equation
`∇ ⋅ (A ∇P) + β⋅∇P = Q`, `A` a symmetric positive definite matrix.

# Examples

```jldoctest
using StaticArrays
A = SMatrix{2,2,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, 3 // 1)
β = SVector{2,Rational{Int64}}(2 // 1, 1 // 1)
Q = Polynomial([(0, 1) => 2 // 1])
P = solve_anisotropic_advect_diffuse(A, β, Q)

# output

-14//25y - 28//25x + 16//25xy + 9//25y² - 4//25x²
```
"""
function solve_anisotropic_advect_diffuse(A::AbstractMatrix, β::AbstractVector,
                                          Q::Polynomial{N,T}) where {N,T}
    @assert length(β) == N "β must be dimensionally consistent with Q"

    n = degree(Q)
    uᵢ = solve_anisotropic_advect(β, deepcopy(Q))
    P = Polynomial{N,T}() + uᵢ
    for i in (n - 1):-1:0
        uᵢ = solve_anisotropic_advect(β, -anisotropic_laplacian(A, uᵢ))
        P = P + uᵢ
    end
    return P
end

"""
    solve_anisotropic_advect(β::AbstractVector, Q::Polynomial)

Return a polynomial `P` satisfying the anisotropic advection equation `β⋅∇P = Q`.

# Examples

```jldoctest
using StaticArrays
β = SVector{2,Rational{Int64}}(2 // 1, 1 // 1)
Q = Polynomial([(0, 0) => 2 // 1])
P = solve_anisotropic_advect(β, Q)

# output

2//5y + 4//5x
```
"""
function solve_anisotropic_advect(β::AbstractVector, Q::Polynomial{N,T}) where {N,T}
    @assert length(β) == N "β must be dimensionally consistent with Q"

    n = degree(Q)
    betagradellq = deepcopy(Q)
    cₗ = big(1) # c₀
    P = cₗ * multiply_by_anisotropic_β_r(β, Q, 1)
    β2 = sum(β[i]^2 for i in 1:N)
    for l in 1:n
        cₗ = -cₗ / ((l + 1) * β2)
        betagradellq = sum(β[i] * gradient(betagradellq)[i] for i in 1:N)
        P = P + cₗ * multiply_by_anisotropic_β_r(β, betagradellq, l + 1)
    end
    return (1 / β2) * P
end

"""
    solve_bilaplace(Q::Polynomial)

Compute a polynomial solution to `Δ²P = Q`. `Q` is required to be homogeneous.

# Examples

```jldoctest
julia> Q = Polynomial((1,0)=>1)
x

julia> P = solve_bilaplace(Q)
1//192x⁵ + 1//96x³y² + 1//192xy⁴
```
"""
function solve_bilaplace(Q::Polynomial{N}) where {N}
    P′ = solve_laplace(Q)
    P = solve_laplace(P′)
    return P
end

"""
    solve_stokes(Q::NTuple{N,Polynomial{N,T}};μ=1)

Compute a vector of polynomials `U` and a polynomial `P` satisfying `μΔU - ∇P =
Q` with `∇ ⋅ U = 0`. `Q` is required to be homogeneous.

# Examples

```jldoctest
julia> Q = (Polynomial((1,0)=>1),Polynomial((0,0)=>1))
(x, 1)

julia> P = solve_stokes(Q;μ=Rational(1))
((-1//8xy + 1//16xy² + 1//48x³, 3//16x² + 1//16y² - 1//48y³ - 1//16x²y), -1//2y - 3//8x² - 1//8y²)
```
"""
function solve_stokes(Q::NTuple{N,Polynomial{N,T}}; μ=1 // 1) where {N,T}
    # u = Δg - ∇ (∇ ⋅ g), p = -μ Δ (∇ ⋅ g), where g solves μΔΔg = Q
    g = 1 / μ .* map(q -> solve_bilaplace(q), Q)
    h = -divergence(g)
    u = laplacian.(g) .+ gradient(h)
    p = μ * laplacian(h)
    return u, p
end

"""
    solve_brinkman(Q::NTuple{N,Polynomial{N,T}};Re=1,α=1)

Compute a vector of polynomials `U` and a polynomial `P` satisfying the
linearized unsteady Navier-Stokes equations, sometimes referred to as the Brinkman equations
or the modified Stokes equations, `(Δ - α²)U - Re ∇P = Q` with `∇⋅U = 0`. Each component of the
polynomial `Q` is required to be individually homogeneous.

The solutions are given by the expressions

    u = (Δ + α²)(Δ - ∇∇⋅)g,

    p = -1/Re (Δ² - α⁴)∇⋅g,

where the vector potential g satisfies

    (Δ³ - α⁴Δ)g = Q.

# Examples

```jldoctest
julia> Q = (Polynomial([(2, 1) => 2 // 1]), Polynomial([(0, 2) => 4 // 1]))
(2//1x²y, 4//1y²)

julia> U, P = solve_brinkman(Q; Re=Rational(1), α=Rational(1))
((0//1y + xy + 5//24y³ - 5//8x²y, 0//1 + 4//1x + 1//2x² - 1//2y² + 5//8xy² + 11//24x³), -7//6y³ - 1//2x²y - 11//24x³y - 5//24xy³)
```
"""
function solve_brinkman(Q::NTuple{N,Polynomial{N,T}}; Re=1 // 1, α=1 // 1) where {N,T}
    g = brinkman_component_solver.(Q, α)
    divg = divergence(g)
    v = laplacian.(g) .- gradient(divg)

    U = laplacian.(v) .+ α^2 .* v
    P = -1 / Re * (laplacian(laplacian(divg)) - α^4 * divg)

    return U, P
end

"""
    brinkman_component_solver(Q::Polynomial{N,T}, α) where {N,T}

Compute a polynomial vector potential `P` satisfying the auxiliary vector PDE.

    (Δ³ - α⁴Δ)P = Q

for the Brinkman (linearized Navier-Stokes) system.
"""
function brinkman_component_solver(Q::Polynomial{N,T}, α) where {N,T}
    n = degree(Q)
    m = cld(n + 1, 4) - 1 # q = 2, r = 6 in paper
    uᵢ = -1 / α^4 * solve_laplace(deepcopy(Q))
    P = Polynomial{N,T}() + uᵢ
    for _ in 0:(m - 1)
        uᵢ = -1 / α^4 * solve_laplace(-laplacian(laplacian(laplacian(uᵢ))))
        P = P + uᵢ
    end
    return P
end

"""
    solve_elastodynamics(Q::NTuple{N,Polynomial{N,T}};ρ=1,μ=1,ν=1/4,ω=1)

Compute a vector of polynomials `U` satisfying `-μ/(1-2ν) ∇(div U) - μ ΔU - μ
k₂² U = Q`.

# Examples

```jldoctest
julia> Q = (Polynomial((2,1)=>1),Polynomial((1,0)=>1))
(x²y, x)

julia> P = solve_elastodynamics(Q;μ=1)
(-6//1y + x²y, -3//1x)
```
"""
function solve_elastodynamics(Q::NTuple{N,Polynomial{N,T}}; ρ=1 // 1, μ=1 // 1, ν=1 // 4,
                              ω=1 // 1) where {N,T}
    k₁² = ω^2 / (2 * μ * (1 - ν) / (ρ * (1 - 2ν)))
    k₂² = ω^2 * ρ / μ
    g = -1 / (2 * μ * (1 - ν)) .* map(q -> solve_helmholtz(solve_helmholtz(q, k₁²), k₂²), Q)
    u = -2 * (1 - ν) .* (laplacian.(g) .+ k₁² .* g) .+ gradient(divergence(g))
    return u
end

"""
    solve_elastostatic(Q::NTuple{N,Polynomial{N,T}};μ=1,ν=1)

Compute a vector of polynomials `U` satisfying `μ/(1-2ν) ∇(div U) + μΔU = Q`. `Q` is required to be homogeneous.

# Examples

```jldoctest
julia> Q = (Polynomial((1,2)=>1), Polynomial((0,0)=>1))
(xy², 1)

julia> P = solve_elastostatic(Q;ν=1//2)
(-1//8xy + 1//480x⁵ + 1//32x³y² + 1//24xy⁴, 3//16x² + 1//16y² - 1//120y⁵ - 1//96x⁴y - 1//32x²y³)
```
"""
function solve_elastostatic(Q::NTuple{N,Polynomial{N,T}}; μ=1, ν=0) where {N,T}
    g = 1 / (2 * μ * (1 - ν)) .* map(q -> solve_bilaplace(q), Q)
    u = 2(1 - ν) .* laplacian.(g) .- gradient(divergence(g))
    return u
end

@doc raw"""
    solve_maxwell(J::NTuple{3,Polynomial{3,T}};ϵ=1,μ=1,ω=1)

Compute a pair of vectors of polynomials `E` and `H` satisfying the Maxwell
system:

```math
\begin{aligned}
  \mathrm{i}\omega\varepsilon\boldsymbol{E} + \operatorname{rot} \boldsymbol{H} &= \boldsymbol{J}, \qquad &
  -\mathrm{i}\omega\mu\boldsymbol{H} + \operatorname{rot}\boldsymbol{E} &= \boldsymbol{0}, \\
  \varepsilon\operatorname{div}\boldsymbol{E} &= \rho, &
  \mu\operatorname{div}\boldsymbol{H} &= 0,
\end{aligned}
```

with the sources being constrained by the charge conservation equation:

```math
\begin{aligned}
  \operatorname{div}\boldsymbol{J} - \mathrm{i}\omega\rho &= 0.
\end{aligned}
```

Returns the pair `(E, H)`.

# Examples

```jldoctest
julia> J = (Polynomial((0,2,1) => 1), Polynomial((0,1,1) => 1), Polynomial((1,0,1) => 1))
(y²z, yz, xz)

julia> E, H = solve_maxwell(J);

julia> E
((-0.0 - 1.0im) + (0.0 + 2.0im)z + (-0.0 - 1.0im)y²z, (-0.0 - 1.0im)yz, (-0.0 - 1.0im) + (-0.0 - 1.0im)xz)

julia> H
(y, 2.0 + z - y², 2.0yz)
```
"""
function solve_maxwell(J::NTuple{3,Polynomial{3,T}}; ϵ=1, μ=1, ω=1) where {T}
    ρ = -im / ω * divergence(J)
    k² = ω^2 * ϵ * μ
    A = -μ .* map(j -> solve_helmholtz(j, k²), J)
    φ = -1 / ϵ * solve_helmholtz(ρ, k²)
    E = im * ω .* A .- gradient(φ)
    H = 1 / μ .* curl(A)
    return drop_zeros!.(E), drop_zeros!.(H)
end

function assemble_fastevaluator(args...; kwargs...)
    return error("assemble_fastevaluator not found. Did you forget to import FixedPolynomials ?")
end

function fast_evaluate_with_jacobian!(args...; kwargs...)
    return error("fast_evaluate_with_jacobian! not found. Did you forget to import FixedPolynomials ?")
end

function fast_evaluate_with_gradient!(args...; kwargs...)
    return error("fast_evaluate_with_gradient! not found. Did you forget to import FixedPolynomials ?")
end

function fast_evaluate!(args...; kwargs...)
    return error("fast_evaluate! not found. Did you forget to import FixedPolynomials ?")
end

export
       Polynomial,
       convert_coefs,
       assemble_fastevaluator,
       fast_evaluate_with_jacobian!,
       fast_evaluate_with_gradient!,
       fast_evaluate!,
       solve_helmholtz,
       solve_laplace,
       solve_anisotropic_laplace,
       solve_anisotropic_advect,
       solve_anisotropic_advect_diffuse,
       solve_bilaplace,
       solve_stokes,
       solve_brinkman,
       solve_elastostatic,
       solve_elastodynamics,
       solve_maxwell

end # module (Polynomials)
