module PolynomialSolutions

using StaticArrays

"""
    struct Polynomial{N,T}

A polynomial in `N` variables with coefficients of type `T`.
"""
struct Polynomial{N,T}
    order2coeff::Dict{NTuple{N,Int},T}
end

# construct a polynomial from a tuple of pairs
Polynomial(t::NTuple{<:Any,Pair{NTuple{N,Int},T}}) where {N,T} = Polynomial{N,T}(Dict(t))

# construct a polynomial from a vector of pairs
Polynomial(v::Vector{Pair{NTuple{N,Int},T}}) where {N,T} = Polynomial{N,T}(Dict(v))

# construct a polynomial from a single pair
Polynomial(p::Pair{NTuple{N,Int},T}) where {N,T} = Polynomial{N,T}(Dict(p))

# construct a monomial with coefficient 1
monomial(θ::NTuple{N,Int}) where {N} = Polynomial(θ=>Rational(BigInt(1)))
monomial(args...) = monomial(NTuple(args))

function is_homogenous(p::Polynomial{N,T}) where {N,T}
    allequal(sum(θ) for θ in keys(p.order2coeff))
end

function Base.iszero(p::Polynomial)
    q = drop_zeros!(deepcopy(p))
    isempty(q.order2coeff)
end

function drop_zeros!(p::Polynomial)
    for (k,v) in p.order2coeff
        if iszero(v)
            delete!(p.order2coeff, k)
        end
    end
    return p
end

"""
    multiply_by_r(p::Polynomial, k::Int = 2)

Mulitply a polynomial `p` by the monomial `r^k`, where `r = |𝐱|` and `k` is an
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

function Base.:+(p1::Polynomial{N,T}, p2::Polynomial{N,T}) where {N,T}
    acc = deepcopy(p1)
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
Base.:-(p1::Polynomial{N,T}, p2::Polynomial{N,T}) where {N,T} = p1 + (-p2)

function Base.:(==)(p1::Polynomial{N,T}, p2::Polynomial{N,T}) where {N,T}
    q1 = deepcopy(p1) |> drop_zeros!
    q2 = deepcopy(p2) |> drop_zeros!
    return q1.order2coeff == q2.order2coeff
end

# multiply a polynomial by a scalar
function Base.:*(c::T, p::Polynomial{N,T}) where {N,T}
    acc = deepcopy(p)
    for (order, coeff) in acc.order2coeff
        acc.order2coeff[order] = c * coeff
    end
    return acc
end

Base.:*(c,p::Polynomial{N,T}) where {N,T} = T(c)*p

Base.:*(p::Polynomial, c) = c * p

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


function Base.show(io::IO, p::Polynomial{N,T}) where {N,T}
    order2coeff = sort(collect(p.order2coeff))
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
            _monomial_print(io, i, o)
        end
    end
end

# verbose code for pretty printing of monomials using unicode
function _monomial_print(io, i, p)
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
    solve_helmholtz(Q::Polynomial,k=1)

Return the unique polynomial `P` satisfying `ΔP + k²P = Q`.
"""
function solve_helmholtz(Q::Polynomial;k=1)
    n = degree(Q)
    m = floor(Int, n/2)
    P = Q
    ΔⁱQ = laplacian(Q)
    for i in 1:m
        P   = P + (-1/k^2)^i*ΔⁱQ
        ΔⁱQ = laplacian(ΔⁱQ) # next laplacian
    end
    return 1/k^2*P
end

"""
    solve_laplace(Q::Polynomial)

Return a polynomial `P` satisfying `ΔP = Q`.
"""
function solve_laplace(Q::Polynomial{N,T}) where {N,T}
    n = degree(Q)
    γ = (k,p) -> 2*(k+1)*(2k+2p+N) # γₖᵖ
    cₖ  = T(1//γ(0,n)) # c₀
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

Compute a polynomial solution to `Δ²P = Q`.
"""
function solve_bilaplace(Q::Polynomial{N}) where {N}
    P′ = solve_laplace(Q)
    P  = solve_laplace(P′)
    return P
end

"""
    solve_stokes(Q::SVector{N,Polynomial{N,T}};μ=1)

Compute a vector of polymomials `U` and a polynomial `P` satisfying `μΔU - ∇P =
Q` with `∇ ⋅ U = 0`.
"""
function solve_stokes(Q::SVector{N,Polynomial{N,T}};μ=1) where {N,T}
    # u = Δg - ∇ (∇ ⋅ g), p = -μ Δ (∇ ⋅ g), where g solves μΔΔg = Q
    g = map(q->solve_bilaplace(q),Q)
    h = -divergence(g)
    u = μ .* (laplacian.(g) .+ gradient(h))
    p = μ*laplacian(h)
    return u,p
end
solve_stokes(Q::NTuple) = solve_stokes(SVector(Q))

# function solve_elastostatic(Q::NTuple{N,Polynomial{N,T}};ρ=1,μ=1,ν=1) where {N,T}
#     g = map(q->solve_helmholtz(solve_helmholtz(q,k2),k1),Q)
#     u = @. 2*(1-ν)*laplacian(g) + k1^2*g - gradient(divergence(g))
#     return P
# end

export
    Polynomial,
    monomial,
    solve_helmholtz,
    solve_laplace,
    solve_bilaplace,
    solve_stokes

end # module (Polynomials)
