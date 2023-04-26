```@meta
CurrentModule = PolynomialSolutions
```

# PolynomialSolutions

Generate polynomial solutions to certain PDEs given a polynomial source term.

## Helmholtz

The function [`solve_helmholtz`](@ref) takes a [`Polynomial`](@ref) `Q` and
returns the unique [`Polynomial`](@ref) `P` satisfying

```math
    \Delta P + P = Q
```

Some solutions in 2D:

```@example
using PolynomialSolutions
for j in 0:4, i in 0:4
    Q = monomial(i,j)
    P = solve_helmholtz(Q)
    @show P, Q
end
```

And some solutions in 3D:

```@example
using PolynomialSolutions
for k in 0:4, j in 0:4, i in 0:4
    Q = monomial(i,j,k)
    P = solve_helmholtz(Q)
    @show P, Q
end
```

## Laplace

The function [`solve_laplace`](@ref) takes a [`Polynomial`](@ref) `Q` and
returns a [`Polynomial`](@ref) `P` satisfying

```math
    \Delta P = Q
```

Some solutions in 2D:

```@example
using PolynomialSolutions
for j in 0:4, i in 0:4
    Q = monomial(i,j)
    P = solve_laplace(Q)
    @show P, Q
end
```

And some solutions in 3D:

```@example
using PolynomialSolutions
for k in 0:4, j in 0:4, i in 0:4
    Q = monomial(i,j,k)
    P = solve_laplace(Q)
    @show P, Q
end
```

## Bi-Laplace

The function [`solve_bilaplace`](@ref) takes a [`Polynomial`](@ref) `Q` and
returns a [`Polynomial`](@ref) `P` satisfying

```math
    \Delta^2 P = Q
```

Some solutions in 2D:

```@example
using PolynomialSolutions
for j in 0:4, i in 0:4
    Q = monomial(i,j)
    P = solve_bilaplace(Q)
    @show P, Q
end
```

And some solutions in 3D:

```@example
using PolynomialSolutions
for k in 0:4, j in 0:4, i in 0:4
    Q = monomial(i,j,k)
    P = solve_bilaplace(Q)
    @show P, Q
end
```

## Stokes

The function [`solve_stokes`](@ref) takes an `NTuple` of [`Polynomial`](@ref)s `Q` and
returns two [`Polynomial`](@ref)s, `U` and `P`, satisfying

```math
    \begin{align*}
        \Delta U - \nabla P &= Q \\
        \nabla \cdot U &= 0
    \end{align*}
```

```@example
using PolynomialSolutions
Q   = (monomial(1,0),monomial(0,0))
U,P = solve_stokes(Q)
@show U,P
```

## Docstrings

```@autodocs
Modules = [PolynomialSolutions]
```
