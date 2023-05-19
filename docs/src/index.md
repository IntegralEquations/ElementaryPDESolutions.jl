```@meta
CurrentModule = PolynomialSolutions
```

# PolynomialSolutions

*Compute polynomial solutions to certain PDEs given a polynomial source term.*

## Overview 

This package provides functionality for solving

```math
    \mathcal{L}P = Q,
```

where `Q` is a source term of [`Polynomial`](@ref) type (or a `NTuple` of
polynomials for vector-valued problems), `P` is the sought polynomial solution, and
$\mathcal{L}$ is a certain (linear) constant coefficient differential operator.

A typical use case consists of creating a polynomial `Q`

```@repl simple-usecase
using PolynomialSolutions;
Q = Polynomial([(1,2)=>2, (3,0)=>-1])
```

and calling the desired method to obtain `P`:

```@repl simple-usecase
P = solve_helmholtz(Q;k=1)
```

Note that `P` can be used as function:

```@repl simple-usecase
P((0.1,0.2)) # functor interface
```

The following PDEs are currently implemented (see their specific [Docstrings](@ref) for
further details):

- [`solve_laplace`](@ref)
- [`solve_helmholtz`](@ref)
- [`solve_bilaplace`](@ref)
- [`solve_stokes`](@ref)
- [`solve_elastostatic`](@ref)
- [`solve_elastodynamics`](@ref)
- [`solve_maxwell`](@ref)

## Coefficient type, precision, and conversions

The type of coefficients in the polynomial solution `P` is inferred from both
the type of the coefficients of the input polynomial `Q`, and from the parameter
types (e.g. the type of `k` in [`solve_helmholtz`](@ref)). In the presence of
floating point errors, this means that the computed coefficients may be inexact:

```@repl coefficients
using PolynomialSolutions;
Q = Polynomial((2,2)=>1)
P = solve_helmholtz(Q;k=3)
```

The recursive algorithm in `solve_helmholtz` performs repeated divisions by
`k²`; thus, passing a rational types yields an exact result in this case:

```@repl coefficients
Q = Polynomial((3,2)=>big(1//1))
P = solve_helmholtz(Q;k=3//1)
```

You can still convert the coefficients back to e.g. a floating point type:

```@repl coefficients
convert(Polynomial{2,Float64},P)
```

Alternatively, you can use intervals to obtain rigorous bounds on the coefficients:

```@repl coefficients
using IntervalArithmetic
Q = Polynomial((3,2)=>1)
P = solve_helmholtz(Q;k=3..3)
```
