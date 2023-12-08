using ElementaryPDESolutions
using StaticArrays
using Test
using ElementaryPDESolutions: laplacian, divergence, gradient, curl, convert_coefs
using Aqua

# run aqua tests. We disable `unbound_args`` because otherwise signatures like
# `solve_stokes(Q::NTuple{N,PolynomialSolution{N,T}})` fail due to the
# possibility that `Q` is an empty tuple, in which case the element type is not
# defined. Not sure how to fix this, so for now we just disable the test.
Aqua.test_all(ElementaryPDESolutions; unbound_args=false)

@testset "Polynomials" begin
    p1 = Polynomial((0, 0) => 1)
    p2 = Polynomial((0, 0) => 2)
    @test p1 + p2 == Polynomial((0, 0) => 3)
    @test 5 * p1 == Polynomial((0, 0) => 5)
    @test @inferred p1((0.1, 0)) == 1

    p = Polynomial([(1, 0) => 1, (2, 1) => -1]) # x - x²y
    @test @inferred p((0.1, 0)) == 0.1
    @test @inferred p((0.1, 2)) == 0.1 - 0.1^2 * 2
    q = convert_coefs(p, Float64)
    @test @inferred q((0.1, 0)) == 0.1
    @test @inferred q((0.1, 2)) == 0.1 - 0.1^2 * 2
    @test typeof(p) == Polynomial{2,Int64}
    @test typeof(q) == Polynomial{2,Float64}

    @test ElementaryPDESolutions.degree(p1) == 0
    @test ElementaryPDESolutions.degree(Polynomial((2, 0) => 1)) == 2
    @test ElementaryPDESolutions.degree(Polynomial((2, 2, 3) => 1)) == 7

    @test 2 * p1 == Polynomial((0, 0) => 2)
    @test (1 + 0.2 * im) * p1 == Polynomial((0, 0) => 1 + 0.2 * im)

    @test ElementaryPDESolutions.derivative(Polynomial((0, 0) => 1), 1) ==
          Polynomial((0, 0) => 0)
    @test ElementaryPDESolutions.derivative(Polynomial((1, 0) => 1), 1) ==
          Polynomial((0, 0) => 1)
    @test ElementaryPDESolutions.derivative(Polynomial((1, 0) => 1), 2) ==
          Polynomial((0, 0) => 0)
    @test ElementaryPDESolutions.derivative(Polynomial((1, 3) => 2), 2) ==
          Polynomial((1, 2) => 6)

    @test ElementaryPDESolutions.laplacian(Polynomial((2, 0) => 1)) ==
          Polynomial((0, 0) => 2)
    @test ElementaryPDESolutions.laplacian(Polynomial((2, 2) => 1)) ==
          Polynomial((2, 0) => 2) + Polynomial((0, 2) => 2)
    @test ElementaryPDESolutions.laplacian(Polynomial((3, 1) => 3)) ==
          Polynomial((1, 1) => 18)

    for P in (Polynomial((2, 0) => 1), Polynomial((2, 2) => 1), Polynomial((3, 1) => 3))
        @test ElementaryPDESolutions.laplacian(P) == divergence(gradient(P))
    end

    P = (Polynomial((3, 1, 1) => 1), Polynomial((1, 2, 5) => 1), Polynomial((1, 2, 3) => 1))
    curlP = (Polynomial(((1, 1, 3) => 2, (1, 2, 4) => -5)),
             Polynomial(((3, 1, 0) => 1, (0, 2, 3) => -1)),
             Polynomial(((0, 2, 5) => 1, (3, 0, 1) => -1)))
    ElementaryPDESolutions.curl(P) == curlP
    for k in 0:4, j in 0:4, i in 0:4
        P = (Polynomial((i, j, k) => 1), Polynomial((j, i, k) => 1),
             Polynomial((k, j, i) => 1))
        @test iszero(ElementaryPDESolutions.divergence(ElementaryPDESolutions.curl(P)))
        Q = Polynomial((i, j, k) => 1)
        @test all(iszero, ElementaryPDESolutions.curl(ElementaryPDESolutions.gradient(Q)))
        @test ElementaryPDESolutions.curl(ElementaryPDESolutions.curl(P)) ==
              ElementaryPDESolutions.gradient(ElementaryPDESolutions.divergence(P)) .-
              laplacian.(P)
    end
end

@testset "Helmholtz" begin
    for k in 0:3, j in 0:3, i in 0:3
        for κ in (1 // 1, 2 // 1, 3 // 1)
            # 2d
            Q = Polynomial((i, j) => 1 // 1)
            P = @inferred solve_helmholtz(Q; k=κ)
            @test ElementaryPDESolutions.laplacian(P) + κ^2 * P == Q
            # 3d
            Q = Polynomial((i, j, k) => 1 // 1)
            P = solve_helmholtz(Q; k=κ)
            @test ElementaryPDESolutions.laplacian(P) + κ^2 * P == Q
        end
    end
end

@testset "Laplace" begin
    for k in 0:3, j in 0:3, i in 0:3
        # 2d
        Q = Polynomial((i, j) => 1 // 1)
        P = solve_laplace(Q)
        @test ElementaryPDESolutions.laplacian(P) == Q
        # 3d
        Q = Polynomial((i, j, k) => 1 // 1)
        P = solve_laplace(Q)
        @test ElementaryPDESolutions.laplacian(P) == Q
    end
    # test that you cannot pass an inhomogenous polynomials
    Q = Polynomial([(0, 0) => 1, (1, 0) => 2, (0, 1) => 3])
    @test_throws AssertionError solve_laplace(Q)
end

@testset "Anisotropic Laplace" begin
    A = SMatrix{2,2,Float64}(2, 1, 1, 3)
    Q = Polynomial([(1, 3) => 1.0, (2, 2) => 1.0])
    rAQ = ElementaryPDESolutions.multiply_by_anisotropic_r(A, Q, 2)
    @test iszero(ElementaryPDESolutions.drop_zeros!(rAQ -
                                                    Polynomial([(3, 3) => 0.2,
                                                                (4, 2) => 0.6,
                                                                (1, 5) => 0.4]), 10^(-15)))
    @test ElementaryPDESolutions.anisotropic_laplacian(A, Q) ==
          Polynomial([(2, 0) => 6.0, (1, 1) => 26.0, (0, 2) => 10.0])
    @test iszero(ElementaryPDESolutions.drop_zeros!(ElementaryPDESolutions.anisotropic_laplacian(A,
                                                                                                 solve_anisotropic_laplace(A,
                                                                                                                           Q)) -
                                                    Q, 10^(-15)))

    # test that you cannot pass an inhomogenous polynomial
    Q = Polynomial([(0, 0) => 1.0, (1, 0) => 2.0, (0, 1) => 3.0])
    @test_throws AssertionError solve_anisotropic_laplace(A, Q)

    # test that you cannot pass an asymmetric matrix
    A = SMatrix{3,3,Float64}(3, 3, 1.5, 2, 4, 1, 1.5, 1, 3)

    # test that you cannot pass a non-invertible matrix
    A = SMatrix{3,3,Float64}(3, 5, 0.0, 5, 4, 0.0, 0.0, 0.0, 0.0)
    Q = Polynomial([(1, 3, 1) => 1.0, (2, 2, 1) => 1.0])
    @test_throws AssertionError solve_anisotropic_laplace(A, Q)
    A = SMatrix{3,3,Rational{Int64}}(3 // 1, 5 // 1, 0 // 1, 5 // 1, 4 // 1, 0 // 1, 0 // 1,
                                     0 // 1, 0 // 1)
    Q = Polynomial([(1, 3, 1) => 1 // 1, (2, 2, 1) => 1 // 1])
    @test_throws AssertionError solve_anisotropic_laplace(A, Q)

    # Test Rationals
    N = 2
    A = SMatrix{N,N,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, 3 // 1)
    g = Polynomial([(1, 3) => 2 // 1])
    v = solve_anisotropic_laplace(A, g)
    @test ElementaryPDESolutions.anisotropic_laplacian(A, v) == g

    N = 3
    A = SMatrix{N,N,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, 1 // 1, 3 // 1, 1 // 1, 1 // 1,
                                     1 // 1, 3 // 1)
    g = Polynomial([(1, 3, 2) => 2 // 1])
    v = solve_anisotropic_laplace(A, g)
    @test ElementaryPDESolutions.anisotropic_laplacian(A, v) == g

    # Test non-positive anisotropic tensor
    N = 2
    A = SMatrix{N,N,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, -3 // 1)
    g = Polynomial([(1, 3) => 2 // 1])
    v = solve_anisotropic_laplace(A, g)
    @test ElementaryPDESolutions.anisotropic_laplacian(A, v) == g
end

@testset "Anisotropic Advection/Diffusion" begin
    # Advection
    N = 2
    β = SVector{N,Float64}(2, 1)
    g = Polynomial([(1, 3) => 1.0])
    v = solve_anisotropic_advect(β, g)
    @test iszero(ElementaryPDESolutions.drop_zeros!(g -
                                                    sum(β[i] *
                                                        ElementaryPDESolutions.gradient(v)[i]
                                                        for i in 1:N), 10^(-15)))

    N = 3
    β = SVector{N,Float64}(2.0, 1.0, 5.0)
    g = Polynomial([(1, 2, 5) => 2.0])
    v = solve_anisotropic_advect(β, g)
    @test iszero(ElementaryPDESolutions.drop_zeros!(g -
                                                    sum(β[i] *
                                                        ElementaryPDESolutions.gradient(v)[i]
                                                        for i in 1:N), 10^(-14)))

    # Advection + Diffusion
    N = 2
    A = SMatrix{N,N,Float64}(2, 1, 1, 3)
    β = SVector{N,Float64}(2, 1)
    f = Polynomial([(1, 3) => 1.0])
    v = solve_anisotropic_advect_diffuse(A, β, f)
    @test iszero(ElementaryPDESolutions.drop_zeros!(f -
                                                    ElementaryPDESolutions.anisotropic_laplacian(A,
                                                                                                 v)
                                                    -
                                                    sum(β[i] *
                                                        ElementaryPDESolutions.gradient(v)[i]
                                                        for i in 1:N), 10^(-14)))

    N = 3
    A = SMatrix{N,N,Float64}(3, 2, 1.5, 2, 4, 1, 1.5, 1, 3)
    β = SVector{N,Float64}(2, 3, 1)
    f = Polynomial([(1, 3, 1) => 1.0, (2, 2, 1) => 2.0])
    v = solve_anisotropic_advect_diffuse(A, β, f)
    @test iszero(ElementaryPDESolutions.drop_zeros!(f -
                                                    ElementaryPDESolutions.anisotropic_laplacian(A,
                                                                                                 v)
                                                    -
                                                    sum(β[i] *
                                                        ElementaryPDESolutions.gradient(v)[i]
                                                        for i in 1:N), 10^(-13)))

    # Rationals
    N = 2
    A = SMatrix{2,2,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, 3 // 1)
    β = SVector{2,Rational{Int64}}(2 // 1, 1 // 1)
    g = Polynomial([(1, 1) => 2 // 1])
    v = solve_anisotropic_advect_diffuse(A, β, g)
    @test ElementaryPDESolutions.anisotropic_laplacian(A, v) + sum(β[i] *
                                                                   ElementaryPDESolutions.gradient(v)[i] for i in 1:N) ==
          g

    # Test non-positive anisotropic tensor
    A = SMatrix{N,N,Rational{Int64}}(2 // 1, 1 // 1, 1 // 1, -3 // 1)
    v = solve_anisotropic_advect_diffuse(A, β, g)
    @test ElementaryPDESolutions.anisotropic_laplacian(A, v) + sum(β[i] *
                                                                   ElementaryPDESolutions.gradient(v)[i] for i in 1:N) ==
          g
end

@testset "Bilaplace" begin
    for k in 0:4, j in 0:4, i in 0:4
        # 2d
        Q = Polynomial((i, j) => 1 // 1)
        P = solve_bilaplace(Q)
        @test laplacian(laplacian(P)) == Q
        # 3d
        Q = Polynomial((i, j, k) => 1 // 1)
        P = solve_bilaplace(Q)
        @test laplacian(laplacian(P)) == Q
    end
end

@testset "Stokes" begin
    μ = 2 // 1
    # 2d
    I = Iterators.product(0:4, 0:4)
    J = Iterators.product(0:4, 0:4)
    for θi in I, θj in J
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1))
        U, P = solve_stokes(Q; μ)
        @test μ .* laplacian.(U) .- gradient(P) == Q
        @test iszero(divergence(U))
    end
    # 3d
    I = Iterators.product(0:1, 0:2, 0:1)
    J = Iterators.product(0:2, 0:1, 0:1)
    K = Iterators.product(0:2, 0:1, 0:1)
    for θi in I, θj in J, θk in K
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1), Polynomial(θk => 1 // 1))
        U, P = solve_stokes(Q; μ)
        @test μ .* laplacian.(U) .- gradient(P) == Q
        @test iszero(divergence(U))
    end
end

@testset "Brinkman" begin
    Re = 314 // 1
    α = 159 // 1
    Q = (Polynomial([(2, 6) => 5 // 1]), Polynomial([(3, 5) => 8 // 1]))
    U, P = solve_brinkman(Q; Re=Re, α=α)
    @test laplacian.(U) .- α^2 .* U .- Re .* gradient(P) == Q
    @test iszero(divergence(U))

    Q = (Polynomial([(2, 1, 3) => 3 // 1]), Polynomial([(3, 1, 2) => 7 // 1]),
         Polynomial([(2, 3, 1) => 3 // 1]))
    U, P = solve_brinkman(Q; Re=Re, α=α)
    @test laplacian.(U) .- α^2 .* U .- Re .* gradient(P) == Q
    @test iszero(divergence(U))
end

@testset "Elastostatics" begin
    μ = 1 // 2
    ν = 3 // 8
    # 2d
    I = Iterators.product(0:4, 0:4)
    J = Iterators.product(0:4, 0:4)
    for θi in I, θj in J
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1))
        U = solve_elastostatic(Q; μ, ν)
        @test μ / (1 - 2ν) .* gradient(divergence(U)) .+ μ .* laplacian.(U) == Q
    end
    # 3d
    I = Iterators.product(0:1, 0:2, 0:1)
    J = Iterators.product(0:2, 0:1, 0:1)
    K = Iterators.product(0:2, 0:1, 0:1)
    for θi in I, θj in J, θk in K
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1), Polynomial(θk => 1 // 1))
        U = solve_elastostatic(Q; μ, ν)
        @test μ / (1 - 2ν) .* gradient(divergence(U)) .+ μ .* laplacian.(U) == Q
    end
end

@testset "Elastodynamics" begin
    ρ = 1 // 3
    μ = 2 // 1
    ν = 1 // 7
    ω = 5 // 1
    k₂² = ω^2 * ρ / μ
    # 2d
    I = Iterators.product(0:4, 0:4)
    J = Iterators.product(0:4, 0:4)
    for θi in I, θj in J
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1))
        U = solve_elastodynamics(Q; ρ, μ, ν, ω)
        @test μ / (1 - 2ν) .* gradient(divergence(U)) .+ μ .* laplacian.(U) .+
              k₂² * μ .* U == Q
    end
    # 3d
    I = Iterators.product(0:3, 0:2, 0:3)
    J = Iterators.product(0:1, 0:2, 0:1)
    K = Iterators.product(0:2, 0:2, 0:4)
    for θi in I, θj in J, θk in K
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1), Polynomial(θk => 1 // 1))
        U = solve_elastodynamics(Q; ρ, μ, ν, ω)
        @test μ / (1 - 2ν) .* gradient(divergence(U)) .+ μ .* laplacian.(U) .+
              k₂² * μ .* U == Q
    end
end

@testset "Maxwell" begin
    μ = 4 // 1
    ϵ = 3 // 2
    ω = 2 // 5
    I = Iterators.product(0:3, 0:2, 0:3)
    J = Iterators.product(0:1, 0:2, 0:1)
    K = Iterators.product(0:2, 0:2, 0:3)
    for θi in I, θj in J, θk in K
        Q = (Polynomial(θi => 1 // 1), Polynomial(θj => 1 // 1), Polynomial(θk => 1 // 1))
        E, H = solve_maxwell(Q; ϵ=ϵ, μ=μ, ω=ω)

        ρ = -im / ω * divergence(Q)
        poly1 = ϵ * divergence(E) - ρ
        @test iszero(poly1)

        poly2 = im * ω * ϵ .* E .+ curl(H) .- Q
        @test all(iszero, poly2)

        poly3 = -im * ω * μ .* H .+ curl(E)
        @test all(iszero, poly3)

        poly4 = μ * divergence(H)
        @test iszero(poly4)
    end
end
