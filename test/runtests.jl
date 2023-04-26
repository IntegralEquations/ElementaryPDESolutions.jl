using PolynomialSolutions
using StaticArrays
using Test
using PolynomialSolutions: laplacian, divergence, gradient

@testset "Polynomials" begin
    p1 = Polynomial((0,0)=>1)
    p2 = Polynomial((0,0)=>2)
    @test p1+p2 == Polynomial((0,0)=>3)
    @test 5*p1 == Polynomial((0,0)=>5)

    @test PolynomialSolutions.degree(p1) == 0
    @test PolynomialSolutions.degree(Polynomial((2,0)=>1)) == 2
    @test PolynomialSolutions.degree(Polynomial((2,2,3)=>1)) == 7

    @test PolynomialSolutions.derivative(Polynomial((0,0)=>1),1) == Polynomial((0,0)=>0)
    @test PolynomialSolutions.derivative(Polynomial((1,0)=>1),1) == Polynomial((0,0)=>1)
    @test PolynomialSolutions.derivative(Polynomial((1,0)=>1),2) == Polynomial((0,0)=>0)
    @test PolynomialSolutions.derivative(Polynomial((1,3)=>2),2) == Polynomial((1,2)=>6)

    @test PolynomialSolutions.laplacian(Polynomial((2,0)=>1)) == Polynomial((0,0)=>2)
    @test PolynomialSolutions.laplacian(Polynomial((2,2)=>1)) == Polynomial((2,0)=>2) + Polynomial((0,2)=>2)
    @test PolynomialSolutions.laplacian(Polynomial((3,1)=>3)) == Polynomial((1,1)=>18)

    for P in (Polynomial((2,0)=>1),Polynomial((2,2)=>1),Polynomial((3,1)=>3))
        @test PolynomialSolutions.laplacian(P) == divergence(gradient(P))
    end
end

@testset "Helmholtz" begin
    for k in 0:3, j in 0:3, i in 0:3
        for κ in (1//1,2//1,3//1)
            # 2d
            Q = monomial(i,j)
            P = solve_helmholtz(Q;k=κ)
            @test PolynomialSolutions.laplacian(P) + κ^2*P == Q
            # 3d
            Q = monomial(i,j,k)
            P = solve_helmholtz(Q,k=κ)
            @test PolynomialSolutions.laplacian(P) + κ^2*P == Q
        end
    end
end

@testset "Laplace" begin
    for k in 0:3, j in 0:3, i in 0:3
        # 2d
        Q = monomial(i,j)
        P = solve_laplace(Q)
        @test PolynomialSolutions.laplacian(P) == Q
        # 3d
        Q = monomial(i,j,k)
        P = solve_laplace(Q)
        @test PolynomialSolutions.laplacian(P) == Q
    end
end

@testset "Bilaplace" begin
    for k in 0:4, j in 0:4, i in 0:4
        # 2d
        Q = monomial(i,j)
        P = solve_bilaplace(Q)
        @test laplacian(laplacian(P)) == Q
        # 3d
        Q = monomial(i,j,k)
        P = solve_bilaplace(Q)
        @test laplacian(laplacian(P)) == Q
    end
end

@testset "Stokes" begin
    # 2d
    I = Iterators.product(0:4,0:4)
    J = Iterators.product(0:4,0:4)
    for θi in I, θj in J
        Q = SVector(monomial(θi),monomial(θj))
        U,P = solve_stokes(Q)
        @test laplacian.(U) - gradient(P) == Q
        @test iszero(divergence(U))
    end
    # 3d
    I = Iterators.product(0:1,0:2,0:1)
    J = Iterators.product(0:2,0:1,0:1)
    K = Iterators.product(0:2,0:1,0:1)
    for θi in I, θj in J, θk in K
        Q = SVector(monomial(θi),monomial(θj),monomial(θk))
        U,P = solve_stokes(Q)
        @test laplacian.(U) - gradient(P) == Q
        @test iszero(divergence(U))
    end
end
