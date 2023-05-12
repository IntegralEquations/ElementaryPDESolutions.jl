using PolynomialSolutions
using StaticArrays
using Test
using PolynomialSolutions: laplacian, divergence, gradient, curl

@testset "Polynomials" begin
    p1 = Polynomial((0,0)=>1)
    p2 = Polynomial((0,0)=>2)
    @test p1+p2 == Polynomial((0,0)=>3)
    @test 5*p1 == Polynomial((0,0)=>5)

    @test PolynomialSolutions.degree(p1) == 0
    @test PolynomialSolutions.degree(Polynomial((2,0)=>1)) == 2
    @test PolynomialSolutions.degree(Polynomial((2,2,3)=>1)) == 7

    @test 2*p1 == Polynomial((0,0)=>2)
    @test (1 + 0.2*im)*p1 == Polynomial((0,0)=>1 + 0.2*im)

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

    @test (PolynomialSolutions.curl(
        SVector((
            Polynomial((3,1,1)=>1),
            Polynomial((1,2,5)=>1),
            Polynomial((1,2,3)=>1))))
        ==
        SVector((
            Polynomial(((1,1,3)=>2,(1,2,4)=>-5)),
            Polynomial(((3,1,0)=>1,(0,2,3)=>-1)),
            Polynomial(((0,2,5)=>1,(3,0,1)=>-1))))
        )
    for k in 0:4, j in 0:4, i in 0:4
        P = SVector((monomial(i, j, k), monomial(j,i,k), monomial(k,j,i)))
        @test iszero(PolynomialSolutions.divergence(PolynomialSolutions.curl(P)))
        Q = monomial(i, j, k)
        @test iszero(PolynomialSolutions.curl(PolynomialSolutions.gradient(Q)))
        @test (PolynomialSolutions.curl(PolynomialSolutions.curl(P)) ==
            PolynomialSolutions.gradient(PolynomialSolutions.divergence(P)) - laplacian.(P))
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
    μ = 2//1
    # 2d
    I = Iterators.product(0:4,0:4)
    J = Iterators.product(0:4,0:4)
    for θi in I, θj in J
        Q = SVector(monomial(θi),monomial(θj))
        U,P = solve_stokes(Q;μ)
        @test μ * laplacian.(U) - gradient(P) == Q
        @test iszero(divergence(U))
    end
    # 3d
    I = Iterators.product(0:1,0:2,0:1)
    J = Iterators.product(0:2,0:1,0:1)
    K = Iterators.product(0:2,0:1,0:1)
    for θi in I, θj in J, θk in K
        Q = SVector(monomial(θi),monomial(θj),monomial(θk))
        U,P = solve_stokes(Q;μ)
        @test μ * laplacian.(U) - gradient(P) == Q
        @test iszero(divergence(U))
    end
end

@testset "Elastostatics" begin
    μ = 1//1
    ν = 3//4
    # 2d
    I = Iterators.product(0:4,0:4)
    J = Iterators.product(0:4,0:4)
    for θi in I, θj in J
        Q = SVector(monomial(θi),monomial(θj))
        U = solve_elastostatic(Q;μ,ν)
        @test μ/(1 - 2ν) * gradient(divergence(U)) + μ * laplacian.(U) == Q
    end
    # 3d
    I = Iterators.product(0:1,0:2,0:1)
    J = Iterators.product(0:2,0:1,0:1)
    K = Iterators.product(0:2,0:1,0:1)
    for θi in I, θj in J, θk in K
        Q = SVector(monomial(θi),monomial(θj),monomial(θk))
        U = solve_elastostatic(Q;μ,ν)
        @test μ/(1 - 2ν) * gradient(divergence(U)) + μ * laplacian.(U) == Q
    end
end

@testset "Elastodynamics" begin
    ρ = 1//3
    μ = 2//1
    ν = 1//7
    ω = 5//1
    k₂² = ω^2*ρ/μ

    # 2d
    I = Iterators.product(0:4,0:4)
    J = Iterators.product(0:4,0:4)
    for θi in I, θj in J
        Q = SVector(monomial(θi),monomial(θj))
        U = solve_elastodynamics(Q;ρ,μ,ν,ω)
        @test -μ/(1-2ν)*gradient(divergence(U)) - μ*laplacian.(U) - k₂²*μ*U == Q
    end
    # 3d
    I = Iterators.product(0:3,0:2,0:3)
    J = Iterators.product(0:1,0:2,0:1)
    K = Iterators.product(0:2,0:2,0:4)
    for θi in I, θj in J, θk in K
        Q = SVector(monomial(θi),monomial(θj),monomial(θk))
        U = solve_elastodynamics(Q;ρ,μ,ν,ω)
        @test -μ/(1-2ν)*gradient(divergence(U)) - μ*laplacian.(U) - k₂²*μ*U == Q
    end
end

@testset "Maxwell" begin
    μ = 4//1
    ϵ = 3//2
    ω = 2//5
    I = Iterators.product(0:3, 0:2, 0:3)
    J = Iterators.product(0:1, 0:2, 0:1)
    K = Iterators.product(0:2, 0:2, 0:3)
    for θi in I, θj in J, θk in K
        Q = (1.0 + 0*im)*SVector(monomial(θi), monomial(θj), monomial(θk))
        # Enforce charge conservation
        ρ = -im/ω*divergence(Q)
        E, H, A, φ = solve_maxwell(Q, ρ; ϵ=ϵ,μ=μ,ω=ω)

        poly1 = ϵ*divergence(E) - ρ
        order2coeff = collect(poly1.order2coeff)
        for (order, coeff) in order2coeff
            @test abs(coeff) < 10^(-14)
        end

        poly2 = im*ω*ϵ*E + curl(H) - Q
        for poly in poly2
            order2coeff = collect(poly.order2coeff)
            for (order, coeff) in order2coeff
                @test abs(coeff) < 10^(-14)
            end
        end

        poly3 = -im*ω*μ*H + curl(E)
        for poly in poly3
            order2coeff = collect(poly.order2coeff)
            for (order, coeff) in order2coeff
                @test abs(coeff) < 10^(-14)
            end
        end

        poly4 = μ*divergence(H)
        order2coeff = collect(poly4.order2coeff)
        for (order, coeff) in order2coeff
            @test abs(coeff) < 10^(-14)
        end

        # Test the gauge condition
        poly5 = divergence(A) - im*ω*ϵ*μ*φ
        polydata = collect(poly5.order2coeff)
        for (order, coeff) in polydata
            @test abs(coeff) < 10^(-14)
        end
    end
end
