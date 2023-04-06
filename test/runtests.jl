using PolynomialSolutions
using Test
using PolynomialSolutions: laplacian

@testset "Polynomials" begin
    p1 = Polynomial((0,0)=>1)
    p2 = Polynomial((0,0)=>2)
    @test p1+p2 == Polynomial((0,0)=>3)
    @test 5*p1 == Polynomial((0,0)=>5)

    @test PolynomialSolutions.degree(p1) == 0
    @test PolynomialSolutions.degree(Polynomial((2,0)=>1)) == 2
    @test PolynomialSolutions.degree(Polynomial((2,2,3)=>1)) == 7

    @test PolynomialSolutions.laplacian(Polynomial((2,0)=>1)) == Polynomial((0,0)=>2)
    @test PolynomialSolutions.laplacian(Polynomial((2,2)=>1)) == Polynomial((2,0)=>2) + Polynomial((0,2)=>2)
    @test PolynomialSolutions.laplacian(Polynomial((3,1)=>3)) == Polynomial((1,1)=>18)

end

@testset "Helmholtz" begin
    for k in 0:3, j in 0:3, i in 0:3
        # 2d
        Q = monomial(i,j)
        P = solve_helmholtz(Q)
        @test PolynomialSolutions.laplacian(P) + P == Q
        # 3d
        Q = monomial(i,j,k)
        P = solve_helmholtz(Q)
        @test PolynomialSolutions.laplacian(P) + P == Q
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
