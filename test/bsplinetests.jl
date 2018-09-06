@testset "B-Spline: Find Span Index Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    n = length(U)-p-1
    i = 5

    span = Splines.getspanindex(n,p,u,U)
    @test span == i

    u = 5
    i = 8
    span = Splines.getspanindex(n,p,u,U)
    @test span == i

    u = 0
    i = 3
    span = Splines.getspanindex(n,p,u,U)
    @test span == i

end #Find Span Tests

@testset "B-Spline: Basis Function Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    i = 5
    n = length(U)-p-1
    bases1 = Splines.basisfunctions(i,u,p,U)
    @test isapprox([1/8,6/8,1/8],bases1,atol=1e15)
    # println("bases: ", bases1)
    derivatives = Splines.basisfunctionsderivatives(i,u,p,n,U)
    # println("derivatives: ")
    # display(derivatives)
    # println()
    # println("actual:")
    # display([1/8 6/8 1/8; -1/2 0 1/2; 1 -2 1])
    # println()
    @test isapprox([1/8, 6/8, 1/8], derivatives[1,:],atol=1e15)
    @test isapprox([-1/2, 1.0], derivatives[2:3,1],atol=1e15)
    @test isapprox([0.0, -2.0], derivatives[2:3,2],atol=1e15)
    @test isapprox([1/2, 1.0], derivatives[2:3,3],atol=1e15)


end #Basis Functions Tests

@testset "B-Spline: Curve Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    n = length(U)-p-1
    P = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]
    d = 1
    curveDerivatives = Splines.curvederivatives1(n, p, U, P, u, d)
    # println("curveDerivatives")
    # display(curveDerivatives)
    # println()
    # println("actual")
    # display(-1/2*P[3,:] + 1/2*P[5,:])
    # println()
    @test isapprox(-1/2*P[3,:] + 1/2*P[5,:], curveDerivatives',atol=1e15)


    U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
    u = 1/2
    i = 4
    d = 1
    p = 3
    P = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]
    n = length(P[:,1])-1
    r1 = 0
    r2 = n
    cprime = Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)

    cp = [15/2*(P[2,:] - P[1,:]), 5*(P[3,:] - P[2,:]), 5*(P[4,:] - P[3,:]), 5*(P[5,:] - P[4,:]), 15/2*(P[6,:] - P[5,:]), 15/2*(P[7,:] - P[6,:])]

    @test cprime[2,:,1] == cp[1,1]
    @test cprime[2,:,2] == cp[2,1]
    @test cprime[2,:,3] == cp[3,1]
    @test cprime[2,:,4] == cp[4,1]
    @test cprime[2,:,5] == cp[5,1]
    @test cprime[2,:,6] == cp[6,1]
end