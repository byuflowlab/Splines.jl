@testset "B-Spline: Find Span Index Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    n = length(U)-p-1
    i = 4

    span = Splines.getspanindex(n,p,u,U)
    @test span == i

    u = 5
    i = 8
    span = Splines.getspanindex(n,p,u,U)
    @test span == i

    u = 0
    i = 2
    span = Splines.getspanindex(n,p,u,U)
    @test span == i

end #Find Span Tests

@testset "B-Spline: Basis Function Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    i = 5
    n = p
    bases1 = Splines.basisfunctions(i,u,p,U)
    # println("bases: ", bases1)
    @test [1/8,6/8,1/8] == bases1
    derivatives = Splines.basisfunctionsderivatives(i,u,p,n,U)
    # println("derivatives: ")
    # display(derivatives)
    # println()
    # println("actual:")
    # display([1/8 6/8 1/8; -1/2 0 1/2; 1 -2 1])
    # println()
    @test [1/8, 6/8, 1/8] == derivatives[1,:]
    @test [-1/2, 1.0] == derivatives[2:3,1]
    @test [0.0, -2.0] == derivatives[2:3,2]
    @test [1/2, 1.0] == derivatives[2:3,3]

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
    @test -1/2*P[3,:] + 1/2*P[5,:] == curveDerivatives[2,:]


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

@testset "B-Spline: Interpolation Tests" begin
    Q = [0 0; 3 4; -1 4; -4 0; -4 -3]
    n= 4
    p = 3
    U = [0 0 0 0 28/51 1 1 1 1]
    ubar = [5/17 9/17 14/17]
    N = zeros(3,4)
    for i=1:3
        span = Splines.getspanindex(n,p,ubar[i],U)
        N[i,:] = Splines.basisfunctions(span+1, ubar[i], p, U)
    end
    A = zeros(5,5)
    A[1,1] = 1
    A[5,5] = 1
    for i = 2:3
        for j = 1:4
            A[i,j]=N[i-1,j]
        end
    end
    for i=2:5
        A[4,i]=N[3,i-1]
    end
    P = inv(A)*Q
    m = length(U)-1

    mprime, Uprime, Pprime =Splines.globalcurveinterpolation(n,Q,2,p;knotplacement="chordlength")
    Uprime = reshape(Uprime,1,9)
    @test mprime == m
    @test isapprox(Uprime,U,atol=1e-5)
    @test isapprox(Pprime,P,atol=1e-5)
end
