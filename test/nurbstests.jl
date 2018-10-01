@testset "NURBS: Curve Point Test" begin
    U = [0,0,0,1,2,3,3,3]
    w = [1,4,1,1,1]
    P = [0 0; 1 1; 3 2; 4 1; 5 -1]
    n = length(P[:,1])-1
    u = 1
    p = 2
    i = 5
    Pw = hcat(w.*P,w)
    # println(Pw)
    C = Splines.curvepoint(n, p, U, Pw, u)
    # println(C)
    @test C == [7/5, 6/5, 1]
end

@testset "NURBS: Rational Curve Derivative Tests" begin
    U = [0,0,0,1,1,1]
    u = 0
    p = 2
    P = [1 0; 1 1; 0 1]
    w = [1 1 2]
    Pw = [1 0 1; 1 1 1; 0 2 2]
    n = length(P[:,1])-1
    d = 2
    ders = Splines.curvederivatives1(n, p, U, Pw, u, d)
    Aders = ders[:,1:end-1]
    wders = ders[:,end]
    CK = Splines.rationalcurvederivatives(Aders, wders, d)
    @test CK[end,:] == [-4, 0]

    u = 1
    ders = Splines.curvederivatives1(n, p, U, Pw, u, d)
    Aders = ders[:,1:end-1]
    wders = ders[:,end]
    CK = Splines.rationalcurvederivatives(Aders, wders, d)
    @test CK[end,:] == [1, -1]

end

@testset "NURBS: Knot Insertion - Unique Knot" begin
    UP = [0,0,0,0,1,2,3,4,5,5,5,5]
    u = 5/2
    p = 3
    P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1 1]
    w = [1 1 1 1 1 1 1 1]
    Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1]
    np = length(P[:,1])-1
    k = 5
    s = 0
    r = 1

    Qwbyhand = zeros(np+r+1,length(Pw[1,:]))
    Qwbyhand[1:3,:] = Pw[1:3,:]
    Qwbyhand[4,:] = 5/6*Pw[4,:] + 1/6*Pw[3,:]
    Qwbyhand[5,:] = 1/2*Pw[5,:] + 1/2*Pw[4,:]
    Qwbyhand[6,:] = 1/6*Pw[6,:] + 5/6*Pw[5,:]
    Qwbyhand[7:end,:] = Pw[6:end,:]

    nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)

    # println("By Hand")
    # display(Qwbyhand)
    # println()
    # println("Calculated")
    # display(Qw)
    # println()

    @test nq == np+r
    @test UQ == [0,0,0,0,1,2,5/2,3,4,5,5,5,5]
    @test isapprox(Qw, Qwbyhand, atol=1e-15)
end

@testset "NURBS: Knot Insertion - Repeated Knot" begin
    UP = [0,0,0,0,1,2,3,4,5,5,5,5]
    u = 2
    p = 3
    P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1 1]
    w = [1 1 1 1 1 1 1 1]
    Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1]
    np = length(P[:,1])-1
    k = 5
    s = 1
    r = 1
    # println(np)
    # println(np+p+1)

    Qwbyhand = zeros(np+r+1,length(Pw[1,:]))
    Qwbyhand[1:3,:] = Pw[1:3,:]
    Qwbyhand[4,:] = 2/3*Pw[4,:] + 1/3*Pw[3,:]
    Qwbyhand[5,:] = 1/3*Pw[5,:] + 2/3*Pw[4,:]
    Qwbyhand[6,:] = Pw[5,:]
    Qwbyhand[7:end,:] = Pw[6:end,:]

    nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)

    # println("By Hand")
    # display(Qwbyhand)
    # println()
    # println("Calculated")
    # display(Qw)
    # println()

    @test nq == np+r
    @test UQ == [0,0,0,0,1,2,2,3,4,5,5,5,5]
    @test isapprox(Qw, Qwbyhand, atol=1e-15)
end

@testset "NURBS: Knot Insertion - Multiple Knots" begin
    U = [0,0,0,0,1,2,3,4,5,5,5,5]
    X = [1.5, 2.5]
    p = 3
    P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1 1]
    w = [1 1 1 1 1 1 1 1]
    Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1]
    n = length(P[:,1])-1
    r = length(X)-1
    # println(n)
    # println(n+p+1)

    Ubar, Qwcalcd = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)

    # println("Calculated Ubar")
    # println(Ubar)
    # println()
    # println("Correct Ubar")
    # println([0,0,0,0,1,1.5,2,2.5,3,4,5,5,5,5])
    # println(Pw)

    @test Ubar == [0.0,0.0,0.0,0.0,1.0,1.5,2.0,2.5,3.0,4.0,5.0,5.0,5.0,5.0]

    nq, UQ, Qw = Splines.curveknotinsertion(n, p, U, Pw, X[1], 4, 0, 1)
    _, UQ, Qw = Splines.curveknotinsertion(nq, p, UQ, Qw, X[2], 6, 0, 1)
    # println(UQ)
    # println(Qw)
    # println(Qwcalcd)

    @test isapprox(Qwcalcd, Qw, atol=1e-15)
end

@testset "NURBS: 1 Degree Elevation" begin
    U = [0,0,0,0,3/10,7/10,1,1,1,1]
    p = 3
    P = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]
    w = [1 1 1 1 1 1]
    Pw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]
    n = length(P[:,1])-1
    t = 1
    # println(n)
    # println(n+p+1)

    nh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)

    # println("Calculated Uh")
    # println(Uh)
    # println()
    # println("Correct Uh")
    # println([0,0,0,0,0,3/10,3/10,7/10,7/10,1,1,1,1,1])

    @test Uh == [0,0,0,0,0,3/10,3/10,7/10,7/10,1,1,1,1,1]


    s = length(unique(U)) - 2 #number of unique internal knots
    @test nh == n + t*(s+1)


    curvepoints = collect(0:0.01:1)
    Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
    for i = 1:length(curvepoints)
    Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])
    end

    Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
    for i = 1:length(curvepoints)
    Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])
    end
    @test isapprox(norm(Cw1-Cw2),0.0,atol=1e-10)
end


@testset "NURBS: 2 Degree Elevation" begin
    U = [0,0,0,0,3/10,7/10,1,1,1,1]
    p = 3
    P = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]
    w = [1 1 1 1 1 1]
    Pw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]
    n = length(P[:,1])-1
    t = 2
    # println(n)
    # println(n+p+1)

    nh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)

    # println("Calculated Uh")
    # println(Uh)
    # println()
    # println("Correct Uh")
    # println([0,0,0,0,0,3/10,3/10,7/10,7/10,1,1,1,1,1])

    @test Uh == [0,0,0,0,0,0,3/10,3/10,3/10,7/10,7/10,7/10,1,1,1,1,1,1]


    s = length(unique(U)) - 2 #number of unique internal knots
    @test nh == n + t*(s+1)


    curvepoints = collect(0:0.01:1)
    Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
    for i = 1:length(curvepoints)
    Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])
    end
    Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
    for i = 1:length(curvepoints)
    Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])
    end

    @test isapprox(norm(Cw1-Cw2),0.0,atol=1e-10)
end