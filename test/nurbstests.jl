import LinearAlgebra

@testset "NURBS: Curve Point Test" begin
    U = [0.0,0,0,1,2,3,3,3]
    w = [1.0,4,1,1,1]
    P = [[0.0, 0], [1, 1], [3, 2], [4, 1], [5, -1]]
    deg = 2

    nurbs = NURBS(deg, U, w, P)

    u = 1

    C = curvepoint(nurbs, u)
    @test C == [7/5, 6/5]
end

@testset "NURBS: Rational Curve Derivative Tests" begin
#Example 4.2 in NURBS Book
    U = [0.0, 0, 0, 1, 1, 1]
    p = 2
    P = [[1.0, 0], [1, 1], [0, 1]]
    w = [1.0, 1, 2]
    nurbs = NURBS(p, U, w, P)

    u = 0
    d = 2
    C = curvepoint(nurbs, u)
    CK = curvederivatives(nurbs, u, d)

    @test C == CK[1]  # function value
    @test CK[2] == [0.0, 2.0]  # first derivative
    @test CK[3] == [-4.0, 0.0]  # second derivative

    u = 1
    C = curvepoint(nurbs, u)
    CK = curvederivatives(nurbs, u, d)
    @test C == CK[1]
    @test CK[2] == [-1, 0.0]
    @test CK[3] == [1, -1]
end

@testset "NURBS: Knot Insertion - Unique Knot" begin
    UP = [0,0,0,0,1,2,3,4,5,5,5,5]
    u = 5/2
    p = 3
    P = [[0,0],[1,1],[2,0],[3,0],[4,1],[3,2],[2,2],[1,1]]
    w = [1;1;1;1;1;1;1;1]
    Pw = [[0,0,1],[1,1,1],[2,0,1],[3,0,1],[4,1,1],[3,2,1],[2,2,1],[1,1,1]]
    np = length(P[:,1])-1
    k = 5
    s = 0
    r = 1

    Qwbyhand = [zeros(length(Pw[1, :][1])) for _ in 1:np+r+1]
    Qwbyhand[1:3,:] = Pw[1:3,:]
    Qwbyhand[4,:] = 5/6*Pw[4,:] + 1/6*Pw[3,:]
    Qwbyhand[5,:] = 1/2*Pw[5,:] + 1/2*Pw[4,:]
    Qwbyhand[6,:] = 1/6*Pw[6,:] + 5/6*Pw[5,:]
    Qwbyhand[7:end,:] = Pw[6:end,:]

    nurbs = Splines.NURBS(p,UP,w,P)
    nq, UQ, Qw = Splines.curveknotinsertion(nurbs, u, r,separated_outputs=true)

    @test nq == np+r
    @test UQ == [0,0,0,0,1,2,5/2,3,4,5,5,5,5]
    @test isapprox(Qw, Qwbyhand, atol=1e-15)

    nurbs2 = Splines.curveknotinsertion(nurbs, u, r)
    nurbs3 = Splines.NURBS(p,UQ,getindex.(Qw,3),getindex.(Qw,[1:2]))

    @test nurbs2.degree == nurbs3.degree
    @test nurbs2.knots == nurbs3.knots
    @test nurbs2.weights == nurbs3.weights
    @test nurbs2.ctrlpts == nurbs3.ctrlpts
end

@testset "NURBS: Knot Insertion - Repeated Knot" begin
    UP = [0,0,0,0,1,2,3,4,5,5,5,5]
    u = 2
    p = 3
    P = [[0,0],[1,1],[2,0],[3,0],[4,1],[3,2],[2,2],[1,1]]
    w = [1;1;1;1;1;1;1;1]
    Pw = [[0.,0,1],[1,1,1],[2,0,1],[3,0,1],[4,1,1],[3,2,1],[2,2,1],[1,1,1]]
    np = length(P[:,1])-1
    k = 5
    s = 1
    r = 1

    Qwbyhand = [zeros(length(Pw[1, :][1])) for _ in 1:np+r+1]
    Qwbyhand[1:3] = Pw[1:3]
    Qwbyhand[4] = 2/3*Pw[4] + 1/3*Pw[3]
    Qwbyhand[5] = 1/3*Pw[5] + 2/3*Pw[4]
    Qwbyhand[6] = Pw[5]
    Qwbyhand[7:end] = Pw[6:end]

    nurbs = Splines.NURBS(p,UP,w,P)
    nq, UQ, Qw = Splines.curveknotinsertion(nurbs, u, r,separated_outputs=true)

    @test nq == np+r
    @test UQ == [0,0,0,0,1,2,2,3,4,5,5,5,5]
    @test isapprox(Qw, Qwbyhand, atol=1e-15)
end

@testset "NURBS: Knot Insertion - Multiple Knots" begin
    U = [0,0,0,0,1,2,3,4,5,5,5,5]
    X = [1.5, 2.5]
    p = 3
    P = [[0,0],[1,1],[2,0],[3,0],[4,1],[3,2],[2,2],[1,1]]
    w = [1;1;1;1;1;1;1;1]
    Pw = [[0,0,1],[1,1,1],[2,0,1],[3,0,1],[4,1,1],[3,2,1],[2,2,1],[1,1,1]]
    n = length(P[:,1])-1
    r = length(X)-1

    nurbs = Splines.NURBS(p,U,w,P)
    Ubar, Qwcalcd = Splines.refineknotvectorcurve(nurbs, X,separated_outputs=true)

    @test Ubar == [0.0,0.0,0.0,0.0,1.0,1.5,2.0,2.5,3.0,4.0,5.0,5.0,5.0,5.0]

    nurbs = Splines.NURBS(p,U,w,P)
    nq, UQ, Qw = Splines.curveknotinsertion(nurbs, X[1], r,separated_outputs=true)
    nurbs = Splines.NURBS(p,UQ,getindex.(Qw,3),getindex.(Qw,[1:2]))
    _, UQ, Qw = Splines.curveknotinsertion(nurbs, X[2], r,separated_outputs=true)

    @test isapprox(Qwcalcd, Qw, atol=1e-15)
end

@testset "NURBS: Unweighted 1 Degree Elevation" begin
    U = [0,0,0,0,3/10,7/10,1,1,1,1]
    p = 3
    P = [[-1,0],[-1.5,1],[-0.5,2],[0.5,2],[1.5,1],[1,0]]
    w = [1.;1;1;1;1;1]
    Pw = [[-1,0,1],[-1.5,1,1],[-0.5,2,1],[0.5,2,1],[1.5,1,1],[1,0,1]]
    n = length(P[:,1])-1
    t = 1

    nurbs = Splines.NURBS(p,U,w,P)
    nh, Uh, Qw = Splines.degreeelevatecurve(nurbs,t,separated_outputs=true)

    @test Uh == [0,0,0,0,0,3/10,3/10,7/10,7/10,1,1,1,1,1]


    s = length(unique(U)) - 2 #number of unique internal knots
    @test nh == n + t*(s+1)


    curvepoints = collect(0:0.01:1)
    Cw1 = zeros(length(curvepoints), length(P[1, :][1]))
    for i = 1:length(curvepoints)
        Cw1[i, :] = Splines.curvepoint(nurbs, curvepoints[i])
    end

    Cw2 = zeros(length(curvepoints), length(P[1, :][1]))
    for i = 1:length(curvepoints)
        nurbs2 = Splines.NURBS(p+t, Uh, getindex.(Qw,3),getindex.(Qw,[1:2]))
        Cw2[i, :] = Splines.curvepoint(nurbs2, curvepoints[i])
    end

    Cw3 = zeros(length(curvepoints), length(P[1, :][1]))
    for i = 1:length(curvepoints)
        nurbs3 = Splines.degreeelevatecurve(nurbs,t)
        Cw3[i, :] = Splines.curvepoint(nurbs3, curvepoints[i])
    end

    @test isapprox(LinearAlgebra.norm(Cw1-Cw2),0.0,atol=1e-14)
    @test isapprox(LinearAlgebra.norm(Cw1-Cw3),0.0,atol=1e-14)
end


@testset "NURBS: Unweighted 2 Degree Elevation" begin
    U = [0,0,0,0,3/10,7/10,1,1,1,1]
    p = 3
    P = [[-1,0],[-1.5,1],[-0.5,2],[0.5,2],[1.5,1],[1,0]]
    w = [1.;1;1;1;1;1]
    Pw = [[-1,0,1],[-1.5,1,1],[-0.5,2,1],[0.5,2,1],[1.5,1,1],[1,0,1]]
    n = length(P[:,1])-1
    t = 2

    nurbs = Splines.NURBS(p,U,w,P)
    nh, Uh, Qw = Splines.degreeelevatecurve(nurbs,t,separated_outputs=true)

    @test Uh == [0,0,0,0,0,0,3/10,3/10,3/10,7/10,7/10,7/10,1,1,1,1,1,1]


    s = length(unique(U)) - 2 #number of unique internal knots
    @test nh == n + t*(s+1)


    curvepoints = collect(0:0.01:1)
    Cw1 = zeros(length(curvepoints), length(P[1, :][1]))
    for i = 1:length(curvepoints)
        Cw1[i, :] = Splines.curvepoint(nurbs, curvepoints[i])
    end

    Cw2 = zeros(length(curvepoints), length(P[1, :][1]))
    for i = 1:length(curvepoints)
        nurbs2 = Splines.NURBS(p+t, Uh, getindex.(Qw,3),getindex.(Qw,[1:2]))
        Cw2[i, :] = Splines.curvepoint(nurbs2, curvepoints[i])
    end

    Cw3 = zeros(length(curvepoints), length(P[1, :][1]))
    for i = 1:length(curvepoints)
        nurbs3 = Splines.degreeelevatecurve(nurbs,t)
        Cw3[i, :] = Splines.curvepoint(nurbs3, curvepoints[i])
    end

    @test isapprox(LinearAlgebra.norm(Cw1-Cw2),0.0,atol=1e-14)
    @test isapprox(LinearAlgebra.norm(Cw1-Cw3),0.0,atol=1e-14)
end


@testset "NURBS: Basis Function Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    P = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
    w = [1;1;1;1;1;1;1]
    u = 5/2
    p = 2
    d = 1
    nurbs = Splines.NURBS(p,U,w,P)
    R, dR = Splines.nurbsbasis(nurbs,u,d)
    @test [1/8,6/8,1/8] == R
    @test [-1/2, 0, 1/2] == dR
end #Basis Functions Tests