@testset "B-Spline: Find Span Index Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    i = 4

    span = Splines.getspanindex(p, U, u)
    @test span == i

    u = 5
    i = 7
    span = Splines.getspanindex(p, U, u)
    @test span == i

    u = 0
    i = 2
    span = Splines.getspanindex(p, U, u)
    @test span == i

end #Find Span Tests

@testset "B-Spline: Basis Function Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    span = Splines.getspanindex(p, U, u)
    bases = Splines.basisfunctions(span, p, U, u)
    @test bases[0] == 1/8
    @test bases[1] == 6/8
    @test bases[2] == 1/8
    
    derivatives = Splines.basisfunctionsderivatives(span, p, U, u, n)
    
    @test derivatives[0, 0] == 1/8
    @test derivatives[0, 1] == 6/8
    @test derivatives[0, 2] == 1/8
    @test derivatives[1, 0] == -1/2
    @test derivatives[1, 1] == 0.0
    @test derivatives[1, 2] == 1/2
    @test derivatives[2, 0] == 1.0
    @test derivatives[2, 1] == -2.0
    @test derivatives[2, 2] == 1.0

end #Basis Functions Tests

@testset "B-Spline: Curve Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5.0]
    u = 5/2
    p = 2
    P = [[0.0, 0], [1/2, 1/2], [1.0, 0], [3/2, 1/2], [2.0, 0], [5/2, 1/2], [3.0, 0]]
    d = 1

    bsp = BSpline(p, U, P)

    ders = Splines.curvederivatives(bsp, u, d)
    @test -1/2*P[3] + 1/2*P[5] == ders[2]


    p = 3
    U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
    P = [[0.0, 0], [1/2, 1/2], [1.0, 0], [3/2, 1/2], [2.0, 0], [5/2, 1/2], [3.0, 0]]
    bsp = BSpline(p, U, P)
    
    # u = 1/2
    # i = 4
    # P = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]
    # n = length(P[:,1])-1
    r1 = 1
    r2 = length(P)
    d = 1
    # cprime = Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)

    PK = curvederivativecontrolpoints(bsp, r1, r2, d)

    cp = [15/2*(P[2] - P[1]), 5*(P[3] - P[2]), 5*(P[4] - P[3]), 5*(P[5] - P[4]), 15/2*(P[6] - P[5]), 15/2*(P[7] - P[6])]

    @test PK[2, 1] == cp[1, 1]
    @test PK[2, 2] == cp[2, 1]
    @test PK[2, 3] == cp[3, 1]
    @test PK[2, 4] == cp[4, 1]
    @test PK[2, 5] == cp[5, 1]
    @test PK[2, 6] == cp[6, 1]
end

@testset "B-Spline: Interpolation Tests" begin
    Q = [[0.0, 0], [3.0, 4], [-1.0, 4], [-4.0, 0], [-4.0, -3]]
    p = 3
    U = [0; 0; 0; 0; 28/51; 1; 1; 1; 1]
    ubar = [5/17 9/17 14/17]
    N = zeros(3,4)
    for i = 1:3
        span = Splines.getspanindex(p, U, ubar[i])
        N[i, :] = Splines.basisfunctions(span, p, U, ubar[i])
    end
    A = zeros(5, 5)
    A[1,1] = 1
    A[5,5] = 1
    for i = 2:3
        for j = 1:4
            A[i,j] = N[i-1,j]
        end
    end
    for i=2:5
        A[4,i] = N[3,i-1]
    end
    Qs = hcat(Q...)'
    P = A\Qs

    bsp = globalcurveinterpolation(Q, p)

    @test isapprox(bsp.knots, U, atol=eps())
    @test isapprox(bsp.ctrlpts[1], P[1, :], atol=1e-14)
    @test isapprox(bsp.ctrlpts[2], P[2, :], atol=1e-14)
    @test isapprox(bsp.ctrlpts[3], P[3, :], atol=1e-14)
    @test isapprox(bsp.ctrlpts[4], P[4, :], atol=1e-14)
    @test isapprox(bsp.ctrlpts[5], P[5, :], atol=1e-14)
end



# @testset "B-Spline: Point Projection onto Curve" begin

# #Some more test sets for auxiliary functions
# # @testset "Point Projection: Loop Function" begin
# # #TODO not sure what to do here, need to think about it more.

# # end

# # @testset "Point Projection: f(u)" begin
# # #just set up some toy function, could be anything, to make sure this works. just do something by hand. (eqn 6.3 in book)
# # end

# ## Full Function Tests
# #gbsairfoil spline parameters
# #knot vector for gbs airfoil:
# knots = [0; 0; 0; 0; 0.5; 0.5; 0.5; 1; 1; 1; 1]
# #controlpoints for gbs airfoil:
# controlpoints = [1.0 0.0
#                  1/3 0.147796
#                  0.0 0.057735
#                  0.0 0.0
#                  0.0 -0.057735
#                  1/3 -0.0291073
#                  1.0 0.0]

# weightedcontrolpoints = [controlpoints ones(length(controlpoints[:,1]))]

# p = 3
# n = length(controlpoints[:,1])-1

# #get some points from this spline definition to use for fit tests and stuff.
# u = collect(range(0,stop=1,length=11))
# Q = zeros(length(u),3)
# for i=1:length(u)
#     Q[i,:] = Splines.curvepoint(n, p, knots, weightedcontrolpoints, u[i])
# end
# Q = Q[:,1:2] #just get the x,z coordinates.

# #Projection points
# projectionpoints = [1.2  0.0;
#                     1.2 -0.125;
#                     1.0 -0.125;
#                     0.8 -0.125;
#                     0.6 -0.125;
#                     0.4 -0.125;
#                     0.2 -0.125;
#                     0.0 -0.125;
#                    -0.1 -0.125;
#                    -0.1  0.0;
#                    -0.1  0.2;
#                     0.0  0.2;
#                     0.2  0.2;
#                     0.4  0.2;
#                     0.6  0.2;
#                     0.8  0.2;
#                     1.0  0.2;
#                     1.2  0.2]

# #The projection function will return the points on the spline where these points are projected to. use the function from the IGAFOIL package that gets the curve normals at such points and make sure the normals are aligned with the vectors from the curve to the projection points.
# uproj, R = Splines.projectpoints(n,p,knots,controlpoints,projectionpoints; eps1=eps(), eps2=eps(), ncheckvals=1000)

# R = R[:,1:end-1]

# normal = zeros(length(uproj),2)
# for i=1:length(uproj)
#     if uproj[i] == 0.0 && projectionpoints[i,1] == 1.0
#         uproj[i] += eps()
#     end
#     ##Code for Obtaining Normals
#     #call curve derivatives function
#     ders = Splines.curvederivatives1(n, p, knots, weightedcontrolpoints, uproj[i], 1)
#     #separate output into parts for NURBS derivatives
#     Aders = ders[:, 1:end-1]
#     wders = ders[:, end]
#     #call NURBS derivatives function
#     dR = Splines.rationalcurvederivatives(Aders, wders, 1)
#     #just use the first derivative (don't need the zeroeth derivative)
#     tangent = dR[end, :]
#     # println("tangent = ", tangent)

#     #take the tangent: rotate and normalize to obtain unit normal vector
#     normal[i,:] = [tangent[2]; -tangent[1]]/LinearAlgebra.norm(tangent) #+ R[i,:]
# end


# #test dot product of normals and projection vectors.

# for i=1:length(uproj)
#     if projectionpoints[i,1] > 1.0
#         @test R[i,:] == [1.0; 0.0]
#     else
#         #normalize projection.
#         projnormed = (projectionpoints[i,:]-R[i,:])/LinearAlgebra.norm(projectionpoints[i,:]-R[i,:])
#         @test isapprox(LinearAlgebra.dot(projnormed,normal[i,:]),1.0,atol=eps())
#     end
# end


# #points on curve:
# pointsoncurve = Q
# uproj, R = Splines.projectpoints(n,p,knots,controlpoints,pointsoncurve; eps1=eps(), eps2=eps(), ncheckvals=1000)
# #the projection function should return the values of u that produced the points Q.
# for i=1:length(pointsoncurve[:,1])
#     #Test that points R, and Q are the same.
#     @test isapprox(R[i,1:end-1], Q[i,:], atol=eps())
# end

# end


# @testset "B-Spline: Knot Removal"
#     #Take the knot insertion test and reverse it. Start with the output of the insertion test, and remove the knot that was inserted. you should end up with the same original spline as the input to the insertion test.

#     #TODO This is the knot insertion test code
#     UP = [0,0,0,0,1,2,3,4,5,5,5,5]
#     u = 5/2
#     p = 3
#     P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1 1]
#     w = [1 1 1 1 1 1 1 1]
#     Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1]
#     np = length(P[:,1])-1
#     k = 5
#     s = 0
#     r = 1

#     Qwbyhand = zeros(np+r+1,length(Pw[1,:]))
#     Qwbyhand[1:3,:] = Pw[1:3,:]
#     Qwbyhand[4,:] = 5/6*Pw[4,:] + 1/6*Pw[3,:]
#     Qwbyhand[5,:] = 1/2*Pw[5,:] + 1/2*Pw[4,:]
#     Qwbyhand[6,:] = 1/6*Pw[6,:] + 5/6*Pw[5,:]
#     Qwbyhand[7:end,:] = Pw[6:end,:]

#     nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)

#     @test nq == np+r
#     @test UQ == [0,0,0,0,1,2,5/2,3,4,5,5,5,5]
#     @test isapprox(Qw, Qwbyhand, atol=1e-15)

# end

# @testset "B-spline: Bounded Knot Removal" begin

#     @testset "Bounded Knot Removal: Error Bounds" begin
#     #could probably use the same test as the plain knot removal test.
#     #TODO this is the repeated knot insertion test code
#     UP = [0,0,0,0,1,2,3,4,5,5,5,5]
#     u = 2
#     p = 3
#     P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1 1]
#     w = [1 1 1 1 1 1 1 1]
#     Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1]
#     np = length(P[:,1])-1
#     k = 5
#     s = 1
#     r = 1

#     Qwbyhand = zeros(np+r+1,length(Pw[1,:]))
#     Qwbyhand[1:3,:] = Pw[1:3,:]
#     Qwbyhand[4,:] = 2/3*Pw[4,:] + 1/3*Pw[3,:]
#     Qwbyhand[5,:] = 1/3*Pw[5,:] + 2/3*Pw[4,:]
#     Qwbyhand[6,:] = Pw[5,:]
#     Qwbyhand[7:end,:] = Pw[6:end,:]

#     nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)

#     @test nq == np+r
#     @test UQ == [0,0,0,0,1,2,2,3,4,5,5,5,5]
#     @test isapprox(Qw, Qwbyhand, atol=1e-15)



#     end

#     @testset "Bounded Knot Removal: Knot distribution" begin
#         #example (9.1) in the book for chordlength formulation.
#         Qk = [0 0; 3 4; -1 4; -4 0; -4 -3]
#         m = length(Qk[:,1])-1
#         ubar = [0; 5/17; 9/17; 14/17; 1]
#         ub = Splines.computeubar(m,Qk,"chordlength")

#         @test ubar==ub

#         #TODO do another test with the centripetal formulation, need to do it by hand.

#     end

# end

@testset "B-Spline: Least Squares Fit" begin
#gbsairfoil spline parameters
#knot vector for gbs airfoil:
knots = [0; 0; 0; 0; 0.5; 0.5; 0.5; 1; 1; 1; 1]
#controlpoints for gbs airfoil:
controlpoints = [[1.0, 0.0],
                 [1/3, -0.0291073],
                 [0.0, -0.057735],
                 [0.0, 0.0],
                 [0.0, 0.057735],
                 [1/3, 0.147796],
                 [1.0, 0.0]]

# weightedcontrolpoints = [controlpoints ones(length(controlpoints[:,1]))]
p = 3
bsp = BSpline(p, knots, controlpoints)

# create some points
u = range(0, 1, length=50)
x1 = zeros(50)
y1 = zeros(50)
for i = 1:50
    x1[i], y1[i] = curvepoint(bsp, u[i])
end

pts = [[p[1], p[2]]  for p in zip(x1, y1)]

# fit bspline to these points
n = 30
bspfit = leastsquarescurve(pts, n, p)

# evalute new bspline on a mesh spaced more towards nose
u1 = sin.(range(0, pi/2, length=50))/2.0
u2 = (sin.(range(-pi/2, 0, length=51)) .+ 1)/2.0 .+ 0.5
u = [u1; u2[2:end]]
x2 = zeros(100)
y2 = zeros(100)
for i = 1:100
    x2[i], y2[i] = curvepoint(bspfit, u[i])
end

# linearly interpolate onto same x values so we can compare
using FLOWMath: linear

idx = findfirst(y1 .> 0)
xu1 = x1[idx:end]
yu1 = y1[idx:end]
xl1 = x1[idx-1:-1:1]
yl1 = y1[idx-1:-1:1]

idx = findfirst(y2 .> 0)
xu2 = x2[idx:end]
yu2 = y2[idx:end]
xl2 = x2[idx-1:-1:1]
yl2 = y2[idx-1:-1:1]

yu3 = linear(xu2, yu2, xu1)
yl3 = linear(xl2, yl2, xl1)

@test all(isapprox.(yu3[4:end], yu1[4:end], atol=1e-4))
@test all(isapprox.(yl3[2:end], yl1[2:end], atol=2e-4))
# more error towards nose just because we interpolated in x not b/c of bspline
@test all(isapprox.(yu3[1:3], yu1[1:3], atol=1e-3))
@test all(isapprox.(yl3[1], yl1[1], atol=3e-3))

end

# @testset "B-Spline: Approximating with Error Bound" begin

# #gbsairfoil spline parameters
# #knot vector for gbs airfoil:
# knots = [0; 0; 0; 0; 0.5; 0.5; 0.5; 1; 1; 1; 1]
# #controlpoints for gbs airfoil:
# controlpoints = [1.0 0.0
#                  1/3 -0.0291073
#                  0.0 -0.057735
#                  0.0 0.0
#                  0.0 0.057735
#                  1/3 0.147796
#                  1.0 0.0]

# weightedcontrolpoints = [controlpoints ones(length(controlpoints[:,1]))]

# p = 3
# n = length(controlpoints[:,1])-1

# #get some points from this spline definition to use for fit tests and stuff.
# u = collect(range(0,stop=1,length=50))
# Q = zeros(length(u),3)
# for i=1:length(u)
#     Q[i,:] = Splines.curvepoint(n, p, knots, weightedcontrolpoints, u[i])
# end
# Q = Q[:,1:2] #just get the x,z coordinates.

# #Test that the output spline at the same U values, has the same Q values, within the specified error.

# end