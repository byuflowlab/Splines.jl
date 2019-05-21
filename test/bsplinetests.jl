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
    n = 4
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
            A[i,j] = N[i-1,j]
        end
    end
    for i=2:5
        A[4,i] = N[3,i-1]
    end
    P = inv(A)*Q
    m = length(U)-1

    mprime, Uprime, Pprime = Splines.globalcurveinterpolation(n,Q,2,p,knotplacement="chordlength")
    Uprime = reshape(Uprime,1,9)
    @test mprime == m
    @test isapprox(Uprime,U,atol=eps())
    @test isapprox(Pprime,P,atol=1e-14)
end



# @testset "B-Spline: Point Projection onto Curve" begin

# #Some more test sets for auxiliary functions
# @testset "Point Projection: Loop Function" begin
# #TODO not sure what to do here, need to think about it more.

# end

# @testset "Point Projection: f(u)" begin
# #just set up some toy function, could be anything, to make sure this works. just do something by hand. (eqn 6.3 in book)
# end

# ## Full Function Tests
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

# #TODO This is code from IGAFoil. Need to put in loop using "normalpoints" that come from projection function outputs.
# ##Code for Obtaining Normals
# #call curve derivatives function
# ders = Splines.curvederivatives1(n, p, knots, weightedcontrolpoints, normalpoint, 1)
# #separate output into parts for NURBS derivatives
# Aders = ders[:, 1:end-1]
# wders = ders[:, end]
# #call NURBS derivatives function
# dR = Splines.rationalcurvederivatives(Aders, wders, 1)
# #just use the first derivative (don't need the zeroeth derivative)
# tangent = dR[end, :]
# #take the tangent: rotate and normalize to obtain unit normal vector
# normal = [tangent[2]; -tangent[1]]/LinearAlgebra.norm(tangent)


# #test dot product of normals and projection vectors.
# @test true


# #points on curve:
# pointsoncurve = Q
# #the projection function should return the values of u that produced the points Q.

# #Test that points R, and Q are the same.
# @test true

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

# @testset "B-Spline: Least Squares Fit"
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

# #Test that the fit is within a specified error, that is, test that at the same U, values, the same Q values are had withing some expected error. One would hope that they would be very close.

# end

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