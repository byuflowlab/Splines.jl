@testset "B-Spline: Find Span Index Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    i = 5

    span = Splines.BSpline.getSpanIndex(p,u,U)
    @test isapprox(span,i,atol=1e15)

    u = 5
    i = 9
    span = Splines.BSpline.getSpanIndex(p,u,U)
    @test isapprox(span,i,atol=1e15)

    u = 0
    i = 1
    span = Splines.BSpline.getSpanIndex(p,u,U)
    @test isapprox(span,i,atol=1e15)

end #Find Span Tests

@testset "B-Spline: Basis Function Tests" begin
    U = [0,0,0,1,2,3,4,4,5,5,5]
    u = 5/2
    p = 2
    i = 5

    bases1 = Splines.BSpline.basisFunctions(i,u,p,U)
    @test isapprox([1/8,6/8,1/8],bases1,atol=1e15)
    # println("bases: ", bases1)
    derivatives = Splines.BSpline.basisFunctionsDerivatives(i,u,p,U)
    # println("derivatives: ")
    # display(derivatives)
    # println()
    # println("actual:")
    # display([1/8 6/8 1/8; -1/2 0 1/2; 1 -2 1])
    # println()
    @test isapprox([1/8, 6/8, 1/8], derivatives[1,:])
    @test isapprox([-1/2, 1.0], derivatives[2:3,1])
    @test isapprox([0.0, -2.0], derivatives[2:3,2])
    @test isapprox([1/2, 1.0], derivatives[2:3,3])
end #Basis Functions Tests