@testset "Bezier Test Set" begin

    checkdata = [0.0 0.0;
    0.00161561 0.00707209;
    0.012681 0.0264594;
    0.0419753 0.054321;
    0.0975461 0.0853528;
    0.186709 0.112788;
    0.316049 0.128395;
    0.491419 0.122481;
    0.717939 0.0838897;
    1.0 0.0]

    P = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition
    u = linspace(0,1,10) #global parameter

    curve1 = Splines.simple_bezier1D(P, u)
    @test isapprox(curve1,checkdata, atol=1e15)


    curve2 = zeros(length(u),2)
    for i = 1:length(u)
        curve2[i,:] = Splines.decasteljau_bezier1D(P, u[i])
    end
    @test isapprox(curve2,checkdata, atol=1e15)

end #Bezier Tests