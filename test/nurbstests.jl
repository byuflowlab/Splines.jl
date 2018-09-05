@testset "NURBS: Curve Point Test" begin
    U = [0,0,0,1,2,3,3,3]
    w = [1,4,1,1,1]
    P = [0 0; 1 1; 3 2; 4 1; 5 -1]
    n = length(P[:,1])-1
    u = 1
    p = 2
    i = 5
    Pw = hcat(w.*P,w)
    println(Pw)
    N = Splines.NURBS.curvePoint(n, p, U, Pw, u)

end