module NURBS
include("BSpline.jl")
"""
"""
function curvePoint(n, p, U, Pw, u)
    span = BSpline.getSpanIndex(n,p,u,U)
    println("span = ", span)
    N = BSpline.basisFunctions(span,u,p,U)
    Cw = 0.0
    for j = 0:p
        println("N[$j] = ", N[j+1])
        println("Pw[$(span-p+j)] = ", Pw[span-p+j,:])
        Cw += N[j+1]*Pw[span-p+j,:] #don't need a +1 on the j here because span is already +1 from getSpanIndex
    end
    println(Cw)
    return Cw[1:end-1]/Cw[end]
end




end #module NURBS