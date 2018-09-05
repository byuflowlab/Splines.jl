module NURBS

"""
"""
function curvePoint(n, p, U, Pw, u)
    span = Bspline.getSpanIndex(n,p,u,U)
    N = Bspline.basisFunctions(span,u,p,U)
    return N
end




end #module NURBS