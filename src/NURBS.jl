export curvePoint

"""
    curvepoint(n, p, U, Pw, u)

Compute point on rational B-Spline curve. (NURBS, A4.1)
"""
function curvepoint(n, p, U, Pw, u)
    span = getSpanIndex(n,p,u,U)
    # println("span = ", span)
    N = basisFunctions(span,u,p,U)
    Cw = 0.0
    for j = 0:p
        # println("N[$j] = ", N[j+1])
        # println("Pw[$(span-p+j)] = ", Pw[span-p+j,:])
        Cw += N[j+1]*Pw[span-p+j,:] #don't need a +1 on the j here because span is already +1 from getSpanIndex
    end
    # println(Cw)
    return Cw[1:end-1]/Cw[end]
end
