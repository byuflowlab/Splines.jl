export binomialcoeff

"""
binomialcoeff(n,i)

    Calculate Binomial Coefficient
"""
function binomialcoeff(n, i)
    return factorial(n)./(factorial(i).*factorial(n-i))
end

"""
    bernsteincoeff(n,i)

Calculate Bernstein Coefficient. (NURBS, eqn 1.8)
"""
function bernsteincoeff(u, n, i)
    return binomialcoeff(n,i) .* u.^i .* (1-u).^(n-i)
end

"""
    CoxdeBoorBezier1D(p, u)

Calculate a point along a Bezier curve at the parametric point, u, based on the control points P using Cox-de Boor Algorithm.
"""
function CoxdeBoorBezier1D(p, u)
    n = length(p[:,1])
    bc = zeros(length(u),n)

    for i=1:n
        bc[:,i] = bernsteincoeff(u, n, i)
    end

    C = bc*p

    return C
end


"""
    deCasteljauBezier1D(p, u)

Calculate a point along a Bezier curve at the parametric point, u, based on the control points P using deCasteljau's Algorithm.
"""
function deCasteljauBezier1D(p, u)
    n = length(p[:,1])
    for k=1:n
        for i=1:n-k
            p[i,:] = (1-u)*p[i,:] + u*p[i+1,:]
        end
    end
    C = p[1,:]

end