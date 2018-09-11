"""
    curvepoint(n, p, U, Pw, u)

Compute point on rational B-Spline curve defined as:

```math
\\mathbf{C}^w(u) = \\sum_{i=0}^n N_{i, p}(u) \\mathbf{P}_i^w
```

where \$ \\mathbf{P}_i^w \$ are the set of weighted control points and weights such that \$ \\mathbf{P}_i^w = (w_ix_i, w_iy_i, w_iz_i, w_i) \$.

(see NURBS, eqn 4.5 and A4.1)

Inputs:
- n : the number of control points minus 1 (the index of the last control point)
- p : the curve order
- U : the knot vector
- Pw : the set of weighted control points and weights
- u : the parametric point of interest

TODO: if u value outside of U vector range is given, function hangs, but doesn't throw error. Need to add a check/error.
"""
function curvepoint(n, p, U, Pw, u)
    span = getspanindex(n, p, u, U)
    # println("span = ", span)
    N = basisfunctions(span+1, u, p, U)
    Cw = 0.0
    for j = 0:p
        # println("N[$j] = ", N[j+1])
        # println("Pw[$(span-p+j)] = ", Pw[span-p+j, :])
        Cw += N[j+1]*Pw[span-p+j+1, :]
    end
    # println(Cw)
    return Cw[1:end-1]/Cw[end]
end

"""
    rationalcurvederivatives(Aders, wders, d)

Compute the point \$ \\mathbf{C}(u) \$ and the derivatives \$ \\mathbf{C}^{(k)}(u) \$  for \$ 1 \\leq k \\leq d \$ where:

```math
\\mathbf{C}^{(k)}(u) = \\frac{ \\mathbf{A}^{(k)}(u) - \\sum_{i=1}^k \\binom{k}{i} w^{(i)}(u) \\mathbf{C}^{(k-1)}(u) }{w(u)}
```

where \$ \\mathbf{A}^{(k)}(u) \$ and \$ w^{(i)}(u) \$ are precomputed using preweighted control points for some parametric point, \$ 0 \\leq u \\leq 1 \$, from ```curvederivatives1``` and are inputs Aders and wders, respectively.

(see NURBS eqn 4.8 and A4.2)
"""
function rationalcurvederivatives(Aders, wders, d)
    CK = zeros(d+1, length(Aders[1, :]))
    for k=0:d
        v = Aders[k+1, :]
        # println("Aders[$k] = ", Aders[k+1, :])
        for i=1:k
            if k >= i
                # println("binom = ",  binomialcoeff(k, i))
                # println("wders[$i] = ", wders[i+1])
                # println("Ck[$(k-i)] = ", CK[k-i+1, :])
                v -= binomialcoeff(k, i)*wders[i+1]*CK[k-i+1, :]
            end
        end
        CK[k+1, :] = v/wders[0+1]
    end
    return CK
end


"""
    curveknotinsertion(np, p, UP, Pw, u, k, s, r)

Compute a new curve from knot insertion. Using the formula:

```math
\\mathbf{Q}_{i, r}^w = \\alpha_{i, r} \\mathbf{Q}_{i, r-1}^w + (1-\\alpha_{i, r}) \\mathbf{Q}_{i-1, r-1}^w
```

where
```math
\\alpha_{i, r} =
\\begin{cases}
      1 & i \\leq k-p+r-1 \\\\
      \\frac{\\bar{u} - u_i}{u_{i+p-r+1} - \\bar{u}_i} & k-p+r \\leq i\\leq k-s \\\\
      0 & i \\geq k-s+1
\\end{cases}
```

(see NURBS eqn 5.15 and A5.1)

Inputs:
- np : the number of control points minus 1 (the index of the last control point) before insertion
- p : the curve order
- UP : the knot vector before insertion
- Pw : the set of weighted control points and weights before insertion
- u : the parametric point of interest
- k : the span index at which the knot is to be inserted.
- s : numer of instances of the new knot alrady present in the knot vector, UP
- r : number of times the new knot is inserted (it is assumed that \$ r+s \\leq p \$)

Outputs:
- nq : the number of control points minus 1 (the index of the last control point) after insertion
- UQ : the knot vecotr after insertion
- Qw : the set of weighted control points and weights after insertion
"""
function curveknotinsertion(np, p, UP, Pw, u, k, s, r)
    mp = np+p+1
    nq = np+r

    #initialize output vectors
    UQ = zeros(length(UP)+r)
    # println("UP, $UP")
    # println("UQ, $UQ")
    Qw = zeros(nq+1, length(Pw[1, :]))
    Rw = zeros(p+1, length(Pw[1, :]))
    #Load new knot vector
    for i=0:k
        UQ[i+1] = UP[i+1]
        # println("UQ$i = ", UQ[i+1])
    end
    # println("UQ, $UQ")
    for i=1:r
        UQ[k+i+1] = u
        # println("UQ$(k+i) = ", UQ[k+i+1])
    end
    # println("UQ, $UQ")
    for i=k+1:mp
        UQ[i+1+r] = UP[i+1]
        # println("UQ$(i+r) = ", UQ[i+1+r])
    end
    # println("UQ, $UQ")
    #Save unaltered control points
    for i=0:k-p
        Qw[i+1, :] = Pw[i+1, :]
    end
    for i=k-s:np
        Qw[i+1+r, :] = Pw[i+1, :]
    end
    for i=0:p-s
        Rw[i+1, :] = Pw[k-p+i+1, :]
    end
    #insert new knot r times
    L = 0.0 #initialize L in this scope
    for j=1:r
        L = k-p+j
        for i=0:p-j-s
            alpha = (u-UP[L+i+1])/(UP[i+1+k+1]-UP[L+i+1])
            Rw[i+1, :] = alpha*Rw[i+1+1, :] + (1.0-alpha)*Rw[i+1, :]
        end
        Qw[L+1, :] = Rw[0+1, :]
        # println("Qw $L = ", Qw[L+1, :])
        Qw[k+r-j-s+1, :] = Rw[p-j-s+1, :]
        # println("Qw $(k+r-j-s) = ", Qw[k+r-j-s+1, :])
    end
    #Load remaining control points
    for i=L+1:k-s
        Qw[i+1, :] = Rw[i+1-L, :]
    end


    return nq, UQ, Qw
end