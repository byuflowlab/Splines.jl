"""
    getspanindex(n, p, u, U)

Complete binary search to find span index of vector, U, in which the parametric point, u, lies. (NURBS A2.1)
"""
function getspanindex(n, p, u, U)
    # println("Getting Span Index, u = $u")
    if u == U[n+1+1]
        # println("Special case")
        # println("Returning Span Index")
        return n #special case
    else
        lo = p+1
        hi = n+2
        mid = div((lo+hi), 2)
        while u<U[mid] || u >= U[mid+1]
            if u<U[mid]
                hi = mid
            else
                lo = mid
            end #if low
            mid = div((lo+hi), 2)
        end #while not found
        # println("Returning Span Index")
        return mid-1
    end #if not end
end

"""
    basisFunctions(i, u, p, U)

Calculate the non-vanishing basis functions of the B-Spline of order p, defined by knots U at parametric point, ``u``.

The formula for the basis functions is:

```math
N_{i,0}(u) =
\\begin{cases}
      1 & \\textrm{if } u_i \\leq u \\leq u_{i+1} \\\\
      0 & \\textrm{otherwise}
\\end{cases}
```

```math
N_{i,p}(u) = \\frac{u-u_i}{u_{i+p} - u_i} N_{i,p-1}(u) - \\frac{u_{i+p+1} - u}{u_{i+p+1} - u_{i+1}} N_{i+1,p-1}(u)
```

Note that the algorithm used in ```basisFunctions``` removes redunant calculation and potential division by zero (see NURBS, eqn 2.5 and A2.2).
"""
function basisfunctions(i, u, p, U)
    N = ones(p+1)
    left = zeros(p) #was p+1, not sure if that was needed
    right = zeros(p) #was p+1, not sure if that was needed
    for j=1:p
        left[j] = u-U[i+1-j]
        # println("left = ", left)
        right[j] = U[i+j]-u
        # println("right = ", right)
        saved = 0.0
        for r=0:j-1
            # println("denom = ", right[r+1] + left[j-r])
            temp = N[r+1]/(right[r+1] + left[j-r])
            N[r+1] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        end
        N[j+1] = saved
    end
    return N
end

"""
    basisfunctionsderivatives(i, u, p, n, U)

Calculate the non-vanishing basis functions and derivatives of the B-Spline of order ``p```, defined by knots U at parametric point, ``u``.

The basis function derivative is given by

```math
N_{i,p}^{'} = \\frac{p}{u_{i+p} - u_i} N_{i,p-1}(u) - \\frac{p}{u_{i+p+1} - u_{i+1}} N_{i+1,p-1}(u)
```

(see NURBS, eqn 2.7 and A2.3)

Inputs:

- i : knot span containing u
- u : parametric point of interest
- p : the curve order
- n : the max derivative order (n ≦ p)
- U : the knot vector
"""
function basisfunctionsderivatives(i, u, p, n, U)
    #Initialize
    # n = length(U)-p-1
    ndu = ones(p+1, p+1)
    a = zeros(p+1, p+2)
    ders = zeros(n+1, p+1)
    left = zeros(p+1)
    right = zeros(p+1)

    #---Compute (and save) Basis Functions and Knot Differences
    for j=1:p
    left[j] = u-U[i+1-j]
    right[j] = U[i+j]-u
    saved = 0.0
        for r=0:j-1
            #upper triangle (basis functions)
            ndu[j+1, r+1] = right[r+1] + left[j-r]
            temp = ndu[r+1, j]/ndu[j+1, r+1]
            #lower triangle (knot differences)
            ndu[r+1, j+1] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        end #for r
    ndu[j+1, j+1] = saved
    end #for j

    #Load Basis Functions
    for j=1:p+1
        ders[1, j] = ndu[j, p+1]
    end #for j

    #---Compute Derivatives
    # println("\n\nDerivative Algorithm:\nFor knot span ", i)
    for r=0:p
        # println("For basis fuction: ", r)
        #Set row indices for coefficient array (swaps back and forth as they're computed/used)
        s1 = 0
        s2 = 1
        a[0+1, 0+1] = 1.0
        for k=1:n
            # println("Computing Derivative ", k)
            # println("Coefficient Array Row Indices: s1 = ", s1, "\ts2 = ", s2)
            d = 0.0
            rk = r-k
            pk = p-k
            # println("rk = ", rk)
            # println("pk = ", pk)
            if r >= k
                # println("r >= k")
                a[s2+1, 0+1] = a[s1+1, 0+1]/ndu[pk+1+1, rk+1]
                d = a[s2+1, 0+1]*ndu[rk+1, pk+1]
            end #if r>=k

            if rk >= -1
                # println("rk >= -1, \t j1 = 1")
                j1 = 1
            else
                # println("rk < -1, \t j1 = -rk")
                j1 = -rk
            end #if rk>=-1

            if r-1 <= pk
                # println("r-1 <= pk, \t j2 = k-1")
                j2 = k-1
            else
                # println("r-1 > pk\t j2 = p-r")
                j2 = p-r
            end
            # println("j1 = ", j1, "\tj2 = ", j2)
            #add in condition so ranges work for julia
            if j2 >= j1
                for j=j1:j2
                    # println("j = ", j)
                    # println("ndu[pk+1, rk+j] = ", ndu[pk+1+1, rk+j+1])
                    # println("a[s1, j] = ", a[s1+1, j+1])
                    # println("a[s1, j-1] = ", a[s1+1, j-1+1])
                    a[s2+1, j+1] = (a[s1+1, j+1]-a[s1+1, j-1+1])/ndu[pk+1+1, rk+j+1]
                    # println("a[s2, j] = ", a[s2+1, j+1])
                    # println("ndu[rk+j, pk+1] = ", ndu[rk+j+1, pk+1])
                    d += a[s2+1, j+1]*ndu[rk+j+1, pk+1]
                    # println("d = ", d)
                end #for j
            end

            if r <= pk
                # println("r <= pk")
                # println("a[s1, k-1] = ", a[s1+1, k-1+1])
                # println("ndu[pk+1, r] = ", ndu[pk+1+1, r+1])
                a[s2+1, k+1] = -a[s1+1, k-1+1]/ndu[pk+1+1, r+1]
                # println("a[s2, k-1] = ", a[s2+1, k+1])
                # println("a[s2, k] = ", a[s2+1, k+1])
                # println("ndu[r, pk+1] = ", ndu[r+1, pk+1])
                d += a[s2+1, k+1]*ndu[r+1, pk+1]
                # println("d = ", d)
            end #if r
            ders[k+1, r+1] = d
            # println("derivatives: N", r+2, p, "^", k)
            # display(ders)
            # println()
            #switch rows
            j = s1
            s1 = s2
            s2 = j
            # println("\n")
        end #for k
        # println("\n")
    end #for r
    #Multiply through by the correct factors

    for k=1:n
        for j=0:p
            ders[k+1, j+1] *= factorial(p)/factorial(p-k)
        end
    end

    return ders

end #function

"""
    curvederivatives1(n, p, U, P, u, d)

Compute a curve point and its derivatives up do the dth derivative at parametric point, ``u``. (NURBS, A3.2)

#### Inputs
- n : the number of control points is n+1
- p : the degree of the curve
- U : the knot vector
- P : the control points
- u : the parametric point of interest
- d : derivative order (0 ≤ k ≦ d)
"""
function curvederivatives1(n, p, U, P, u, d)
    du = min(d, p)
    CK = zeros(d+1, length(P[1, :]))
    span = getspanindex(n, p, u, U)
    # println("span: ", span)
    nders = basisfunctionsderivatives(span+1, u, p, du, U)
    # println("ders:")
    # display(nders)
    # println()
    for k=0:du
        for j=0:p
            # println("nders[k, j] = ",  nders[k+1, j+1])
            # println("P[span-p+j] = ", P[span-p+j+1, :])
            CK[k+1, :] += nders[k+1, j+1]*P[span-p+j+1, :]
        end
    end

    return CK

end

"""
    curvederivativecontrolpoints(n, p, U, P, d, r1, r2)

Compute control points of curve derivatives:

```math
\\mathbf{C}^{(k)}(u) = \\sum_{i=0}^{n-k}N_{i,p-k}(u) \\mathbf{P}_i^{(k)}
```
with
```math
\\mathbf{P}_i^{(k)} =
\\begin{cases}
    \\mathbf{P}_i & k=0 \\\\
    \\frac{p-k+1}{u_{i+p+1}-u_{i+k}}\\left(\\mathbf{P}_{i+1}^{(k)} - \\mathbf{P}_i^{(k)} \\right) & k > 0
\\end{cases}
```

(see NURBS, eqn 3.8 and A3.3)

#### Inputs
- n : the number of control points is n+1
- p : the degree of the curve
- U : the knot vector
- P : the control points
- u : the parametric point of interest
- d : derivative order (0 ≤ k ≦ d)
- r1 : first control point index
- r2 : last control point index
"""
function curvederivativecontrolpoints(n, p, U, P, d, r1, r2)
    r = r2-r1
    # println("r2-r1 = $r")
    PK = zeros(d+1, length(P[1, :]), r+1)
    # println("size of PK = ", size(PK))
    for i=0:r
        PK[0+1, :, i+1] = P[r1+i+1, :]
    end
    # println("PK[0, :] = ")
    # display(PK)
    # println()
    if d >= 1
        for k=1:d
            # println("k = $k")
            tmp = p-k+1
            # println("tmp = $tmp")
            for i=0:r-k
                # println("i = $i")
                PK[k+1, :, i+1] = tmp*(PK[k-1+1, :, i+1+1] - PK[k-1+1, :, i+1])/(U[r1+i+p+1+1] - U[r1+i+k+1])
                # println("PKi1 = ", PK[k-1+1, :, i+1+1])
                # println("PKi = ", PK[k-1+1, :, i+1])
                # println("Uip1 = ", U[r1+i+p+1+1])
                # println("Uik = ", U[r1+i+k+1])
                # println("PK[k, i] = ", PK[k+1, :, i+1])
            end
        end
    end

    return PK

end

#There is another curvederivatives algorithm in the book (Algorithm 3.4)

"""
    globalcurveinterpolation(n,Q,r,p; knotplacement)

Interpolate points Q, with a B-Spline of degree p. (NURBS A9.1)

Inputs:
- n : n+1 is number of data points to be interpolated
- Q : coordinates of data points to be interpolated
- r : the number of coordinates per Q (the spacial dimension)
- p : degree of interpalatory spline
- knotplacement : the knot placement scheme; either chordlength (common, uniform parameterization) or centripetal (good for data that takes sharp turns).

Outputs:
- m : number of knots
- U : knot vector
- P : control points
"""
function globalcurveinterpolation(n,Q,r,p; knotplacement="centripetal")
    m = n+p+1

    ##-- Get knot vector
    ubar = zeros(n+1)
    U = zeros(m+1)
    P = zeros(n+1,r)

    # find \bar{u}_k

    ubar[1] = 0
    ubar[end] = 1

    if knotplacement == "centripetal"  #eqn 9.6
        global d = 0.0
        for i=2:n+1
            global d += sqrt(LinearAlgebra.norm(Q[i-1,:] - Q[i,:]))
        end
        for i = 2:n
            ubar[i] = ubar[i-1] + sqrt(LinearAlgebra.norm(Q[i,:]-Q[i-1,:]))/d
        end
    elseif knotplacement == "chordlength" #eqn 9.5
        global d = 0.0
        for i=2:n+1
            global d += LinearAlgebra.norm(Q[i-1,:] - Q[i,:])
        end
        for i = 2:n
            ubar[i] = ubar[i-1] + LinearAlgebra.norm(Q[i,:]-Q[i-1,:])/d
        end
    end

    if knotplacement != "centripetal" && knotplacement != "chordlength" #eqn 9.3
        warn("No valid knot placement scheme selected, using equidistant...")
        ubar = collect(range(0,stop=1,length=n+1))
        U[1:p+1] .= 0
        U[m+1-p:m+1] .= 1
        if 2==n+1-p
            U[p+1+1] = 1/2
        else
            U[p+1+1:m+1-p-1] = range(0,stop=1,length=m+1-p-1-(p+1+1))
        end
    else
        #from \bar{u}_k get the knot vector
        U[1:p+1] .= 0
        U[m+1-p:m+1] .= 1
        for i=2:n+1-p
            U[i+p] = sum(ubar[i:i+p-1])/p
        end
    end

    A = zeros(n+1,n+1)
    for i=1:n+1
        span = Splines.getspanindex(n,p,ubar[i],U)
        N = Splines.basisfunctions(span+1, ubar[i], p, U)
        A[i,span+1-p:span+1] = N
    end

    for i=1:r
        P[:,i] = A\Q[:,i]
    end

    return m, U, P
end


"""
    leastsquarescurve(Q,r,n,p, Wq=[], D=[], s=[], I=[], Wd=[]; knotplacement)

Compute the weighted, constrained, least squares curve fit. (NURBS A9.6)

Inputs:

- Q : Data points to be approximated
- r : number of datapoints.
- Wq : weights of "tightness" of approximation to each data point (values greater than zero indicate unconstrained, values less than zero indicat constraint.)
- D : optional derivatives at any of the points, Q
- s : number of derivatives in D is s+1
- I : Maps the derivatives in D to the corresponding points in Q
- Wd : weights associated with derivatives. Values greater than zero indicate
unconstrained, values less than zero indicate constraint
- n : n+1 control points are used for the fit.
- p : the degree of the curve to fit.
- knotplacement : knot placement scheme ("centripital" or "chordlength")

Outputs:

- U : knot vector
- P : control points
"""
function leastsquarescurve(Q,r,n,p, Wq=[], D=[], s=-1, I=[], Wd=[]; knotplacement="centripetal")

    #initialize output
    m = n+p+1
    U = zeros(m+1)
    P = zeros(n+1,length(Q[1,:]))

    #do some setting up of weights based on inputs.
    if isempty(Wq)
        Wq = ones(length(Q[:,1]))
        ru = length(Q[:,1]) - 1
        rc = -1
    else
        ru = -1
        rc = -1
        for i=1:r+1
            if Wq[i] > 0.0
                ru += 1
            else
                rc += 1
            end
        end
    end

    if isempty(Wd)
        Wd = ones(length(Q[:,1]))
        I = -ones(length(Q[:,1]))
        su = length(Q[:,1])-1
        sc = -1
    else
        su = -1
        sc = -1
        for i=1:s+1
            if Wq[i] > 0.0
                su += 1
            else
                sc += 1
            end
        end
    end

    mu = ru+su+1
    mc = rc+sc+1

    if mc >= n || mc+n >= mu+1
        error("That's not going to work. (see NURBS eqn 9.70)")
    else
        #initialize all the local matrices
        N = zeros(mu+1,n+1)
        M = zeros(mc+1,n+1)
        S = zeros(mu+1)
        T = zeros(mc+1)
        A = zeros(mc+1)
        W = zeros(mu+1,mu+1)
    end

    ##-- set up knots
    ubar = zeros(r+1)
    # find \bar{u}_k
    if knotplacement == "centripetal"  #eqn 9.6
        global d = 0.0
        for i=2:r+1
            global d += sqrt(LinearAlgebra.norm(Q[i-1,:] - Q[i,:]))
        end

        ubar[1] = 0
        ubar[end] = 1

        for i = 2:r
            ubar[i] = ubar[i-1] + sqrt(LinearAlgebra.norm(Q[i,:]-Q[i-1,:]))/d
        end

    elseif knotplacement == "chordlength" #eqn 9.5
        global d = 0.0
        for i=2:r+1
            global d += LinearAlgebra.norm(Q[i-1,:] - Q[i,:])
        end

        ubar[1] = 0
        ubar[end] = 1

        for i = 2:r
            ubar[i] = ubar[i-1] + LinearAlgebra.norm(Q[i,:]-Q[i-1,:])/d
        end

    end

    if knotplacement != "centripetal" && knotplacement != "chordlength" #eqn 9.3
        warn("No valid knot placement scheme selected, using equidistant...")
        U = collect(range(0,stop=1,length=m))
    else
        #from \bar{u}_k get the knot vector
        d = (r+1)/(n-p+1)
        #from \bar{u}_k get the knot vector
        U[1:p+1] .= 0
        U[m+1-p:m+1] .= 1
        for j=2:n+1-p
            i = Int(floor(j*d))
            alpha = j*d-i
            U[j+p] = (1-alpha)*ubar[i-1] + alpha*ubar[i]
        end
     end

    ##-- Set up arrays: N, W, S, T, M

    for k=1:length(Q[1,:]) #do each dimension separately.
        j = 1
        mu2 = 1
        mc2 = 1
        for i=1:r+1
            span = Splines.getspanindex(n,p,ubar[i],U)
            dflag = 0 #derivative flag
            #check for derivative at point
            if j <= s
                if i == I[j]
                    dflag = 1
                end
            end

            if dflag == 0 #if derivative not present
                funs = Splines.basisfunctions(span+1, ubar[i], p, U)
            else
                funs = Splines.basisfunctionsderivatives(span+1, ubar[i], p, 1, U)
            end

            #if point is unconstrained
            if Wq[i] > 0
                W[mu2,mu2] = Wq[i]
                N[mu2,span+1-p:span+1] = funs
                S[mu2] = W[mu2,mu2]*Q[i,k]
                mu2 += 1
            else #if point is constrained
                M[mc2,span+1-p:span+1] = funs
                T[mc2] = Q[i,k]
                mc2 += 1
            end #if unconstrained

            #if derivative given for this point
            if dflag == 1
                if Wd[j] > 0 #unconstrained derivative
                    W[mu2,mu2] = Wd[j]
                    N[mu2,span+1-p:span+1] = funs[2,:]
                    S[mu2] = W[mu2,mu2]*D[j]
                    mu2 += 1
                else #constrained derivative
                    M[mc2,span+1-p:span+1] = funs[2,:]
                    T[mc2] = D[j]
                    mc2 += 1
                end
                j += 1
            end #if dflag
        end #for r

        NtransWN = N'*W*N
        NtransWNinv = LinearAlgebra.inv(NtransWN)
        NtransWS = N'*W*S

        if mc < 0 #if no constraints
            P[:,k] = NtransWN\NtransWS
        else
            A = (M*NtransWNinv*M') \ (M*NtransWNinv*NtransWS-T)
            P[:,k] = NtransWN \ (NtransWS - M'*A)
        end

    end #for dimension

    return U, P

end