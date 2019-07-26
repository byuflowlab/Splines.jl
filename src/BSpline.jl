"""
    getspanindex(n, p, u, U) [A2.1]

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
            # println("lo: ", lo)
            # println("hi: ", hi)
            mid = div((lo+hi), 2)
            # println("mid: ", mid)
            if U[mid]==U[mid+1] && abs(lo-hi)==1
                break
            end
        end #while not found
        # println("Returning Span Index")
        return mid-1
    end #if not end
end

"""
    basisfunctions(i, u, p, U)

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
    singlebasisfunction(p,m,U,i,u)

Compute single basis function N_i^p. (NURBS A2.4)

Inputs:
- p : the basis degree up to the the curve order
- m : there are m+1 knots in U
- U : the knot vector
- i : the index of the basis (the i in N_i)
- u : parametric point of interest

Outputs:
- Nip: the N_i^p basis function value
"""
function singlebasisfunction(p,m,U,i,u)
    if ( i==1 && u == U[1] ) || ( i == m-p-1+1 && u == U[m+1] )  #special case
        return 1.0
    end

    if ( u < U[i] ) || ( u >= U[i+p+1] ) #local property
        return 0.0
    end

    N = zeros(p+1)
    for j=1:p+1 #initialize zeroth-degree functions
        if u>=U[i+j] && u < U[i+j+1]
            N[j] = 1.0
        else
            N[j] = 0.0
        end
    end

    for k=1:p #compute triangular table.
        if N[1]==0
            saved = 0.0
        else
            saved = ((u-U[i])*N[1])/(U[i+k]-U[i])
        end

        for j=1:p-k+1
            Uleft = U[i+j+1]
            Uright = U[i+j+k+1]
            if N[j+1]==0.0
                N[j] = saved
                saved = 0.0
            else
                temp = N[j+1]/(Uright-Uleft)
                N[j] = saved+(Uright-u)*temp
                saved = (u-Uleft)*temp
            end #if
        end #for
    end #for

    return N[1]
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
    computeubar(r,Q,knotplacement)

Compute ubar from points, Q, given knotplacement scheme.

Inputs:
- r : there are r+1 datapoints
- Q : array of data points
- knotplacement : type of placement, centripetal or chordlength

Outputs:
- ubar : well distributed parametric points.
"""
function computeubar(r,Q,knotplacement)
    ubar = zeros(r+1)
    ubar[1] = 0
    ubar[end] = 1
    # find \bar{u}_k
    if knotplacement == "centripetal"  #eqn 9.6
        global d = 0.0
        for i=2:r+1
            global d += sqrt(LinearAlgebra.norm(Q[i-1,:] - Q[i,:]))
        end

        for i = 2:r
            ubar[i] = ubar[i-1] + sqrt(LinearAlgebra.norm(Q[i,:]-Q[i-1,:]))/d
        end

    elseif knotplacement == "chordlength" #eqn 9.5
        global d = 0.0
        for i=2:r+1
            global d += LinearAlgebra.norm(Q[i-1,:] - Q[i,:])
        end

        for i = 2:r
            ubar[i] = ubar[i-1] + LinearAlgebra.norm(Q[i,:]-Q[i-1,:])/d
        end

    else
        ubar = collect(range(0,stop=1,length=r+1))
    end

    return ubar

end

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
- m : number of knots-1
- U : knot vector
- P : control points
"""
function globalcurveinterpolation(n,Q,r,p; knotplacement="centripetal")
    m = n+p+1

    ##-- Get knot vector
    ubar = computeubar(length(Q[:,1])-1,Q,knotplacement)
    U = zeros(m+1)
    P = zeros(n+1,r)

    if knotplacement != "centripetal" && knotplacement != "chordlength" #eqn 9.3
        warn("No valid knot placement scheme selected, using equidistant...")
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
        # println("ubar: ", ubar)
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
    surfacepoint(n, p, U, m, q, V, P, u, v) [A3.5]
Compute the z component of a surface from u,v.

Inputs:
    n - number of u control points
    p - degree of u direction curve
    U - u direction knot vector
    m - number of v control points
    q - degree of v direction curve
    V - v knot vector
    P - Control points matrix (n)x(m)? Possibly different.
        Is P a 3D matrix?
    u - u coordinate
    v - v coordinate
Outputs:
    S - Point on surface
"""

function surfacepoint(n, p, U, m, q, V, P, u, v)
    #MUST obey these identities.
    # r = n + p + 1 # U has r + 1 knots
    # s = m + q + 1 # V has s + 1 knots
    r = length(U)-1
    s = length(V)-1
    if r!= n + p + 1
        error("U Knot and control point length not consistent with degree.")
    end
    if s!= m + q + 1
        error("V Knot and control point length not consistent with degree.")
    end

    #TODO: Calculate r, s, n, m based on inputs.

    #QUESTION: Do I need to adjust the span?
    #Find Span and basis of u and v directions
    uspan = Splines.getspanindex(n, p, u, U)
    Nu = Splines.basisfunctions(uspan+1, u, p, U)
    vspan = Splines.getspanindex(m, q, v, V)
    Nv = Splines.basisfunctions(vspan+1, v, q, V)
    uind = uspan-p #TODO: needs currently goes to 0. However,
    # we add k below, so it may not be a problem.
    S = 0.0
    # println("uspan: ", uspan)
    # println("Nu: ")
    # display(Nu)
    # println("vspan: ", vspan)
    # println("Nv: ")
    # display(Nv)
    # println("uind: ", uind)
    #QUESTION: Am I going to have to iterate through dimensions?
    for l=1:q+1
        temp = 0.0
        vind = vspan-q+l #QUESTION: I added a -1, but I don't know if that's the problem...
        # println("vind: ", vind)
        # println("l: ", l)
        for k=1:p+1
            temp = temp + Nu[k]*P[uind+k, vind]
            # println("temp: ", temp)
            # println("k: ", k)
        end
        S = S + Nv[l]*temp
    end
    return S
end

# """
#     leastsquarescurve(Q,r,n,p, Wq=[], D=[], s=[], I=[], Wd=[]; knotplacement)

# Compute the weighted, constrained, least squares curve fit. (NURBS A9.6)

# Inputs:

# - Q : Data points to be approximated
# - r : there are r+1 datapoints.
# - Wq : weights of "tightness" of approximation to each data point (values greater than zero indicate unconstrained, values less than zero indicat constraint.)
# - D : optional derivatives at any of the points, Q
# - s : number of derivatives in D is s+1
# - I : Maps the derivatives in D to the corresponding points in Q
# - Wd : weights associated with derivatives. Values greater than zero indicate
# unconstrained, values less than zero indicate constraint
# - n : n+1 control points are used for the fit.
# - p : the degree of the curve to fit.
# - knotplacement : knot placement scheme ("centripital" or "chordlength")

# Outputs:

# - U : knot vector
# - P : control points
# """
# function leastsquarescurve(Q,r,n,p,ubar=[],U=[], Wq=[], D=[], s=-1, I=[], Wd=[]; knotplacement="centripetal")

#     #initialize output
#     m = n+p+1
#     P = zeros(n+1,length(Q[1,:]))

#     #do some setting up of weights based on inputs.
#     if isempty(Wq)
#         Wq = ones(length(Q[:,1]))
#         ru = length(Q[:,1]) - 1
#         rc = -1
#     else
#         ru = -1
#         rc = -1
#         for i=1:r+1
#             if Wq[i] > 0.0
#                 ru += 1
#             else
#                 rc += 1
#             end
#         end
#     end

#     if isempty(Wd)
#         Wd = ones(length(Q[:,1]))
#         I = -ones(length(Q[:,1]))
#         su = length(Q[:,1])-1
#         sc = -1
#     else
#         su = -1
#         sc = -1
#         for i=1:s+1
#             if Wq[i] > 0.0
#                 su += 1
#             else
#                 sc += 1
#             end
#         end
#     end

#     mu = ru+su+1
#     mc = rc+sc+1

#     if mc >= n || mc+n >= mu+1
#         error("That's not going to work. (see NURBS eqn 9.70)")
#     else
#         #initialize all the local matrices
#         N = zeros(mu+1,n+1)
#         M = zeros(mc+1,n+1)
#         S = zeros(mu+1)
#         T = zeros(mc+1)
#         A = zeros(mc+1)
#         W = zeros(mu+1,mu+1)
#     end

#     ##-- set up knots
#     if isempty(U)
#         U = zeros(m+1)
#         ubar = computeubar(r,Q,knotplacement)

#         if knotplacement != "centripetal" && knotplacement != "chordlength" #eqn 9.3
#             warn("No valid knot placement scheme selected, using equidistant...")
#             U = collect(range(0,stop=1,length=m))
#         else
#             #from \bar{u}_k get the knot vector
#             d = (r+1)/(n-p+1)
#             #from \bar{u}_k get the knot vector
#             U[1:p+1] .= 0
#             U[m+1-p:m+1] .= 1
#             for j=2:n+1-p
#                 i = Int(floor(j*d))
#                 alpha = j*d-i
#                 U[j+p] = (1-alpha)*ubar[i-1] + alpha*ubar[i]
#             end
#         end
#     end

#     ##-- Set up arrays: N, W, S, T, M
# #! YOU ARE HERE, BIG MATRICES ARE NOT WORKING, LOTS OF ZEROS AND NANS...
#     for k=1:length(Q[1,:]) #do each dimension separately.
#         j = 1
#         mu2 = 1
#         mc2 = 1
#         for i=1:r
#             println("i: ", i)
#             span = Splines.getspanindex(n,p,ubar[i],U)
#             dflag = 0 #derivative flag
#             #check for derivative at point
#             if j <= s
#                 if i == I[j]
#                     dflag = 1
#                 end
#             end

#             if dflag == 0 #if derivative not present
#                 funs = Splines.basisfunctions(span+1, ubar[i], p, U)
#             else
#                 funs = Splines.basisfunctionsderivatives(span+1, ubar[i], p, 1, U)
#             end

#             #if point is unconstrained
#             if Wq[i] > 0
#                 W[mu2,mu2] = Wq[i]
#                 N[mu2,span+1-p:span+1] = funs
#                 S[mu2] = W[mu2,mu2]*Q[i,k]
#                 mu2 += 1
#             else #if point is constrained
#                 M[mc2,span+1-p:span+1] = funs
#                 T[mc2] = Q[i,k]
#                 mc2 += 1
#             end #if unconstrained

#             #if derivative given for this point
#             if dflag == 1
#                 if Wd[j] > 0 #unconstrained derivative
#                     W[mu2,mu2] = Wd[j]
#                     N[mu2,span+1-p:span+1] = funs[2,:]
#                     S[mu2] = W[mu2,mu2]*D[j]
#                     mu2 += 1
#                 else #constrained derivative
#                     M[mc2,span+1-p:span+1] = funs[2,:]
#                     T[mc2] = D[j]
#                     mc2 += 1
#                 end
#                 j += 1
#             end #if dflag
#         end #for r

#         NtransWN = N'*W*N
#         NtransWNinv = LinearAlgebra.inv(NtransWN)
#         NtransWS = N'*W*S

#         if mc < 0 #if no constraints
#             P[:,k] = NtransWN\NtransWS
#         else
#             A = (M*NtransWNinv*M') \ (M*NtransWNinv*NtransWS-T)
#             P[:,k] = NtransWN \ (NtransWS - M'*A)
#         end

#     end #for dimension

#     return U, P

# end

# """
#     getremovalboundcurve(n,p,U,P,u,r,s)

# Compute the knot removal error bounds, \$B_r\$. (NURBS A9.8)

# Inputs:
# - n : there are n+1 control points
# - p : curve degree
# - U : knot vector
# - P : control points
# - u : internal knot
# - r : index of u in U (assuming zero indexing)
# - s : multiplicity of u

# Outputs:
# - Br : distances
# """
# function getremovalboundcurve(n,p,U,P,u,r,s)
#     println("knot: ",r)
#     ord = p+1
#     println("ord: ", ord)
#     lastidx = r-s +1
#     println("last: ", lastidx)
#     firstidx = r-p +1
#     println("first: ", firstidx)
#     off = firstidx-1 #difference in index between temp and P.
#     println("off: ", off)
#     temp = zeros(size(P))
#     println("P[off,:]: ", P[off,:])
#     temp[1,:] = P[off,:]
#     println("P[lastidx+1,:]: ", P[lastidx+1,:])
#     temp[lastidx+1-off,:] = P[lastidx+1,:]
#     i = firstidx
#     j = lastidx
#     ii = 1 +1
#     jj = lastidx-off +1
#     while j-i > 0
#         #compute new control points for one removal step
#         alfi = (u-U[i])/(U[i+ord]-U[i])
#         alfj = (u-U[j])/(U[j+ord]-U[j])
#         temp[ii] = (P[i]-(1-alfi)*temp[ii-1])/alfi
#         temp[jj] = (P[j]-alfj*temp[jj+1])/(1-alfj)
#         i += 1
#         ii += 1
#         j -= 1
#         jj -= 1
#     end #while j-1>0

#     if j-i < 0 #now get bound
#         Br = LinearAlgebra.norm(temp[ii-1]-temp[jj+1]) #eqn 9.82
#     else
#         alfi = (u-U[i])/(U[i+ord]-U[i])
#         Br = LinearAlgebra.norm( P[i] - (alfi*temp[ii+1]+(1-alfi)*temp[ii-1]) ) #eqn 9.80
#     end

#     return Br
# end

# """
#     removeknotsboundcurve(n,p,U,P,ub,ek,E)

# Remove knots from bounded curve. (NURBS A9.9)

# Inputs:
# - n : there are n+1 control points
# - p : curve degree
# - U : knot vector
# - P : control points
# - u : internal knot
# - ub : knot to remove?
# - ek : accumulated error?
# - E : max allowable error

# Outputs:
# - ek : accumulated error
# - nh : new n
# - Uh : new U
# - Ph : new P
# """
# function removeknotsboundcurve(n,p,U,P,ub,ek,E)

#     #get Br values for all distinct interior knots
#     uniqueknot = unique(U)
#     Br = zeros(length(uniqueknot)-2)

#     multiplicity = [(count(U->U==i,U)) for i in uniqueknot]
#     multiplicity = multiplicity[2:end-1]
#     firstidx = zeros(Int,length(Br))
#     for i=1:length(firstidx)
#         firstidx[i] = (LinearIndices(U))[findall(U->U==uniqueknot[i+1],U)][1]
#     end


#     for i=1:length(uniqueknot)-2
#         Br[i] = Splines.getremovalboundcurve(n,p,U,P,uniqueknot[i+1],firstidx[i]-1,multiplicity[i])
#     end

#     #for each basis function, get range of parameter indices.
#     basisindices = zeros(length(P[:,1]),p+2)
#     for i=1:length(P[:,1])
#         basisindices[i,:] = collect(i:p+i+1)
#     end

#     while true
#         m = length(U)-1

#         #find knot with smallest Br bound
#         minBr, minidx = findmin(Br)

#         #set r and s
#         #get unique knots and their mulitplicities
#         uniqueknot = unique(U)
#         multiplicity = [(count(U->U==i,U)) for i in uniqueknot]
#         multiplicity = multiplicity[2:end-1]
#         s = multiplicity[minidx]

#         firstidx = zeros(Int,length(Br))
#         for i=1:length(uniqueknot)-2
#             firstidx[i] = (LinearIndices(U))[findall(U->U==uniqueknot[i+1],U)][1]
#         end
#         r = firstidx[minidx]

#         if minBr==Inf
#             break #finished, you've done all the interior knots
#         end
#         #using eqns 9.81 and 9.83 and A2.4, compute NewError[k], form temp[k] = ek[k] + NewError[k] at all ub[k] values falling within the relevant domain
#         temp = zeros(size(ek))
#         NewError = 0.0
#         for i=1:length(ub)
#             if ub[i] >= U[r]-Br[minidx] && ub[i] <= U[r]+Br[minidx]
#                 N = Splines.singlebasisfunction(p,m,U,i,ub[i])
#                 if mod(p+s,2) == 0
#                     k = Int((p+s)/2)
#                     NewError = N*minBr
#                 else
#                     k = Int(floor((p+s+1)/2))
#                     alpha = (U[r] - U[r-k+1])/(U[r-k+p+2] - U[r-k+1])
#                     NewError = (1-alpha)*N*minBr
#                 end
#                 temp[k] = ek[k] + NewError
#             end
#         end

#         #if knot is removable (all temp[k] <= E)
#         if all(temp .<= E)
#             println("REMOVING")
#             #update ek for relevant range
#             for i=1:length(ub)
#                 if ub[i] >= U[r]-Br[minidx] && ub[i] <= U[r]+Br[minidx]
#                     ek[k] = temp[k]
#                 end
#             end

#             #remove knot (A5.8 without tolerance check)
#             Pw = [P ones(length(P[:,1]))]
#             t, U, Pw = Splines.removecurveknot(n,p,U,Pw,U[r],r,s,num,tolcheck=false)

#             #if no more knots, break
#             if length(U) == 2*(p+1)
#                 break #no more knots to remove.
#             end

#             #using Eqn 9.84, compute new index ranges for affected basis functions
#             for i=r-p-1:r-s
#                 basisindices[i,:] = collect(i:p+i+1)
#             end

#             #using Eqn 9.85 compute new error bounds for the relevant knots
#             for i=1:length(uniqueknot)-2
#                 firstidx[i] = (LinearIndices(U))[findall(U->U==uniqueknot[i+1],U)][1]
#                 if i < max(r-p,p+1) && i > min(r+p-s+1,n)
#                     Br[i] = getremovalboundcurve(n,p,U,P,uniqueknot[i+1],firstidx[i]-1,multiplicity[i+1])
#                 end
#             end

#         else
#             #set this Br to Inf
#             Br[minidx] = Inf
#         end
#     end

#     nh = length(P[:,1])-1
#     return ek,nh,U,P
# end

"""
    f(u,Q,C,dC,ddC)

Auxiliar function for point projection: Solves for f(u) and f'(u) used in eqn 6.3 in NURBS.

Inputs:
- u : knot value
- Q : point to project
- C : curve point at u
- dC : first curve derivative at u
- ddC : second curve derivative at u

Outputs:
- fu : f(u)
- fprimeu : f'(u)
"""
function f(u,Q,C,dC,ddC)

    fu = LinearAlgebra.dot(dC,C-Q)
    fprimeu = LinearAlgebra.dot(ddC,(C-Q)) + LinearAlgebra.norm(dC)^2

    return fu, fprimeu

end

"""
    projectionloop(n,p,U,P,ui,Q,eps1,eps2,a=0.0,b=1.0; closed=false)

Run convergence criteria loop (Newton iteration), returning the knot value that satisfies the criteria.

Inputs:
- n : there are n+1 control points, P
- p : curve degree
- U : knot vector
- P : control points
- ui : the parametric start point
- Q : the point to project
- eps1 : tolerance for points being on the curve
- eps2 : tolerance for points being projected on the curve
- a : the lower bound of the span we're checking
- b : the upper bound of the span we're checking
- closed : boolean whether spline is closed or open

Outputs:
- uip1 : the final knot calculated after criteria are met
"""
function projectionloop(n,p,U,P,ui,Q,eps1,eps2,a=0.0,b=1.0; closed=false)
    local uip1
    while true

        ## Check criteria
        CK = Splines.curvederivatives1(n, p, U, P, ui, 1)
        C = CK[1,:]
        dC = CK[2,:]
        #criteria 1: are they the same point on the spline
        crit1 = LinearAlgebra.norm(C-Q)
        #criteria 2: are they aligned
        num2 = LinearAlgebra.norm(LinearAlgebra.dot(dC,(C-Q)))
        den2 = LinearAlgebra.norm(dC)*LinearAlgebra.norm(C-Q)
        crit2 = num2/den2
        #if criteria 1 or 2 is not met
        if crit1 > eps1 || crit2 > eps2
            ##compute a new uip1
            # update u to be next value
            # println("ui= ", ui)
            #calculate C(u_i) and C'(u_i)
            CK = Splines.curvederivatives1(n, p, U, P, ui, 2)
            C = CK[1,:]
            dC = CK[2,:]
            ddC = CK[3,:]
            # println("ui= ", ui)
            # println("Q= ", Q)
            # println("C= ", C)
            # println("dC= ", dC)
            # println("ddC= ", ddC)

            #calculate f(u_i) and f'(u_i)
            fu, fprimeu = Splines.f(ui,Q,C,dC,ddC)
            #calculate u_{i+1}
            uip1 = ui - fu/fprimeu
            # println("uip1: ", uip1)

            ##check criteria 3 and 4
            #criteria 3: is the point beyond the spline ends? if so, adjust the point
            if closed == false
                if uip1 < a
                    uip1 = a
                elseif uip1 > b
                    uip1 = b
                end
            else
                if uip1 < a
                    uip1 = b - (a - uip1)
                elseif uip1 > b
                    uip1 = a + (uip1 - b)
                end
            end

            #criteria 4: Was the change in points very small?
            crit4 = LinearAlgebra.norm((uip1-ui)*dC)
            # println("it's: ", crit1 <= eps1 , crit2 <= eps2 , crit4 <= eps1)
            #if any of criteria 1, 2, or 4 is satisfied, break.
            if crit1 <= eps1 || crit2 <= eps2 || crit4 <= eps1
                #crit1 means the point is on the spline
                #crit2 means the point is properly projected
                #crit4 means the point is off the end point of the spline.
                break
            end

            ui = uip1

        else #if both were satisfied, then we're done.
            uip1 = ui
            break
        end

    end

    return uip1

end

"""
    projectpoints(n,p,U,P,Q,eps1=1e-10,eps2=1e-10,ncheckvals=100)

Project points, Q, onto spline, returning projected points, R.

Inputs:
- n : there are n+1 control points, P
- p : curve degree
- U : knot vector
- P : control points
- Q : points to project
- eps1 : tolerance for points being on the curve
- eps2 : tolerance for points being projected on the curve
- ncheckvals : number of points used to look for good starting points.

Outputs:
- uproj : knot values on curve where Q has been projected
- R : the projected points.
"""
function projectpoints(n,p,U,P,Q; eps1=eps(),eps2=eps(),ncheckvals=1000) #TODO check what tolerances are good.
    #weight control points
    Pw = [P ones(length(P[:,1]))]

    ##find start values, u0
    #setup a somewhat fine distribution of candidates
    ucheck = collect(range(0,stop=1,length=ncheckvals))

    #get curve points at knot values we're comparing
    Ccheck = zeros(ncheckvals,length(Pw[1,:]))
    for i=1:ncheckvals
        Ccheck[i,:] = Splines.curvepoint(n, p, U, Pw, ucheck[i])
    end
    Ccheck = Ccheck[:,1:end-1]
    #get distances for each candidate start points and find minimum for each point, Q
    u0 = zeros(length(Q[:,1]))
    for i=1:length(Q[:,1])
        mindist = Inf
        for j=1:ncheckvals
            dist = LinearAlgebra.norm(Q[i,:]-Ccheck[j,:])
            if dist < mindist
                mindist = dist
                u0[i] = ucheck[j]
            end
        end
    end

    #for each u0, find the parametric value associated with the point's projection onto the curve.
    uproj = zeros(length(Q[:,1]))
    for i=1:length(uproj)
        uproj[i] = Splines.projectionloop(n,p,U,P,u0[i],Q[i,:],eps1,eps2)
    end

    #return projected curve points, R_k
    Pw = [P ones(length(P[:,1]))]
    R = zeros(length(Q[:,1]),length(Pw[1,:]))
    for i=1:length(uproj)
        R[i,:] = curvepoint(n, p, U, Pw, uproj[i])
    end

    return uproj, R

end


# """
#     globalcurveapproximation(m,Q,p,E; knotplacement)

# Compute global curve approximation of datapoints within bound, E. (NURBS A9.10)

# Inputs:
# - m : there are m+1 data points to approximate
# - Q : data points to approximate
# - p : degree of approximating curve
# - E : maximum error allowance for approximating curve
# - knotplacement : knot placement scheme (centripetal, or chordlength)

# Outputs:
# - n : there are n+1 control points for approximating curve
# - U : knot vector for approximating curve
# - P : Control points of approximating curve
# """
# function globalcurveapproximation(m,Q,p,E; knotplacement="centripetal")

#     #compute ubar and load into ub[]
#     ub = Splines.computeubar(m,Q,knotplacement)
#     #set U and P to be the degree 1 curve interpolating Q
#     U = [0.0; ub; 1.0]
#     P = copy(Q)
#     ek = zeros(m+1)
#     n = m

#     for deg=1:p+1
#         ek,n,U,P = Splines.removeknotsboundcurve(n,deg,U,P,ub,ek,E)
#         if deg == p
#             break
#         end

#         #let U be the knot vector obtained by degree elevating Uh from deg to deg+1 (this is simply increasing the mulitplicities of each knot by 1. No need for degreeelevatecurve function)
#         uniqueknot = unique(U)

#         insert!(U,1,0)
#         knot2add = 0
#         nextidx = 1
#         for i=2:length(uniqueknot)
#             nextidx = findnext(y->y!=knot2add,U,nextidx)
#             knot2add = uniqueknot[i]
#             insert!(U,nextidx,knot2add)
#         end

#         #reset n
#         n = length(U) - deg - 1

#         #fit a least squares curve to the Q_k, using n, ub, degree=deg+1, and new knot vector U to get new control pointsopen
#         #! THIS IS BEING BROKEN RIGHT NOW, NEED TO FIGURE IT OUT...
#         U, P = Splines.leastsquarescurve(Q,m,n,deg+1,ub,U,knotplacement=knotplacement)

#         #project all Q_k to current curve to get R_k = C(u_k).
#         #! YOU ARE HERE
#         R = projectpoints(n,p,U,P,Q)

#         #Update ek and ub
#         for i=1:length(ek)
#             ek[i] = LinearAlgebra.norm(Q[i,:]-R[i,:])
#         end
#         ub = U #? Not sure how this works...


#     end

#     return n, U, P
# end
