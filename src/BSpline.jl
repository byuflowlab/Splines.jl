"""
    BSpline(degree, knots, ctrlpts)

Construct a b-spline object

# Arguments
- `deg::Integer`: degree
- `knots::Vector{Float64}`:: a knot vector (u_0, ... u_n+1) 
- `ctrlpts::Vector{Vector{Float64}}`:: control points.  outer index is number of control points, inner index of dimensionality of point.
"""
struct BSpline{TF, TI}
    degree::TI
    knots::Vector{TF}
    ctrlpts::Vector{Vector{TF}}
end


"""
    getspanindex(deg, knots, u)

(private function) binary search to find span index of vector, knots, in which the parametric point, u, lies. (NURBS A2.1)

# Arguments
- `deg::Integer`: degree
- `knots::Vector{Float64}`:: a knot vector (u_0, ... u_n+1) 
- `u`::Float64`: nondimensional location we are searching for

# Returns
- `span::Integer`: corresponding index i for u between knots_i and knots_i+1
"""
function getspanindex(deg, knots, u)
    
    n = length(knots) - deg - 2
    U = OffsetArray(knots, 0:length(knots)-1)

    if u >= U[n+1]
        return n
    elseif u <= U[deg]
        return deg
    end

    low = deg
    high = n+1
    mid = (low + high) ÷ 2
    while (u < U[mid] || u >= U[mid+1])
        if u < U[mid]
            high = mid
        else
            low = mid
        end
        mid = (low + high) ÷ 2
    end
    return mid
end


"""
    basisfunctions(span, deg, knots, u)

(private function) Compute nonvanishing basis functions (NURBS A2.2)

# Arguments
- `deg::Integer`: degree
- `knots::Vector{Float64}`:: a knot vector (u_0, ... u_n+1) 
- `u`::Float64`: nondimensional location we are searching for
- `span::Integer`: corresponding index i for u between knots_i and knots_i+1 (computed from getspanindex)

# Returns
- `N::Vector{Float64}`: vector of length N_0 ... N_deg
"""
function basisfunctions(span, deg, knots, u)

    N = OffsetArray(zeros(deg+1), 0:deg)
    left = OffsetArray(zeros(deg+1), 0:deg)
    right = OffsetArray(zeros(deg+1), 0:deg)
    U = OffsetArray(knots, 0:length(knots)-1)

    N[0] = 1.0
    for j = 1:deg
        left[j] = u - U[span+1-j]
        right[j] = U[span+j] - u
        saved = 0.0
        for r = 0:j-1
            temp = N[r] / (right[r+1] + left[j-r])
            N[r] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        end
        N[j] = saved
    end
    return N
end


"""
    singlebasisfunction(i, deg, knots, u)

(private function) Compute single basis function N_i^p. (NURBS A2.4)

# Arguments
- `i::Integer` : the index of the basis (the i in N_i)
- `deg::Integer` : the basis degree up to the the curve order
- `knots::Vector{Float}` : the knot vector
- `u::Float` : parametric point of interest

# Returns
- `Nip::Float`: the N_i^p basis function value
"""
function singlebasisfunction(i, deg, knots, u)
    m = length(knots) - 1
    U = OffsetArray(knots, 0:m)
    p = deg

    if (i == 0 && u == U[0]) || (i == m-p-1 && u == U[m])
        return 1.0
    elseif u < U[i] || u >= U[i+p+1]
        return 0.0
    end
    N = OffsetArray(zeros(p+1), 0:p)
    for j = 0:p
        if u >= U[i+j] && u < U[i+j+1]
            N[j] = 1.0
        else
            N[j] = 0.0
        end
    end
    for k = 1:p
        if N[0] == 0.0
            saved = 0.0
        else
            saved = ((u-U[i])*N[0]) / (U[i+k]-U[i])
        end
        for j = 0:p-k
            Uleft = U[i+j+1]
            Uright = U[i+j+k+1]
            if N[j+1] == 0.0
                N[j] = saved
                saved = 0.0
            else
                temp = N[j+1] / (Uright-Uleft)
                N[j] = saved + (Uright-u)*temp
                saved = (u-Uleft)*temp
            end
        end
    end
    return N[0]
end

"""
    basisfunctionsderivatives(span, deg, knots, u, n)

(private function) Calculate the non-vanishing basis functions and derivatives of the B-Spline of order ``p```, 
defined by knots U at parametric point, ``u``. (NURBS A 2.3)

# Arguments
- `span::Integer`: knot span containing u
- `deg::Integer`: the curve order
- `knots::Vector{Float}`: the knot vector
- `u::Float`: parametric point of interest
- `n::Integer` : the max derivative order (n ≦ p)

# Returns
- `ders::Matrix{Float}`: [0..n, 0..p]  ders[0, :] function values, ders[1: :], first derivatives, etc.
"""
function basisfunctionsderivatives(span, deg, knots, u, n)
    
    U = OffsetArray(knots, 0:length(knots)-1)
    p = deg
    i = span
    
    ndu = OffsetArray(zeros(p+1, p+1), 0:p, 0:p)
    a = OffsetArray(zeros(2, p+1), 0:1, 0:p)
    ders = OffsetArray(zeros(n+1, p+1), 0:n, 0:p)
    left = OffsetArray(zeros(p+1), 0:p)
    right = OffsetArray(zeros(p+1), 0:p)
    ndu[0, 0] = 1.0

    for j = 1:p
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved = 0.0
        for r = 0:j-1
            ndu[j, r] = right[r+1] + left[j-r]
            temp = ndu[r, j-1] / ndu[j, r]

            ndu[r, j] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        end
        ndu[j, j] = saved
    end

    for j = 0:p
        ders[0, j] = ndu[j, p]
    end

    for r = 0:p
        s1 = 0; s2 = 1
        a[0, 0] = 1.0
        for k = 1:n
            d = 0.0
            rk = r-k; pk = p-k
            if r >= k
                a[s2, 0] = a[s1, 0] / ndu[pk+1, rk]
                d = a[s2, 0] * ndu[rk, pk]
            end
            j1 = (rk >= -1) ? 1 : -rk
            j2 = (r-1 <= pk) ? k-1 : p-r
            for j = j1:j2
                a[s2, j] = (a[s1, j] - a[s1,j-1]) / ndu[pk+1, rk+j]
                d += a[s2, j] * ndu[rk+j, pk]
            end
            if r <= pk
                a[s2, k] = -a[s1, k-1] / ndu[pk+1, r]
                d += a[s2, k] * ndu[r, pk]
            end
            ders[k, r] = d
            j=s1; s1=s2; s2=j
        end
    end

    r = p
    for k = 1:n
        for j = 0:p
            ders[k, j] *= r
        end
        r *= (p-k)
    end

    return ders
end 



"""
Evaluate point on B-spline curve (NURBS, A3.1)

# Arguments
- `bspline::BSpline: bspline object`
- `u::Float`: point on spline to evaluate at

# Returns
- `C::Vector{Float}`: point in ND space
"""
function curvepoint(bspline::BSpline, u)
    
    P = OffsetArray(bspline.ctrlpts, 0:length(bspline.ctrlpts)-1)
    p = bspline.degree

    span = getspanindex(p, bspline.knots, u)
    N = basisfunctions(span, p, bspline.knots, u)
    C = zeros(length(bspline.ctrlpts[1]))
    for i = 0:p
        C += N[i]*P[span-p+i]
    end
    return C
end



"""
    curvederivatives(bspline, u, d)

Compute a curve point and its derivatives up do the dth derivative at parametric point, ``u``. 
(NURBS, A3.2)

# Arguments
- `bspline::BSpline: bspline object`
- `u::Float`: point on spline to evaluate at
- `d::Float`: derivative order (0 ≤ k ≤ d)

# Returns
- `CK::Vector{Vector{Float}}`.  where CK[0] is the point, CK[1] the first derivative, and so on.
"""
function curvederivatives(bspline::BSpline, u, d)
    p = bspline.degree
    P = OffsetArray(bspline.ctrlpts, 0:length(bspline.ctrlpts)-1)

    du = min(d, p)
    ndim = length(bspline.ctrlpts[1])
    # CK = OffsetArray(zeros(du+1, ndim), 0:du, 1:ndim)
    CK = OffsetArray(Vector{Vector{Float64}}(undef, du+1), 0:du)

    for k = p+1:d
        CK[k] = zeros(ndim)
    end

    span = getspanindex(p, bspline.knots, u)
    nders = basisfunctionsderivatives(span, p, bspline.knots, u, du)

    for k = 0:du
        CK[k] = zeros(ndim)
        for j = 0:p
            CK[k] += nders[k, j] * P[span-p+j]
        end
    end

    return CK  #onebased(CK)
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

# Arguments
- `bspline::BSpline: bspline object`
- `r1::Integer : first control point index to take derivatives at
- `r2::Integer : last control point index to take derivatives at
- `d::Float`: derivative order (0 ≤ k ≦ d)

# Returns
- `PK::Matrix{Vector{Float}}`: PK[i, j] is the ith derivative of the jth control point
"""
function curvederivativecontrolpoints(bspline, r1, r2, d)
    p = bspline.degree
    P = OffsetArray(bspline.ctrlpts, 0:length(bspline.ctrlpts)-1)
    U = OffsetArray(bspline.knots, 0:length(bspline.knots)-1)
    r1 -= 1  # 0-based indexing
    r2 -= 1
    ndim = length(bspline.ctrlpts[1])

    r = r2 - r1

    PK = OffsetArray(Matrix{Vector{Float64}}(undef, d+1, r+1), 0:d, 0:r)
    for i = 0:r
        for k = 0:d
            PK[k, i] = zeros(ndim)
        end
    end

    for i = 0:r
        PK[0, i] = P[r1 + i]
    end

    for k = 1:d
        tmp = p-k+1
        for i = 0:r-k
            PK[k, i] = tmp*(PK[k-1, i+1] - PK[k-1, i]) / (U[r1+i+p+1] - U[r1+i+k])
        end
    end

    return OffsetArray(PK, 0:d, 1:r+1)
end


#There is another curvederivatives algorithm in the book (Algorithm 3.4)


"""
    globalcurveinterpolation(pts, deg)

Interpolate ctrl points pts, with a B-Spline of degree deg. (NURBS A9.1)

# Arguments
- `pts::Vector{Vector{Float}}`: outer vector of length npts, inner vector of length dimension of space
- `deg::Integer`: degree of B-spline

Outputs:
- `bspline::BSpline`: bspline object
"""
function globalcurveinterpolation(pts, deg)

    n = length(pts) - 1
    Q = OffsetArray(pts, 0:n)
    m = n + deg + 1

    d = 0.0
    for i = 1:n
        d += norm(Q[i] - Q[i-1])
    end

    u = OffsetArray(zeros(n+1), 0:n)
    u[0] = 0.0
    u[n] = 1.0
    for i = 1:n-1
        u[i] = u[i-1] + norm(Q[i] - Q[i-1]) / d
    end

    U = OffsetArray(zeros(m+1), 0:m)
    for i = 1:deg
        U[i] = 0.0
    end
    for i = m-deg:m
        U[i] = 1.0
    end

    for j = 1:n-deg
        for i = j:j+deg-1
            U[j+deg] += u[i]
        end
        U[j+deg] /= deg
    end
    knots = onebased(U)

    A = OffsetArray(zeros(n+1, n+1), 0:n, 0:n)
    for i = 0:n
        span = getspanindex(deg, knots, u[i])
        N = basisfunctions(span, deg, knots, u[i])
        A[i, span-deg:span] = N
    end
    A1 = onebased(A)

    ndim = length(Q[0])
    npts = length(Q)
    P = Vector{Vector{Float64}}(undef, npts)
    for j = 1:npts
        P[j] = zeros(ndim)
    end
    for i = 1:ndim
        rhs = [q[i] for q in Q]
        x = A1\onebased(rhs)
        for j = 1:npts
            P[j][i] = x[j]
        end
    end

    return BSpline(deg, knots, P)
end


"""
    leastsquarescurve(pts, nctrl, deg)

Least squares fit to provided points.  (NURBS section 9.4.1)
TODO: Currently hard-coded to 2D data.

# Arguments
- `pts::Vector{Vector{Float}}`: data points
- `nctrl::Integer`: number of control points to use in fit
- `deg::Integer`: degree of bspline used in fit

# Returns
- `bspline::BSpline`: a BSpline object
"""
function leastsquarescurve(pts, nctrl, deg)

    m = length(pts) - 1
    Q = OffsetArray(pts, 0:m)
    n = nctrl - 1

    # --- compute ubar_k -----
    dist = 0.0
    for i = 1:m
        dist += norm(Q[i] - Q[i-1])
    end

    u = OffsetArray(zeros(m+1), 0:m)
    u[0] = 0.0
    u[m] = 1.0
    for i = 1:m-1
        u[i] = u[i-1] + norm(Q[i] - Q[i-1]) / dist
    end

    # ---- compute knot vector -----
    r = n + deg + 1
    U = OffsetArray(zeros(r+1), 0:r)
    for i = 1:deg
        U[i] = 0.0
    end
    for i = r-deg:r
        U[i] = 1.0
    end

    d = (m + 1) / (n - deg + 1)
    for j = 1:n-deg
        i = floor(Int, j * d)
        alpha = j * d - i
        U[deg + j] = (1 - alpha)*u[i - 1] + alpha*u[i]
    end
    knots = OffsetArrays.no_offset_view(U)

    # --------- setup N matrices -------
    Nm = zeros(m-1, n-1)
    for j = 1:n-1
        for i = 1:m-1
            Nm[i, j] = singlebasisfunction(j, deg, knots, u[i])
        end
    end

    # ------ Rk --------
    Rk = zeros(m-1, length(Q[0]))
    for k = 1:m-1
        N0 = singlebasisfunction(0, deg, knots, u[k])
        Nn = singlebasisfunction(n, deg, knots, u[k])
        Rk[k, :] = Q[k] - N0*Q[0] - Nn*Q[m]
    end

    # ------ R ---------
    ndim = length(Q[0])
    R = zeros(n-1, ndim)
    # Rx = zeros(n-1)
    # Ry = zeros(n-1)
    for i = 1:m-1
        for j = 1:n-1
            N = singlebasisfunction(j, deg, knots, u[i])
            for k = 1:ndim
                R[j, k] += N * Rk[i, k]
            end
            # Rx[j] += N * Rk[i, 1]
            # Ry[j] += N * Rk[i, 2]
        end
    end    

    # --- ctrlpts -----
    NN = Nm'*Nm

    P = Vector{Vector{Float64}}(undef, nctrl)
    for j = 1:nctrl
        P[j] = zeros(ndim)
    end
    for i = 1:ndim
        Px = NN \ R[:, i]
        Px = [Q[0][i]; Px; Q[m][i]]
        for j = 1:nctrl
            P[j][i] = Px[j]
        end
    end

    # Px = NN \ Rx
    # Py = NN \ Ry

    # Px = [Q[0][1]; Px; Q[m][1]]
    # Py = [Q[0][2]; Py; Q[m][2]]
    # P = [[p[1], p[2]] for p in zip(Px, Py)]

    return BSpline(deg, knots, P)
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

# """
#     f(u,Q,C,dC,ddC)

# Auxiliar function for point projection: Solves for f(u) and f'(u) used in eqn 6.3 in NURBS.

# Inputs:
# - u : knot value
# - Q : point to project
# - C : curve point at u
# - dC : first curve derivative at u
# - ddC : second curve derivative at u

# Outputs:
# - fu : f(u)
# - fprimeu : f'(u)
# """
# function f(u,Q,C,dC,ddC)

#     fu = LinearAlgebra.dot(dC,C-Q)
#     fprimeu = LinearAlgebra.dot(ddC,(C-Q)) + LinearAlgebra.norm(dC)^2

#     return fu, fprimeu

# end

# """
#     projectionloop(n,p,U,P,ui,Q,eps1,eps2,a=0.0,b=1.0; closed=false)

# Run convergence criteria loop (Newton iteration), returning the knot value that satisfies the criteria.

# Inputs:
# - n : there are n+1 control points, P
# - p : curve degree
# - U : knot vector
# - P : control points
# - ui : the parametric start point
# - Q : the point to project
# - eps1 : tolerance for points being on the curve
# - eps2 : tolerance for points being projected on the curve
# - a : the lower bound of the span we're checking
# - b : the upper bound of the span we're checking
# - closed : boolean whether spline is closed or open

# Outputs:
# - uip1 : the final knot calculated after criteria are met
# """
# function projectionloop(n,p,U,P,ui,Q,eps1,eps2,a=0.0,b=1.0; closed=false)
#     local uip1
#     while true

#         ## Check criteria
#         CK = Splines.curvederivatives1(n, p, U, P, ui, 1)
#         C = CK[1,:]
#         dC = CK[2,:]
#         #criteria 1: are they the same point on the spline
#         crit1 = LinearAlgebra.norm(C-Q)
#         #criteria 2: are they aligned
#         num2 = LinearAlgebra.norm(LinearAlgebra.dot(dC,(C-Q)))
#         den2 = LinearAlgebra.norm(dC)*LinearAlgebra.norm(C-Q)
#         crit2 = num2/den2
#         #if criteria 1 or 2 is not met
#         if crit1 > eps1 || crit2 > eps2
#             ##compute a new uip1
#             # update u to be next value
#             # println("ui= ", ui)
#             #calculate C(u_i) and C'(u_i)
#             CK = Splines.curvederivatives1(n, p, U, P, ui, 2)
#             C = CK[1,:]
#             dC = CK[2,:]
#             ddC = CK[3,:]
#             # println("ui= ", ui)
#             # println("Q= ", Q)
#             # println("C= ", C)
#             # println("dC= ", dC)
#             # println("ddC= ", ddC)

#             #calculate f(u_i) and f'(u_i)
#             fu, fprimeu = Splines.f(ui,Q,C,dC,ddC)
#             #calculate u_{i+1}
#             uip1 = ui - fu/fprimeu
#             # println("uip1: ", uip1)

#             ##check criteria 3 and 4
#             #criteria 3: is the point beyond the spline ends? if so, adjust the point
#             if closed == false
#                 if uip1 < a
#                     uip1 = a
#                 elseif uip1 > b
#                     uip1 = b
#                 end
#             else
#                 if uip1 < a
#                     uip1 = b - (a - uip1)
#                 elseif uip1 > b
#                     uip1 = a + (uip1 - b)
#                 end
#             end

#             #criteria 4: Was the change in points very small?
#             crit4 = LinearAlgebra.norm((uip1-ui)*dC)
#             # println("it's: ", crit1 <= eps1 , crit2 <= eps2 , crit4 <= eps1)
#             #if any of criteria 1, 2, or 4 is satisfied, break.
#             if crit1 <= eps1 || crit2 <= eps2 || crit4 <= eps1
#                 #crit1 means the point is on the spline
#                 #crit2 means the point is properly projected
#                 #crit4 means the point is off the end point of the spline.
#                 break
#             end

#             ui = uip1

#         else #if both were satisfied, then we're done.
#             uip1 = ui
#             break
#         end

#     end

#     return uip1

# end

# """
#     projectpoints(n,p,U,P,Q,eps1=1e-10,eps2=1e-10,ncheckvals=100)

# Project points, Q, onto spline, returning projected points, R.

# Inputs:
# - n : there are n+1 control points, P
# - p : curve degree
# - U : knot vector
# - P : control points
# - Q : points to project
# - eps1 : tolerance for points being on the curve
# - eps2 : tolerance for points being projected on the curve
# - ncheckvals : number of points used to look for good starting points.

# Outputs:
# - uproj : knot values on curve where Q has been projected
# - R : the projected points.
# """
# function projectpoints(n,p,U,P,Q; eps1=eps(),eps2=eps(),ncheckvals=1000) #TODO check what tolerances are good.
#     #weight control points
#     Pw = [P ones(length(P[:,1]))]

#     ##find start values, u0
#     #setup a somewhat fine distribution of candidates
#     ucheck = collect(range(0,stop=1,length=ncheckvals))

#     #get curve points at knot values we're comparing
#     Ccheck = zeros(ncheckvals,length(Pw[1,:]))
#     for i=1:ncheckvals
#         Ccheck[i,:] = Splines.curvepoint(n, p, U, Pw, ucheck[i])
#     end
#     Ccheck = Ccheck[:,1:end-1]
#     #get distances for each candidate start points and find minimum for each point, Q
#     u0 = zeros(length(Q[:,1]))
#     for i=1:length(Q[:,1])
#         mindist = Inf
#         for j=1:ncheckvals
#             dist = LinearAlgebra.norm(Q[i,:]-Ccheck[j,:])
#             if dist < mindist
#                 mindist = dist
#                 u0[i] = ucheck[j]
#             end
#         end
#     end

#     #for each u0, find the parametric value associated with the point's projection onto the curve.
#     uproj = zeros(length(Q[:,1]))
#     for i=1:length(uproj)
#         uproj[i] = Splines.projectionloop(n,p,U,P,u0[i],Q[i,:],eps1,eps2)
#     end

#     #return projected curve points, R_k
#     Pw = [P ones(length(P[:,1]))]
#     R = zeros(length(Q[:,1]),length(Pw[1,:]))
#     for i=1:length(uproj)
#         R[i,:] = curvepoint(n, p, U, Pw, uproj[i])
#     end

#     return uproj, R

# end

