"""
    NURBS(degree, knots, weights, ctrlpts)

Construct a NURBS object

# Arguments
- `deg::Integer`: degree
- `knots::Vector{Float}`:: a knot vector (u_0, ... u_n+1)
- `weights::Vector{Float}`:: a corresponding vector of weights
- `ctrlpts::Vector{Vector{Float}}`:: control points.  outer index is number of control points, inner index of dimensionality of point.
"""
struct NURBS{TI, TF1, TF2, TF3, TA<:AbstractVector{TF3}}
    degree::TI
    knots::Vector{TF1}
    weights::Vector{TF2}
    ctrlpts::Vector{TA}
end

"""
    get_degree(spline)

return degree of spline
"""
function get_degree(nurbs)
    return nurbs.degree
end

"""
    get_knots(spline)

return knot vector of spline
"""
function get_knots(nurbs)
    return nurbs.knots
end

"""
    get_weights(nurbs)

return weights of NURBS spline
"""
function get_weights(nurbs)
    return nurbs.weights
end

"""
    get_ctrlpts(spline)

return control points of spline
"""
function get_ctrlpts(nurbs)
    return nurbs.ctrlpts
end

"""
    nurbsbasis(nurbs, u, d)

Get rational basis functions and derivatives.

```math
R_{i,p}(u) = \\frac{N_{i,p}(u)w_i}{\\sum_{j=0}^n N_{j,p}(u)w_j}
```

where `` N_{i,p}(u ) `` are B-Spline Basis Functions and `` w_i `` are weights associated with the NURBS control points.

(see NURBS eqn 4.2)

Inputs:

- `nurbs::NURBS`: nurbs object
- `u::Float`: parametric point of interest
- `d::Integer`: the max derivative order (n ≦ p)

Outputs:

- `R::Vector{Float}`: array of basis function values at the point u.
- `dR::Vector{Float}`: vector of basis function derivative values at point u.

"""
function nurbsbasis(nurbs,u,d)

    # Get spline components
    p = nurbs.degree
    U = zerobased(nurbs.knots)
    w = zerobased(nurbs.weights)

    #get the span index
    i = Splines.getspanindex(p, U, u)

    #get B-Spline basis functions and derivatives for numerator
    bases = Splines.basisfunctionsderivatives(i,p,U,u, d)
    #separate out the bases and their derivatives
    N = bases[0,:]
    dN = bases[1,:] #only doing first derivates right now

    #number of non-zero basis functions at parametric point, u.
    numbasisfunctions = length(N)-1
    Nw = 0
    dNw = 0

    #calculate the denomenator values
    for j=0:numbasisfunctions
        Nw += N[j] .* w[j]
        dNw += dN[j] .* w[j]
    end

    #Calculate each of the non-zero Rational Basis functions and first derivatives
    R = zerobased(zeros(length(N)))
    dR = zerobased(zeros(size(N)))
    for j=0:numbasisfunctions
        R[j] = N[j]/Nw .* w[j]
        dR[j] = w[j] .* (Nw*dN[j] - dNw*N[j]) ./ Nw^2
    end

    return onebased(R), onebased(dR)
end



"""
Evaluate point on rational b-spline curve (NURBS, A4.1)

# Arguments
- `nurbs::NURBS`: NURBS object
- `u::Float`: point on spline to evaluate at

# Returns
- `C::Vector{Float}`: point in ND space
"""
function curvepoint(nurbs::NURBS, u)

    if u > nurbs.knots[end] || u < nurbs.knots[1]
        error("parametric point, u, is outside of knot range, U")
    end

    n = length(nurbs.ctrlpts) - 1
    # P = OffsetArray(nurbs.ctrlpts, 0:n)
    w = nurbs.weights
    Pw = [[p*w[i]; w[i]] for (i, p) in enumerate(nurbs.ctrlpts)]
    Pw = OffsetArray(Pw, 0:n)

    p = nurbs.degree

    span = getspanindex(p, nurbs.knots, u)
    N = basisfunctions(span, p, nurbs.knots, u)
    Cw = zeros(length(Pw[1]))
    for i = 0:p
        Cw += N[i]*Pw[span-p+i]
    end
    return Cw[1:end-1] / Cw[end]
end


"""
    curvederivatives(nurbs, u, d)

Compute a curve point and its derivatives up do the dth derivative at parametric point, ``u``.
(NURBS, A3.2)

# Arguments
- `nurbs::NURBS: bspline object`
- `u::Float`: point on spline to evaluate at
- `d::Float`: derivative order (0 ≤ k ≤ d)

# Returns
- `CK::Vector{Vector{Float}}` where CK[1] is the point, CK[2] the first derivative, and so on.
"""
function curvederivatives(nurbs::NURBS{<:Any, T1, T2, T3}, u::T4, d) where {T1, T2, T3, T4}

    if u > nurbs.knots[end] || u < nurbs.knots[1]
        error("parametric point, u, is outside of knot range, U")
    end

    T = promote_type(T1, T2, T3, T4)

    # evaluate derivatives of A and w
    w = nurbs.weights
    Pw = [[p*w[i]; w[i]] for (i, p) in enumerate(nurbs.ctrlpts)]
    bsp = BSpline(nurbs.degree, nurbs.knots, Pw)
    CK = curvederivatives(bsp, u, d)
    Aders = [c[1:end-1] for c in CK]
    wders = [c[end] for c in CK]

    # initialize
    du = min(d, nurbs.degree)
    CK = OffsetArray(Vector{Vector{T}}(undef, du+1), 0:du)

    # algorithm A 4.2
    for k = 0:d
        v = Aders[k]
        for i = 1:k
            v -= binomial(k, i)*wders[i]*CK[k-i]
        end
        CK[k] = v / wders[0]
    end

    return onebased(CK)
end


"""
    curveknotinsertion(nurbs::NURBS, u, r)

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
- `nurbs::NURBS`: nurbs object
- `u::Float`: the knot to be added
- `r::Integer`: number of times the new knot is inserted (it is assumed that `` r+s \\leq p `` )

Outputs:

if separated_outputs == true
    - `nq::Integer`: the number of control points minus 1 (the index of the last control point) after insertion
    - `UQ::Vector{Float}`: the knot vector after insertion
    - `Qw::Vector{Vector{Float}}`: the set of weighted control points and weights after insertion
else
    -`nurbsout::NURBS`: new nurbs object
"""
function curveknotinsertion(nurbs::NURBS, u, r; separated_outputs=false)

    # Get spline components
    p = nurbs.degree
    UP = zerobased(nurbs.knots)
    P = zerobased(nurbs.ctrlpts)
    w = zerobased(nurbs.weights)
    #get weighted control points
    Pw = zerobased([zeros(length(P[1, :][1])) for _ in 1:length(P)])
    for i=0:length(Pw[:,1])-1
        Pw[i] = [P[i].*w[i]; w[i]]
    end

    #define iteration bounds for algorithm
    np = length(P[:,1])-1
    mp = np+p+1
    nq = np+r

    # get index where to insert
    k = getspanindex(p, UP, u)
    # count current number of repeated instances of knot to insert
    s = count(i->i==u,UP)

    #initialize output vectors
    UQ = zerobased(zeros(length(UP)+r))
    Qw = zerobased([zeros(length(Pw[1, :][1])) for _ in 1:nq+1])
    Rw = zerobased([zeros(length(Pw[1, :][1])) for _ in 1:p+1])


    #TODO: below this point, change all i+1 to just i's
    #Load new knot vector
    for i=0:k
        UQ[i] = UP[i]
    end

    for i=1:r
        UQ[k+i] = u
    end

    for i=k+1:mp
        UQ[i+r] = UP[i]
    end

    #Save unaltered control points
    for i=0:k-p
        Qw[i] = Pw[i]
    end

    for i=k-s:np
        Qw[i+r] = Pw[i]
    end

    for i=0:p-s
        Rw[i] = Pw[k-p+i]
    end
    #insert new knot r times
    L = 0.0 #initialize L in this scope
    for j=1:r
        L = k-p+j
        for i=0:p-j-s
            alpha = (u-UP[L+i])/(UP[i+k+1]-UP[L+i])
            Rw[i] = alpha*Rw[i+1] + (1.0-alpha)*Rw[i]
        end
        Qw[L] = Rw[0]
        Qw[k+r-j-s] = Rw[p-j-s]
    end
    #Load remaining control points
    for i=L+1:k-s
        Qw[i] = Rw[i-L]
    end

    if separated_outputs
        return nq, onebased(UQ), onebased(Qw)
    else
        nurbsout = NURBS(p,onebased(UQ),getindex.(onebased(Qw),3),getindex.(onebased(Qw),[1:2]))
        return nurbsout
    end
end

"""
    refineknotvectorcurve(nurbs::NURBS, X)

Refine curve knot vector using NURBS A5.4.

This algorithm is simply a knot insertion algorithm that allows for multiple knots to be added simulataneously, i.e., a knot refinement procedure.

Inputs:
- `nurbs::NURBS`: nurbs object
- `X::Vector{Float}`: elements, in ascending order, to be inserted

Outputs:

if separated_outputs == true
    - Ubar : the knot vector after insertion
    - Qw : the set of weighted control points and weights after insertion
else
    -`nurbsout::NURBS`: new nurbs object
"""
function refineknotvectorcurve(nurbs::NURBS, X;separated_outputs=false)


     # Get spline components
     p = nurbs.degree
     U = zerobased(nurbs.knots)
     P = zerobased(nurbs.ctrlpts)
     w = zerobased(nurbs.weights)
     #get weighted control points
     Pw = zerobased([zeros(length(P[1, :][1])) for _ in 1:length(P)])
     for i=0:length(Pw)-1
         Pw[i] = [P[i].*w[i]; w[i]]
     end

     #zero base input knots
    X = zerobased(X)
    r = length(X)-1

    #set up indices
    n = length(P)-1
    m = n+p+1
    a = Splines.getspanindex(p,U,X[0])
    b = Splines.getspanindex(p,U,X[r])
    b += 1

    #initialize outputs
    Qw = zerobased([zeros(length(Pw[1, :][1])) for _ in 1:length(Pw)+r+1])
    Ubar = zerobased(zeros(length(U)+r+1))

    for j=0:a-p
        Qw[j] = Pw[j]
    end

    for j=b-1:n
        Qw[j+r+1] = Pw[j]
    end

    for j=0:a
        Ubar[j] = U[j]
    end

    for j=b+p:m
        Ubar[j+r+1] = U[j]
    end
    i = b+p-1
    k = b+p+r

    for j=r:-1:0
        while X[j]<=U[i] && i>a
            Qw[k-p-1,:] = Pw[i-p-1,:]
            Ubar[k] = U[i]
            k -=1
            i -=1
        end
        Qw[k-p-1,:] = Qw[k-p,:]
        for ell=1:p
            ind = k-p+ell
            alpha = Ubar[k+ell] - X[j]
            if abs(alpha)==0.0
                Qw[ind-1,:] = Qw[ind,:]
            else
                alpha /= Ubar[k+ell] - U[i-p+ell]
                Qw[ind-1,:] = alpha*Qw[ind-1,:] + (1.0-alpha)*Qw[ind,:]
            end
        end
        Ubar[k] = X[j]
        k -= 1
    end

    if separated_outputs
        return onebased(Ubar), onebased(Qw)
    else
        nurbsout = NURBS(p,onebased(Ubar),getindex.(onebased(Qw),3),getindex.(onebased(Qw),[1:2]))
        return nurbsout
    end

    return
end

# """
#     removecurveknot(n,p,U,Pw,u,r,s,num;d,tolcheck)

# Remove knot u, num number of times (if possible). (NURBS A5.8)
# Inputs:

# - n : there are n+1 control points
# - p : degree of curve
# - U : knot vector
# - Pw : weighted control points
# - u : knot to be removed
# - r : index of first repetition of knot in U
# - s : multiplicity of u in U
# - num : the number of times desired to remove knot.
# - d : bound on deviation
# - tolcheck: boolean for whether to do tolerance check (we don't do check when this function is being called by the curve fitting algorithm)

# Outputs:
# - t : number of times the knot was actually removed
# - Uhat : new knot vector
# - Phat : new control points
# """
# function removecurveknot(n,p,U,Pw,u,r,s,num; d=1e-6, tolcheck=true)

#     Phat = copy(Pw)
#     Uhat = copy(U)

#     local maxP, minw, i, j, t
#     maxP = 0.0
#     minw = 999.9
#     for i=1:length(Phat[:,1])
#         Pdist = LinearAlgebra.norm(Phat[i,1:end-1])
#         if Pdist > maxP
#             maxP = Pdist
#         end
#         if Phat[i,end] < minw
#             minw = Phat[i,end]
#         end
#     end
#     # println("minw = ", minw)
#     # println("maxP = ", maxP)
#     tol = d*minw/(1+maxP)
#     # println("tol = ", tol)
#     m = n+p+1
#     # println("m = ", m)
#     ord = p+1
#     # println("ord = ", ord)
#     fout = Int(round((2*r-s-p)/2)) #first control point output
#     # println("fout = ", fout)
#     last = r-s+1
#     # println("last = ", last)
#     first = r-p+1
#     # println("first = ", first)

#     #this loop is eqn 5.28
#     temp=zeros(size(Phat))
#     remflag = false
#     for outer t = 0:num
#         # println("t = ", t)
#         off = first-1 #diff in index between temp and P
#         # println("off = ", off)
#         temp[1,:] = Phat[off,:]
#         # println("temp[1,:] = ", temp[1,:])
#         temp[last+1-off,:] = Phat[last+1,:]
#         # println("temp[$(last+1-off),:] = ", temp[last+1-off])
#         i = first
#         # println("i = ", i)
#         j = last
#         # println("j = ", j)
#         ii = 1+1
#         # println("ii = ", ii)
#         jj = last-off+1
#         # println("jj = ", jj)
#         remflag = false

#         while j-i>t
#             # println("j-i>t")
#             #compute new control points for one removal step
#             alfi = (u-Uhat[i])/(Uhat[i+ord+t]-Uhat[i])
#             # println("alfi = ", alfi)
#             alfj = (u-Uhat[j-t])/(Uhat[j+ord]-Uhat[j-t])
#             # println("alfj = ", alfj)
#             temp[ii,:] = (Phat[i,:] - (1.0-alfi)*temp[ii-1,:])/alfi
#             # println("temp[ii,:] = ", temp[ii,:])
#             temp[jj,:] = (Phat[j,:]-alfj*temp[jj+1,:])/(1.0-alfj)
#             # println("temp[jj,:] = ", temp[jj,:])
#             i += 1
#             # println("i = ", i)
#             ii += 1
#             # println("ii = ", ii)
#             j -= 1
#             # println("j = ", j)
#             jj -= 1
#             # println("jj = ", jj)
#         end #while
#         if j-i < t #check if knot is removable
#             # println("j-i < t")
#             if LinearAlgebra.norm(temp[ii-1,:]-temp[jj+1,:]) <= tol || tolcheck==false
#                 remflag = true
#             end
#         else
#             alfi = (u-Uhat[i])/(Uhat[i+ord+t]-Uhat[i])
#             # println("alfi = ", alfi)
#             if LinearAlgebra.norm(Phat[i,:]-(alfi*temp[ii+t+1,:]+(1.0-alfi)*temp[ii-1,:])) <= tol || tolcheck==false
#                 remflag = true
#             end
#         end
#         # println(remflag)
#         if remflag == false #cannot remove any more knots
#             break #get out of for loop
#         else
#             #successful removal. save new cont. pts
#             i = first
#             # println("i = ", i)
#             j = last
#             # println("j = ", j)
#             while j-i > t
#                 # println("j-i > t")
#                 Phat[i,:] = temp[i-off,:]
#                 Phat[j,:] = temp[j-off,:]
#                 i += 1
#                 j -= 1
#             end #while
#         end #if
#         first -= 1
#         last += 1
#     end #for

#     #if no knots removed, end
#     if remflag==false
#         # println("No Knots Removed.")
#         return 0, Uhat, Phat
#     end

#     #shift knots
#     for k = r+1:m
#         Uhat[k-t] = Uhat[k]
#         j = fout
#         i = j #Pj through Pi will be overwritten
#         for k=1:t
#             if mod(k,2) == 1 #k modulo 2
#                 i += 1
#             else
#                 j -= 1
#             end
#         end
#     end

#     #shift
#     for k=i+1:n
#         Phat[j,:] = Phat[k,:]
#         j += 1
#     end

#     return t, Uhat, Phat

# end


"""
    degreeelevatecurve(nurbs::NURBS,t)

Raise degree of spline from p to p `` +t``, `` t \\geq 1 `` by computing the new control point vector and knot vector.

Knots are inserted to divide the spline into equivalent Bezier Curves. These curves are then degree elevated using the following equation.

```math
\\mathbf{P}^t_i = \\sum^{\\textrm{min}(p,i)}_{j=\\textrm{max}(0,i-t)} \\frac{\\binom{p}{j} \\binom{t}{i-j} \\mathbf{P}_j}{\\binom{p+t}{i}}~,~~~~~i=0,...,p+t
```
where `` \\mathbf{P}^t_i `` are the degree elevated control points after `` t `` -degree elevations

Finally, the excess knots are removed and the degree elevated spline is returned.

(see NURBS eqn 5.36, A5.9)

Inputs:
- `nurbs::NURBS`: nurbs object
- `t::Integer`: the number of degrees to elevate, i.e. the new curve degree is p+t

Outputs:

if separated_outputs == true
    - nh : the number of control points minus 1 (the index of the last control point) after degree elevation
    - Uh : the knot vector after degree elevation
    - Qw : the set of weighted control points and weights after degree elevation

else
    -`nurbsout::NURBS`: new nurbs object
"""
function degreeelevatecurve(nurbs::NURBS,t;separated_outputs=false)

    # Get spline components
    p = nurbs.degree
    U = zerobased(nurbs.knots)
    P = zerobased(nurbs.ctrlpts)
    w = zerobased(nurbs.weights)
    #get weighted control points
    Pw = zerobased([zeros(length(P[1, :][1])) for _ in 1:length(P)])
    for i=0:length(Pw)-1
        Pw[i] = [P[i].*w[i]; w[i]]
    end

    #set up indices
    n = length(P)-1
    m = n+p+1
    ph = p+t
    ph2 = Int(floor(ph/2))

    #initialize in outer scope:
    oldr = 0
    ub = 0

    #initialize local arrays
    bezalfs = zerobased([zerobased(zeros(p+1)) for _ in 1:p+t+1]) #coefficients for degree elevating the bezier segments
    bpts = zerobased([zeros(length(Pw[1,:][1])) for _ in 1:p+1]) #pth-degree bezier control points of the current section
    ebpts = zerobased([zeros(length(Pw[1,:][1])) for _ in 1:p+t+1]) #(p+t)th-degree bezier control points of the current section
    nextbpts = zerobased([zeros(length(Pw[1,:][1])) for _ in 1:p-1]) #leftmost control points of next bezier segement
    alfs = zerobased(zeros(p-1)) #knot insertion alphas.

    #initialize outputs:
    s = length(unique(U)) - 2 #number of unique internal knots


    Uh = zerobased(ones( length(U) + t*(s+2) )) #mhat eqn 5.33 is m + (s+2)xt
    Qw = zerobased([zeros(length(Pw[1, :][1])) for _ in 1:length(Pw[:,1]) + t*(s+1)])
    #nhat eqn 5.32 is n+(s+1)xt

    #compute bezier degree elevation coefficients
    bezalfs[ph][p] = 1.0 #bezalfs are coefficients for degree elevating the Bezier segments
    bezalfs[0][0] = 1.0
    for i=1:ph2
        inv = 1.0/Splines.binomialcoeff(ph,i)
        mpi = Int(min(p,i))
        for j=Int(max(0,i-t)):mpi
            bezalfs[i][j] = inv*Splines.binomialcoeff(p,j)*Splines.binomialcoeff(t,i-j) #from eqn 5.36
        end
    end
    for i = ph2+1:ph-1
        mpi = Int(min(p,i))
        for j=Int(max(0,i-t)):mpi
            bezalfs[i][j] = bezalfs[ph-i][p-j]
        end
    end
    mh = ph
    kind = ph+1
    r = -1
    a = p
    b = p+1
    cind = 1
    ua = U[0]

    #first new control point is first old control point
    Qw[0] = copy(Pw[0])

    #fill in first p+t new U vector points (same as old vector with +t multiplicity)
    for i=0:ph
        Uh[i] = ua
    end

    #initialize first bezier segment
    for i=0:p
        bpts[i] = Pw[i] #bpts are the pth-degree bezier control points of the current segment
    end

    #big loop through knot vector starting after first elements of Qw and Uh vectors.
    while b<m
        i = b
        while b<m && U[b] == U[b+1] #incrememnt b for knots with multiplicities.
            b += 1
        end
        mul = b-i+1
        mh += mul + t
        ub = U[b]
        oldr = r
        r = p-mul

        #insert knot u(b) r times
        if oldr > 0
            lbz = Int(floor((oldr + 2)/2))
        else
            lbz = 1
        end
        if r>0
            rbz = Int(ph-floor((r+1)/2))
        else
            rbz = ph
        end

        #insert knot to get bezier segment
        if r>0
            numer = ub - ua
            for k=p:-1:mul+1
                alfs[k-mul-1] = numer/(U[a+k]-ua) #knot insertion alphas
            end
            for j=1:r
                save = r-j
                s = mul+j
                for k=p:-1:s
                    bpts[k] = alfs[k-s]*bpts[k] + (1.0-alfs[k-s])*bpts[k-1]
                end #for
                nextbpts[save] = bpts[p] #nextbpts are the leftmost control points of the next bezier segment
            end #for
        end #insert knot if
        #degree elevate bezier
        for i=lbz:ph
            #only points lbz,...,ph are used below
            ebpts[i] .= 0.0 #ebpts are the (p+t)th-degree bezier control points of the current segment.
            mpi = Int(min(p,i))
            for j = Int(max(0,i-t)):mpi
                ebpts[i] += bezalfs[i][j]*bpts[j]
            end
        end #for; degree elevation

        #must remove knot u=U[a] oldr times
        if oldr > 1
            front = kind-2
            back = kind
            den = ub-ua
            bet = (ub-Uh[kind-1])/den

            #knot removal loop
            if oldr-1 >= 1 #need to make sure this is true since julia does a loop no matter what, when C evaluates condition first.
                for tr=1:oldr-1
                    i = front
                    j = back
                    kj = j-kind+1

                    #loop and compute the new control points for one removal step
                    while j-i > tr
                        if i < cind
                            alf = (ub - Uh[i])/(ua-Uh[i])
                            Qw[i] = alf.*Qw[i] + (1.0-alf).*Qw[i-1]
                        end #if
                        if j >= lbz
                            if j-tr <= kind-ph+oldr
                                gam = (ub-Uh[j-tr])/den
                                ebpts[kj] = gam.*ebpts[kj]+(1.0-gam).*ebpts[kj+1]
                            else
                                ebpts[kj] = bet.*ebpts[kj]+(1.0-bet).*ebpts[kj+1]
                            end #if
                        end #if
                        i += 1
                        j -= 1
                        kj -= 1
                    end #while
                    front -= 1
                    back += 1
                end #for tr
            end #if oldr is big enough
        end #if; remove knot u=U[a]

        #load the knot ua
        if a != p
            for i=0:ph-oldr-1
                Uh[kind] = ua
                kind += 1
            end #for
        end #if

        #load ctrl pts into Qw
        for j=lbz:rbz
            Qw[cind] = copy(ebpts[j])
            cind += 1
        end

        #set up for next pass through loop
        if b < m
            for j=0:r-1
                bpts[j] = nextbpts[j]
            end
            for j=r:p
                bpts[j] = Pw[b-p+j]
            end
            a = b
            b += 1
            ua = ub
        else #end knot
            for i=0:ph
                Uh[kind+i] = ub
            end
        end

    end #while b<m
    nh = mh-ph-1

    if separated_outputs
        return nh, onebased(Uh), onebased(Qw)
    else
        nurbsout = NURBS(p+t,onebased(Uh),getindex.(onebased(Qw),3),getindex.(onebased(Qw),[1:2]))
        return nurbsout
    end


end


# # """
# # """
# # function degreereducecurve(n,p,U,Qw,nh,Uh,Pw)
# #     ph = p-1
# #     mh = ph
# #     kind = pu+1
# #     r = -1
# #     a = p
# #     b = p+1
# #     ind = 1
# #     mult = p
# #     m = n+p+2
# #     Pw[0] = Qw[0]
# #     for i=0:ph #compute left end of knot vecotr
# #         Uh[i] = U[0]
# #     end
# #     for i=0:p #initialize first bezier segment
# #         bpts[i] = Qw[i]
# #     end
# #     for i=0:m-1 #initialize error vector
# #         e[i] = 0.0
# #     end
# #     #loop through the knot vector
# #     while b<m
# #         #first compute knot multiplicity
# #         i = b
# #         while b<m && U[b] == U[b+1]
# #         b += 1
# #         end
# #     mult = b-i+1
# #     mh += mult-1
# #     oldr = r
# #     r = p-mult
# #     if older > 0
# #         lbz = (oldr+2)/2
# #     else
# #         lbz = 1
# #     end
# #     #insert knot U[b] r times
# #     if r>0
# #         numer = U[b] - U[a]
# #         for k=p:-1:mult
# #             alphas[k-mult-1] = numer/(U[a+k]-U[a])
# #         end
# #         for j=1:r
# #             save = r-j
# #             s= mult+j
# #             for k=p:-1:s
# #                 bpts[k] = alphas[k-s]*bpts[k] + (1.0-alphas[k-s])*bpts[k-1]
# #             end
# #             nextbpts[save] = bpts[p]
# #         end
# #     end
# #     #degree reduce bezier segement
# #     BezDegreeReduce(bpts,rbpts,MaxErr)
# #     e[a] = e[a]+MaxErr
# #     if e[a] > tol
# #         return 1 #curve not degree reducible
# #         #remove knot U[a] oldr times
# #     end
# #     if oldr > 0
# #         first = kind
# #         last = kind
# #         for k=0:oldr-1
# #             i = first
# #             j = last
# #             kj = j-kind
# #             while j-i > k
# #                 alfa = (U[a]-Uh[i-1])/(U[b] - Uh[j-k-1])
# #                 Pw[i-1] = (Pw[i-1]-(1-alfa)*Pw[i-2])/alfa
# #                 rbpts[kj] = (rbpts[kj] - beta*rbpts[kj+1])/(1.0-beta)
# #                 i += 1
# #                 j -= 1
# #                 kj -= 1
# #             end
# #             #compute knot removal error bounds (Br)
# #             if j-i <k
# #                 Br = distance4D(Pw[i-2],rbpts[kj+1])
# #             else
# #                 delta = (U[a]-Uh[i-1])/(U[b]-Uh[i-1])
# #                 A = delta*rbpts[kj+1]+(1.0-delta)*Pw[i-2]
# #                 Br = distance4D(Pw[i-1],A)
# #             end
# #             #update the error vector
# #             K = a+oldr-k
# #             q = (2*p-k+1)/2
# #             L = K-q
# #             for ii=L:a
# #                 #these knot spans were affected
# #                 e[ii] = e[ii] + Br
# #                 if e[ii] > tol
# #                     return 1 # curve not degree reducible
# #                 end
# #             end
# #             first -= 1
# #             last += 1
# #         end #for k=0:oldr-1
# #         cind = i-1
# #     end #if oldr>0
# #     #load knot vector and control points
# #     if a != p
# #         for i=0:ph-oldr-1
# #             Uh[kind] = U[a]
# #             kind += 1
# #         end
# #         for i=lbz:ph
# #             Pw[cind] = rbpts[i]
# #             cind += 1
# #         end
# #         #set up for next pass through
# #         if b<m
# #             for i=0:r-1
# #                 bpts[i] = nextbpts[i]
# #             end
# #             for i=r:p
# #                 bpts[i] = Qw[b-p+i]
# #                 a = b
# #                 b += 1
# #             end
# #         else
# #             for i=0:ph
# #                 Uh[kind+i] = U[b]
# #             end
# #         end
# #     end #while b<m
# #     nh = mh-ph-1
# #     return 0
# # end
