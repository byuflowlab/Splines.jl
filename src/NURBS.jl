"""
    NURBS(degree, knots, weights, ctrlpts)

Construct a NURBS object

# Arguments
- `deg::Integer`: degree
- `knots::Vector{Float64}`:: a knot vector (u_0, ... u_n+1)
- `weights::Vector{Float64}`:: a corresponding vector of weights
- `ctrlpts::Vector{Vector{Float64}}`:: control points.  outer index is number of control points, inner index of dimensionality of point.
"""
struct NURBS{TF, TI}
    degree::TI
    knots::Vector{TF}
    weights::Vector{TF}
    ctrlpts::Vector{Vector{TF}}
end

"""
get_degree

return degree of spline
"""
function get_degree(nurbs)
    return nurbs.degree
end

"""
get_knots

return knot vector of spline
"""
function get_knots(nurbs)
    return nurbs.knots
end

"""
get_weights

return weights of spline
"""
function get_weights(nurbs)
    return nurbs.weights
end

"""
get_ctrlpts

return control points of spline
"""
function get_ctrlpts(nurbs)
    return nurbs.ctrlpts
end

"""
Evaluate point on rational b-spline curve (NURBS, A4.1)

# Arguments
- `nurbs::NURBS: NURBS object`
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
function curvederivatives(nurbs::NURBS, u, d)

    if u > nurbs.knots[end] || u < nurbs.knots[1]
        error("parametric point, u, is outside of knot range, U")
    end

    # evaluate derivatives of A and w
    w = nurbs.weights
    Pw = [[p*w[i]; w[i]] for (i, p) in enumerate(nurbs.ctrlpts)]
    bsp = BSpline(nurbs.degree, nurbs.knots, Pw)
    CK = curvederivatives(bsp, u, d)
    Aders = [c[1:end-1] for c in CK]
    wders = [c[end] for c in CK]

    # initialize
    du = min(d, nurbs.degree)
    CK = OffsetArray(Vector{Vector{Float64}}(undef, du+1), 0:du)

    # algorithm A 4.2
    for k = 0:d
        v = Aders[k]
        for i = 1:k
            v -= binomial(k, i)*wders[i]*CK[k-i]
        end
        CK[k] = v / wders[0]
    end

    return CK  #onebased(CK)
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
- u : the knot to be added
- k : the span index at which the knot is to be inserted.
- s : numer of instances of the new knot alrady present in the knot vector, UP
- r : number of times the new knot is inserted (it is assumed that `` r+s \\leq p `` )

Outputs:
- nq : the number of control points minus 1 (the index of the last control point) after insertion
- UQ : the knot vector after insertion
- Qw : the set of weighted control points and weights after insertion
"""
function curveknotinsertion(nurbs, u, r)

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

    return nq, onebased(UQ), onebased(Qw)
end

# """
#     refineknotvectorcurve(n, p, U, Pw, X, r)

# Refine curve knot vector using NURBS A5.4.

# This algorithm is simply a knot insertion algorithm that allows for multiple knots to be added simulataneously, i.e., a knot refinement procedure.

# Inputs:
# - n : the number of control points minus 1 (the index of the last control point) before insertion
# - p : the curve order
# - U : the knot vector before insertion
# - Pw : the set of weighted control points and weights before insertion
# - X : elements, in ascending order, to be inserted into U (elements should be repeated according to their multiplicities, e.g., if x and y have multiplicites 2 and 3, X = [x,x,y,y,y])
# - r : length of X vector - 1

# Outputs:
# - Ubar : the knot vector after insertion
# - Qw : the set of weighted control points and weights after insertion
# """
# function refineknotvectorcurve(n, p, U, Pw, X, r)

#     m = n+p+1
#     a = getspanindex(n,p,X[0+1],U)
#     b = getspanindex(n,p,X[r+1],U)
#     b += 1
#     Qw = zeros(length(Pw[:,1])+r+1,length(Pw[1,:]))
#     Ubar = zeros(length(U)+r+1)
#     for j=0:a-p
#         Qw[j+1,:] = Pw[j+1,:]
#     end

#     for j=b-1:n
#         Qw[j+r+1+1,:] = Pw[j+1,:]
#     end

#     for j=0:a
#         Ubar[j+1] = U[j+1]
#     end

#     for j=b+p:m
#         Ubar[j+r+1+1] = U[j+1]
#     end
#     i = b+p-1
#     k = b+p+r
#     # if r < 0
#         for j=r:-1:0
#             while X[j+1]<=U[i+1] && i>a
#                 Qw[k-p-1+1,:] = Pw[i-p-1+1,:]
#                 Ubar[k+1] = U[i+1]
#                 k -=1
#                 i -=1
#             end
#             Qw[k-p-1+1,:] = Qw[k-p+1,:]
#             for ell=1:p
#                 ind = k-p+ell
#                 alpha = Ubar[k+ell+1] - X[j+1]
#                 if abs(alpha)==0.0
#                     Qw[ind-1+1,:] = Qw[ind+1,:]
#                 else
#                     alpha /= Ubar[k+ell+1] - U[i-p+ell+1]
#                     Qw[ind-1+1,:] = alpha*Qw[ind-1+1,:] + (1.0-alpha)*Qw[ind+1,:]
#                 end
#             end
#             Ubar[k+1] = X[j+1]
#             k -= 1
#         end
#     # end

#     return Ubar, Qw
# end

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


# """
#     degreeelevatecurve(n,p,U,Pw,t)

# Raise degree of spline from p to p `` +t``, `` t \\geq 1 `` by computing the new control point vector and knot vector.

# Knots are inserted to divide the spline into equivalent Bezier Curves. These curves are then degree elevated using the following equation.

# ```math
# \\mathbf{P}^t_i = \\sum^{\\textrm{min}(p,i)}_{j=\\textrm{max}(0,i-t)} \\frac{\\binom{p}{j} \\binom{t}{i-j} \\mathbf{P}_j}{\\binom{p+t}{i}}~,~~~~~i=0,...,p+t
# ```
# where `` \\mathbf{P}^t_i `` are the degree elevated control points after `` t `` -degree elevations

# Finally, the excess knots are removed and the degree elevated spline is returned.

# (see NURBS eqn 5.36, A5.9)

# Inputs:
# - n : the number of control points minus 1 (the index of the last control point) before degree elevation
# - p : the curve order
# - U : the knot vector before degree elevation
# - Pw : the set of weighted control points and weights before degree elevation
# - t : the number of degrees to elevate, i.e. the new curve degree is p+t

# Outputs:
# - nh : the number of control points minus 1 (the index of the last control point) after degree elevation
# - Uh : the knot vector after degree elevation
# - Qw : the set of weighted control points and weights after degree elevation
# """
# function degreeelevatecurve(n,p,U,Pw,t)
#     m = n+p+1
#     ph = p+t
#     ph2 = Int(floor(ph/2))

#     #initialize in outer scope:
#     oldr = 0
#     ub = 0

#     #initialize local arrays
#     bezalfs = zeros(p+t+1,p+1) #coefficients for degree elevating the bezier segments
#     bpts = zeros(p+1,length(Pw[1,:])) #pth-degree bezier control points of the current section
#     ebpts = zeros(p+t+1,length(Pw[1,:])) #(p+t)th-degree bezier control points of the current section
#     nextbpts = zeros(p-1,length(Pw[1,:])) #leftmost control points of next bezier segement
#     alfs = zeros(p-1) #knot insertion alphas.

#     #initialize outputs:
#     s = length(unique(U)) - 2 #number of unique internal knots
#     Uh = ones( length(U) + t*(s+2) ) #mhat eqn 5.33 is m + (s+2)xt
#     Qw = ones( length(Pw[:,1]) + t*(s+1), length(Pw[1,:]) ) #nhat eqn 5.32 is n+(s+1)xt

#     #compute bezier degree elevation coefficients
#     bezalfs[ph+1,p+1] = 1.0 #bezalfs are coefficients for degree elevating the Bezier segments
#     bezalfs[0+1,0+1] = 1.0
#     for i=1:ph2
#         inv = 1.0/binomialcoeff(ph,i)
#         for j=Int(max(0,i-t)):Int(min(p,i))
#             bezalfs[i+1,j+1] = inv*binomialcoeff(p,j)*binomialcoeff(t,i-j) #from eqn 5.36
#         end
#     end
#     for i = ph2+1:ph-1
#         for j=Int(max(0,i-t)):Int(min(p,i))
#             bezalfs[i+1,j+1] = bezalfs[ph-i+1,p-j+1]
#         end
#     end
#     mh = ph
#     kind = ph+1
#     r = -1
#     a = p
#     b = p+1
#     cind = 1
#     ua = U[0+1]

#     #first new control point is first old control point
#     Qw[0+1,:] = Pw[0+1,:]

#     #fill in first p+t new U vector points (same as old vector with +t multiplicity)
#     for i=0:ph
#         Uh[i+1] = ua
#     end

#     #initialize first bezier segment
#     for i=0:p
#         bpts[i+1,:] = Pw[i+1,:] #bpts are the pth-degree bezier control points of the current segment
#     end

#     #big loop through knot vector starting after first elements of Qw and Uh vectors.
#     while b<m
#         i = b
#         while b<m && U[b+1] == U[b+1+1] #incrememnt b for knots with multiplicities.
#             b += 1
#         end
#         mul = b-i+1
#         mh += mul + t
#         ub = U[b+1]
#         oldr = r
#         r = p-mul

#         #insert knot u(b) r times
#         if oldr > 0
#             lbz = Int(floor((oldr + 2)/2))
#         else
#             lbz = 1
#         end
#         if r>0
#             rbz = Int(ph-floor((r+1)/2))
#         else
#             rbz = ph
#         end

#         #insert knot to get bezier segment
#         if r>0
#             numer = ub - ua
#             for k=p:-1:mul+1
#                 alfs[k-mul-1+1] = numer/(U[a+k+1]-ua) #knot insertion alphas
#             end
#             for j=1:r
#                 save = r-j
#                 s = mul+j
#                 for k=p:-1:s
#                     bpts[k+1,:] = alfs[k-s+1]*bpts[k+1,:] + (1.0-alfs[k-s+1])*bpts[k-1+1,:]
#                 end #for
#                 nextbpts[save+1,:] = bpts[p+1,:] #nextbpts are the leftmost control points of the next bezier segment
#             end #for
#         end #insert knot if

#         #degree elevate bezier
#         for i=lbz:ph
#             #only points lbz,...,ph are used below
#             ebpts[i+1,:] .= 0.0 #ebpts are the (p+t)th-degree bezier control points of the current segment.
#             for j = Int(max(0,i-t)):Int(min(p,i))
#                 ebpts[i+1,:] += bezalfs[i+1,j+1]*bpts[j+1,:]
#             end
#         end #for; degree elevation

#         #must remove knot u=U[a] oldr times
#         if oldr > 1
#             first = kind-2
#             last = kind
#             den = ub-ua
#             bet = (ub-Uh[kind-1+1])/den

#             #knot removal loop
#             if oldr-1 >= 1 #need to make sure this is true since julia does a loop no matter what, when C evaluates condition first.
#                 for tr=1:oldr-1
#                     i = first
#                     j = last
#                     kj = j-kind+1

#                     #loop and compute the new control points for one removal step
#                     while j-i > tr
#                         if i < cind
#                             alf = (ub - Uh[i+1])/(ua-Uh[i+1])
#                             Qw[i+1,:] = alf*Qw[i+1,:] + (1.0-alf)*Qw[i-1+1,:]
#                         end #if
#                         if j >= lbz
#                             if j-tr <= kind-ph+oldr
#                                 gam = (ub-Uh[j-tr+1])/den
#                                 ebpts[kj+1,:] = gam*ebpts[kj+1,:]+(1.0-gam)*ebpts[kj+1+1,:]
#                             else
#                                 ebpts[kj+1,:] = bet*ebpts[kj+1,:]+(1.0-bet)*ebpts[kj+1+1,:]
#                             end #if
#                         end #if
#                         i += 1
#                         j -= 1
#                         kj -= 1
#                     end #while
#                     first -= 1
#                     last += 1
#                 end #for tr
#             end #if oldr is big enough
#         end #if; remove knot u=U[a]

#         #load the knot ua
#         if a != p
#             for i=0:ph-oldr-1
#                 Uh[kind+1] = ua
#                 kind += 1
#             end #for
#         end #if

#         #load ctrl pts into Qw
#         for j=lbz:rbz
#             Qw[cind+1,:] = ebpts[j+1,:]
#             cind += 1
#         end

#         #set up for next pass through loop
#         if b < m
#             for j=0:r-1
#                 bpts[j+1,:] = nextbpts[j+1,:]
#             end
#             for j=r:p
#                 bpts[j+1,:] = Pw[b-p+j+1,:]
#             end
#             a = b
#             b += 1
#             ua = ub
#         else #end knot
#             for i=0:ph
#                 Uh[kind+i+1] = ub
#             end
#         end
#     end #while b<m
#     nh = mh-ph-1

# return nh, Uh, Qw
# end


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

# """
#     nurbsbasis(u,p,d,U,w)

# Get rational basis functions and derivatives.

# ```math
# R_{i,p}(u) = \\frac{N_{i,p}(u)w_i}{\\sum_{j=0}^n N_{j,p}(u)w_j}
# ```

# where `` N_{i,p}(u ) `` are B-Spline Basis Functions and `` w_i `` are weights associated with the NURBS control points.

# (see NURBS eqn 4.2)

# Inputs:

# - u : parametric point of interest
# - p : the curve order
# - d : the max derivative order (n ≦ p)
# - U : the knot vector
# - w : control point weights

# """
# function nurbsbasis(u,p,d,U,w)
#     #get the span index
#     n = length(U)-1-p
#     i = getspanindex(n, p, u, U)

#     #get B-Spline basis functions and derivatives for numerator
#     bases = Splines.basisfunctionsderivatives(i+1,u,p,d,U)
#     #separate out the bases and their derivatives
#     N = bases[1,:]
#     dN = bases[2,:] #only doing first derivates right now

#     #number of non-zero basis functions at parametric point, u.
#     numbasisfunctions = length(N)
#     Nw = 0
#     dNw = 0

#     #calculate the denomenator values
#     for j=1:numbasisfunctions
#         Nw += N[j] .* w[j]
#         dNw += dN[j] .* w[j]
#     end

#     #Calculate each of the non-zero Rational Basis functions and first derivatives
#     R = zeros(length(N))
#     dR = zeros(size(N))
#     for j=1:numbasisfunctions
#         R[j] = N[j]/Nw .* w[j]
#         dR[j] = w[j] .* (Nw*dN[j] - dNw*N[j]) ./ Nw^2
#     end

#     return R, dR
# end

