"""
    curvepoint(n, p, U, Pw, u)

Compute point on rational B-Spline curve defined as:

```math
\\mathbf{C}^w(u) = \\sum_{i=0}^n N_{i, p}(u) \\mathbf{P}_i^w
```

where `` \\mathbf{P}_i^w `` are the set of weighted control points and weights such that `` \\mathbf{P}_i^w = (w_ix_i, w_iy_i, w_iz_i, w_i)``.

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
    if u>U[end] || u<U[1]
        error("Nurbs.curvepoint: parametric point, u, is outside of knot range, U.")
    end
    span = getspanindex(n, p, u, U)
    # println("span = ", span)
    N = basisfunctions(span+1, u, p, U)
    Cw = zeros(length(Pw[1,:]))
    for j = 0:p
        # println("N[$j] = ", N[j+1])
        # println("Pw[$(span-p+j)] = ", Pw[span-p+j+1, :])
        # println("span = ", span)
        # println("p = ", p)
        # println("j = ", j)
        # println("j+1 = ", j+1)
        # println("span-p+j+1 = ", span-p+j+1)
        Cw .+= N[j+1]*Pw[span-p+j+1, :]
    end
    # println(Cw)
    return Cw/Cw[end]
end

"""
    rationalcurvederivatives(Aders, wders, d)

Compute the point `` \\mathbf{C}(u) `` and the derivatives `` \\mathbf{C}^{(k)}(u) `` for `` 1 \\leq k \\leq d `` where:

```math
\\mathbf{C}^{(k)}(u) = \\frac{ \\mathbf{A}^{(k)}(u) - \\sum_{i=1}^k \\binom{k}{i} w^{(i)}(u) \\mathbf{C}^{(k-1)}(u) }{w(u)}
```

where `` \\mathbf{A}^{(k)}(u) `` and `` w^{(i)}(u) `` are precomputed using preweighted control points for some parametric point, `` 0 \\leq u \\leq 1``, from ` ```curvederivatives1``` and are inputs Aders and wders, respectively.

(see NURBS eqn 4.8 and A4.2)
"""
function rationalcurvederivatives(Aders, wders, d)
    CK = zeros(d+1, length(Aders[1, :]))
    for k=0:d
        v = Aders[k+1, :]
        # println("Aders[ `` k] = ", Aders[k+1, :])
        for i=1:k
            if k >= i
                # println("binom = ", binomialcoeff(k, i))
                # println("wders[ `` i] = ", wders[i+1])
                # println("Ck[ `` (k-i)] = ", CK[k-i+1, :])
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
- r : number of times the new knot is inserted (it is assumed that `` r+s \\leq p `` )

Outputs:
- nq : the number of control points minus 1 (the index of the last control point) after insertion
- UQ : the knot vector after insertion
- Qw : the set of weighted control points and weights after insertion
"""
function curveknotinsertion(np, p, UP, Pw, u, k, s, r)
    mp = np+p+1
    nq = np+r

    #initialize output vectors
    UQ = zeros(length(UP)+r)
    # println("UP, `` UP")
    # println("UQ, `` UQ")
    Qw = zeros(nq+1, length(Pw[1, :]))
    Rw = zeros(p+1, length(Pw[1, :]))
    #Load new knot vector
    for i=0:k
        UQ[i+1] = UP[i+1]
        # println("UQ `` i = ", UQ[i+1])
    end
    # println("UQ, `` UQ")
    for i=1:r
        UQ[k+i+1] = u
        # println("UQ `` (k+i) = ", UQ[k+i+1])
    end
    # println("UQ, `` UQ")
    for i=k+1:mp
        UQ[i+1+r] = UP[i+1]
        # println("UQ `` (i+r) = ", UQ[i+1+r])
    end
    # println("UQ, `` UQ")
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
        # println("Qw `` L = ", Qw[L+1, :])
        Qw[k+r-j-s+1, :] = Rw[p-j-s+1, :]
        # println("Qw `` (k+r-j-s) = ", Qw[k+r-j-s+1, :])
    end
    #Load remaining control points
    for i=L+1:k-s
        Qw[i+1, :] = Rw[i+1-L, :]
    end


    return nq, UQ, Qw
end

"""
    refineknotvectorcurve(n, p, U, Pw, X, r)

Refine curve knot vector using NURBS A5.4.

This algorithm is simply a knot insertion algorithm that allows for multiple knots to be added simulataneously, i.e., a knot refinement procedure.

Inputs:
- n : the number of control points minus 1 (the index of the last control point) before insertion
- p : the curve order
- U : the knot vector before insertion
- Pw : the set of weighted control points and weights before insertion
- X : elements, in ascending order, to be inserted into U (elements should be repeated according to their multiplicities, e.g., if x and y have multiplicites 2 and 3, X = [x,x,y,y,y])
- r : length of X vector

Outputs:
- Ubar : the knot vector after insertion
- Qw : the set of weighted control points and weights after insertion
"""
function refineknotvectorcurve(n, p, U, Pw, X, r)
    m = n+p+1
    a = getspanindex(n,p,X[0+1],U)
    # println("a = ",a)
    b = getspanindex(n,p,X[r+1],U)
    b += 1
    # println("b = ",b)
    Qw = zeros(length(Pw[:,1])+r+1,length(Pw[1,:]))
    Ubar = zeros(length(U)+r+1)
    for j=0:a-p
        Qw[j+1,:] = Pw[j+1,:]
    end
    for j=b-1:n
        # println(j+r+1+1)
        Qw[j+r+1+1,:] = Pw[j+1,:]
    end
    for j=0:a
        Ubar[j+1] = U[j+1]
        # println("Ubar[ `` j] = ", U[j+1])
    end
    for j=b+p:m
        Ubar[j+r+1+1] = U[j+1]
        # println("Ubar[ `` (j+r+1)] = ", U[j+1])
    end
    i = b+p-1
    k = b+p+r
    # if r < 0
        for j=r:-1:0
            # println("r = ",r)
            while X[j+1]<=U[i+1] && i>a
                Qw[k-p-1+1,:] = Pw[i-p-1+1,:]
                Ubar[k+1] = U[i+1]
                # println("Ubar[ `` k] = ", U[i+1])
                k -=1
                i -=1
            end
            Qw[k-p-1+1,:] = Qw[k-p+1,:]
            for ell=1:p
                ind = k-p+ell
                # println("ind = ", ind)
                alpha = Ubar[k+ell+1] - X[j+1]
                if abs(alpha)==0.0
                    Qw[ind-1+1,:] = Qw[ind+1,:]
                else
                    alpha /= Ubar[k+ell+1] - U[i-p+ell+1]
                    Qw[ind-1+1,:] = alpha*Qw[ind-1+1,:] + (1.0-alpha)*Qw[ind+1,:]
                end
            end
            Ubar[k+1] = X[j+1]
            # println("Ubar[ `` k] = ", X[j+1])
            k -= 1
        end
    # end

    return Ubar, Qw
end


# """
# """
# function removecurveknot(n,p,U,Pw,u,r,s,num)
#     m = n+p+1
#     ord = p+1
#     fout = (2*r-s-p)/2 #first control point output
#     last = r-s
#     first = r-p
#     for t = 0:num
#         #this loop is eqn 5.28
#         off = first-1 #diff in index between temp and P
#         temp[0] = P2[off]
#         temp[last+1-off] = Pw[last+1]
#         i = first
#         j = last
#         ii = 1
#         jj = last-off
#         remflag = 0
#         while j-i>t
#             #compute new control points for one removal step
#             alfi = (u-U[i])/(U[i+ord+t]-U[i])
#             alfj = (u-U[j-t])/(U[j+ord]-U[j-t])
#             temp[ii] = (Pw[i] - (1.0-alfi)*temp[ii-1])/alfi
#             temp[jj] = (Pw[j]-alfj*temp[jj+1])/(1.0-alfj)
#             i += 1
#             ii += 1
#             j -= 1
#             jj -= 1
#         end #while
#         if j-i <t #check if knot is removable
#             if distance4D(temp[ii-1],temp[jj-1]) < tol
#                 remflag = 1
#             end
#         else
#             alfi = (u-U[i])/(U[i+ord+t]-U[i])
#             if distance4D(Pw[i],alfi*temp[ii+t+1]+(1.0-alfi)*temp[ii-1]) <= tol
#                 remflag = 1
#             end
#         end
#         if remflag == 0 #cannot remove any more knots
#             break #get out of for loop
#         else
#             #successful removal. save new cont pts
#             i = first
#             j = last
#             while j-1 > t
#                 Pw[i] = temp[i-off]
#                 Pw[j] = temp[j-off]
#                 i += 1
#                 j -= 1
#             end #while
#         end #if
#         first -= 1
#         last += 1
#     end #for
#     if t==0
#         return
#     end
#     for k = r+1:m
#         U[k-t] = U[k] #shift knots
#         j = fout
#         i = j #Pj through Pi will be overwritten
#         for k=1:t-1
#             if mod(k,2) == 1 #k modulo 2
#                 i += 1
#             else
#                 j -= 1
#             end
#         end
#         for k=i+1:n #shift
#             Pw[j] = Pw[k]
#             j += 1
#         end
#         return
#     end
# end


"""
    degreeelevatecurve(n,p,U,Pw,t)

Raise degree of spline from p to p `` +t``, `` t \\geq 1 `` by computing the new control point vector and knot vector.

Knots are inserted to divide the spline into equivalent Bezier Curves. These curves are then degree elevated using the following equation.

```math
\\mathbf{P}^t_i = \\sum^{\\textrm{min}(p,i)}_{j=\\textrm{max}(0,i-t)} \\frac{\\binom{p}{j} \\binom{t}{i-j} \\mathbf{P}_j}{\\binom{p+t}{i}}~,~~~~~i=0,...,p+t
```
where `` \\mathbf{P}^t_i `` are the degree elevated control points after `` t `` -degree elevations

Finally, the excess knots are removed and the degree elevated spline is returned.

(see NURBS eqn 5.36, A5.9)

Inputs:
- n : the number of control points minus 1 (the index of the last control point) before degree elevation
- p : the curve order
- U : the knot vector before degree elevation
- Pw : the set of weighted control points and weights before degree elevation
- t : the number of degrees to elevate, i.e. the new curve degree is p+t

Outputs:
- nh : the number of control points minus 1 (the index of the last control point) after degree elevation
- Uh : the knot vector after degree elevation
- Qw : the set of weighted control points and weights after degree elevation
"""
function degreeelevatecurve(n,p,U,Pw,t)
    m = n+p+1
    ph = p+t
    ph2 = Int(floor(ph/2))

    #initialize in outer scope:
    oldr = 0
    ub = 0

    #initialize local arrays
    bezalfs = zeros(p+t+1,p+1) #coefficients for degree elevating the bezier segments
    bpts = zeros(p+1,length(Pw[1,:])) #pth-degree bezier control points of the current section
    ebpts = zeros(p+t+1,length(Pw[1,:])) #(p+t)th-degree bezier control points of the current section
    nextbpts = zeros(p-1,length(Pw[1,:])) #leftmost control points of next bezier segement
    alfs = zeros(p-1) #knot insertion alphas.

    #initialize outputs:
    s = length(unique(U)) - 2 #number of unique internal knots
    Uh = ones( length(U) + t*(s+2) ) #mhat eqn 5.33 is m + (s+2)xt
    Qw = ones( length(Pw[:,1]) + t*(s+1), length(Pw[1,:]) ) #nhat eqn 5.32 is n+(s+1)xt

    #compute bezier degree elevation coefficients
    bezalfs[ph+1,p+1] = 1.0 #bezalfs are coefficients for degree elevating the Bezier segments
    bezalfs[0+1,0+1] = 1.0
    for i=1:ph2
        inv = 1.0/binomialcoeff(ph,i)
        for j=Int(max(0,i-t)):Int(min(p,i))
            bezalfs[i+1,j+1] = inv*binomialcoeff(p,j)*binomialcoeff(t,i-j) #from eqn 5.36
        end
    end
    for i = ph2+1:ph-1
        for j=Int(max(0,i-t)):Int(min(p,i))
            bezalfs[i+1,j+1] = bezalfs[ph-i+1,p-j+1]
        end
    end
    mh = ph
    kind = ph+1
    r = -1
    a = p
    b = p+1
    cind = 1
    ua = U[0+1]

    #first new control point is first old control point
    Qw[0+1,:] = Pw[0+1,:]

    #fill in first p+t new U vector points (same as old vector with +t multiplicity)
    for i=0:ph
        Uh[i+1] = ua
    end

    #initialize first bezier segment
    for i=0:p
        bpts[i+1,:] = Pw[i+1,:] #bpts are the pth-degree bezier control points of the current segment
    end

    #big loop through knot vector starting after first elements of Qw and Uh vectors.
    while b<m
        i = b
        while b<m && U[b+1] == U[b+1+1] #incrememnt b for knots with multiplicities.
            b += 1
        end
        mul = b-i+1
        mh += mul + t
        ub = U[b+1]
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
                alfs[k-mul-1+1] = numer/(U[a+k+1]-ua) #knot insertion alphas
            end
            for j=1:r
                save = r-j
                s = mul+j
                for k=p:-1:s
                    bpts[k+1,:] = alfs[k-s+1]*bpts[k+1,:] + (1.0-alfs[k-s+1])*bpts[k-1+1,:]
                end #for
                nextbpts[save+1,:] = bpts[p+1,:] #nextbpts are the leftmost control points of the next bezier segment
            end #for
        end #insert knot if

        #degree elevate bezier
        for i=lbz:ph
            #only points lbz,...,ph are used below
            ebpts[i+1,:] .= 0.0 #ebpts are the (p+t)th-degree bezier control points of the current segment.
            for j = Int(max(0,i-t)):Int(min(p,i))
                ebpts[i+1,:] += bezalfs[i+1,j+1]*bpts[j+1,:]
            end
        end #for; degree elevation

        #must remove knot u=U[a] oldr times
        if oldr > 1
            first = kind-2
            last = kind
            den = ub-ua
            bet = (ub-Uh[kind-1+1])/den

            #knot removal loop
            if oldr-1 >= 1 #need to make sure this is true since julia does a loop no matter what, when C evaluates condition first.
                for tr=1:oldr-1
                    i = first
                    j = last
                    kj = j-kind+1

                    #loop and compute the new control points for one removal step
                    while j-i > tr
                        if i < cind
                            alf = (ub - Uh[i+1])/(ua-Uh[i+1])
                            Qw[i+1,:] = alf*Qw[i+1,:] + (1.0-alf)*Qw[i-1+1,:]
                        end #if
                        if j >= lbz
                            if j-tr <= kind-ph+oldr
                                gam = (ub-Uh[j-tr+1])/den
                                ebpts[kj+1,:] = gam*ebpts[kj+1,:]+(1.0-gam)*ebpts[kj+1+1,:]
                            else
                                ebpts[kj+1,:] = bet*ebpts[kj+1,:]+(1.0-bet)*ebpts[kj+1+1,:]
                            end #if
                        end #if
                        i += 1
                        j -= 1
                        kj -= 1
                    end #while
                    first -= 1
                    last += 1
                end #for tr
            end #if oldr is big enough
        end #if; remove knot u=U[a]

        #load the knot ua
        if a != p
            for i=0:ph-oldr-1
                Uh[kind+1] = ua
                kind += 1
            end #for
        end #if

        #load ctrl pts into Qw
        for j=lbz:rbz
            Qw[cind+1,:] = ebpts[j+1,:]
            cind += 1
        end

        #set up for next pass through loop
        if b < m
            for j=0:r-1
                bpts[j+1,:] = nextbpts[j+1,:]
            end
            for j=r:p
                bpts[j+1,:] = Pw[b-p+j+1,:]
            end
            a = b
            b += 1
            ua = ub
        else #end knot
            for i=0:ph
                Uh[kind+i+1] = ub
            end
        end
    end #while b<m
    nh = mh-ph-1

return nh, Uh, Qw
end


# """
# """
# function degreereducecurve(n,p,U,Qw,nh,Uh,Pw)
#     ph = p-1
#     mh = ph
#     kind = pu+1
#     r = -1
#     a = p
#     b = p+1
#     ind = 1
#     mult = p
#     m = n+p+2
#     Pw[0] = Qw[0]
#     for i=0:ph #compute left end of knot vecotr
#         Uh[i] = U[0]
#     end
#     for i=0:p #initialize first bezier segment
#         bpts[i] = Qw[i]
#     end
#     for i=0:m-1 #initialize error vector
#         e[i] = 0.0
#     end
#     #loop through the knot vector
#     while b<m
#         #first compute knot multiplicity
#         i = b
#         while b<m && U[b] == U[b+1]
#         b += 1
#         end
#     mult = b-i+1
#     mh += mult-1
#     oldr = r
#     r = p-mult
#     if older > 0
#         lbz = (oldr+2)/2
#     else
#         lbz = 1
#     end
#     #insert knot U[b] r times
#     if r>0
#         numer = U[b] - U[a]
#         for k=p:-1:mult
#             alphas[k-mult-1] = numer/(U[a+k]-U[a])
#         end
#         for j=1:r
#             save = r-j
#             s= mult+j
#             for k=p:-1:s
#                 bpts[k] = alphas[k-s]*bpts[k] + (1.0-alphas[k-s])*bpts[k-1]
#             end
#             nextbpts[save] = bpts[p]
#         end
#     end
#     #degree reduce bezier segement
#     BezDegreeReduce(bpts,rbpts,MaxErr)
#     e[a] = e[a]+MaxErr
#     if e[a] > tol
#         return 1 #curve not degree reducible
#         #remove knot U[a] oldr times
#     end
#     if oldr > 0
#         first = kind
#         last = kind
#         for k=0:oldr-1
#             i = first
#             j = last
#             kj = j-kind
#             while j-i > k
#                 alfa = (U[a]-Uh[i-1])/(U[b] - Uh[j-k-1])
#                 Pw[i-1] = (Pw[i-1]-(1-alfa)*Pw[i-2])/alfa
#                 rbpts[kj] = (rbpts[kj] - beta*rbpts[kj+1])/(1.0-beta)
#                 i += 1
#                 j -= 1
#                 kj -= 1
#             end
#             #compute knot removal error bounds (Br)
#             if j-i <k
#                 Br = distance4D(Pw[i-2],rbpts[kj+1])
#             else
#                 delta = (U[a]-Uh[i-1])/(U[b]-Uh[i-1])
#                 A = delta*rbpts[kj+1]+(1.0-delta)*Pw[i-2]
#                 Br = distance4D(Pw[i-1],A)
#             end
#             #update the error vector
#             K = a+oldr-k
#             q = (2*p-k+1)/2
#             L = K-q
#             for ii=L:a
#                 #these knot spans were affected
#                 e[ii] = e[ii] + Br
#                 if e[ii] > tol
#                     return 1 # curve not degree reducible
#                 end
#             end
#             first -= 1
#             last += 1
#         end #for k=0:oldr-1
#         cind = i-1
#     end #if oldr>0
#     #load knot vector and control points
#     if a != p
#         for i=0:ph-oldr-1
#             Uh[kind] = U[a]
#             kind += 1
#         end
#         for i=lbz:ph
#             Pw[cind] = rbpts[i]
#             cind += 1
#         end
#         #set up for next pass through
#         if b<m
#             for i=0:r-1
#                 bpts[i] = nextbpts[i]
#             end
#             for i=r:p
#                 bpts[i] = Qw[b-p+i]
#                 a = b
#                 b += 1
#             end
#         else
#             for i=0:ph
#                 Uh[kind+i] = U[b]
#             end
#         end
#     end #while b<m
#     nh = mh-ph-1
#     return 0
# end

"""
    nurbsbasis(i,p,u,U,w)

Get rational basis functions and derivatives. see eqn 4.2

```math
R_{i,p}(u) = \\frac{N_{i,p}(u)w_i}{\\sum_{j=0}^n N_{j,p}(u)w_j}
```

where `` N_{i,p}(u ) `` are B-Spline Basis Functions and `` w_i `` are weights associated with the NURBS control points.

Inputs:

- u : parametric point of interest
- p : the curve order
- n : the max derivative order (n ≦ p)
- U : the knot vector
- weights : control point weights

"""
function nurbsbasis(u,p,d,U,weights)
    #get the span index
    # println("gettin index")
    # println("d = ", d)
    # println("p = ", p)
    # println("u = ", u)
    # println("U = ", U)
    n = length(U)-1-p
    i = getspanindex(n, p, u, U)
    # println("i = ", i)
    #!this may not be right, make sure that these weights are right (doesn't matter if weights are always 1, but they may not be)
    w = weights[i:i+p]
    # println("weights = ", w)
    #get B-Spline basis functions and derivatives for numerator
    bases = Splines.basisfunctionsderivatives(i+1,u,p,d,U)
    # println("my bases: ", bases)
    #separate out the bases and their derivatives
    N = bases[1,:]
    # println("N: ", N)
    dN = bases[2,:] #only doing first derivates right now
    # println("dN: ", dN)
    #number of non-zero basis functions at parametric point, u.
    numbasisfunctions = length(N)
    Nw = 0
    dNw = 0

    #calculate the denomenator values
    for j=1:numbasisfunctions
        Nw += N[j] .* w[j]
        dNw += dN[j] .* w[j]
    end
    # println("Nw = ", Nw)
    # println("dNw = ", dNw)
    #Calculate each of the non-zero Rational Basis functions and first derivatives
    R = zeros(length(N))
    dR = zeros(size(N))
    for j=1:numbasisfunctions
        R[j] = N[j]/Nw .* w[j]
        dR[j] = w[j] .* (Nw*dN[j] - dNw*N[j]) ./ Nw^2
    end
    # println("R = ", R)
    # println("dR = ", dR)

    return R, dR
end