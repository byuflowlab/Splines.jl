module BSpline
"""
    getSpanIndex(n::Int64,p::Int64,u,U)

Complete binary search to find span index of vector, U, in which knot, u, lies.
"""
function getSpanIndex(p::Int64,u,U)
    n = length(U)-p-1
    if u == U[n+1]
        return n #special case
    else
        lo = p+1
        hi = n+2
        mid = div((lo+hi),2)
        while u<U[mid] || u >= U[mid+1]
            if u<U[mid]
                hi = mid
            else
                lo = mid
            end #if low
            mid = div((lo+hi),2)
        end #while not found
        return mid+1
    end #if not end
end

"""
    basisFunctions(i,u,p,U)

Calculate the non-vanishing basis function of the B-Spline of order p, defined by knots U at knot u.
"""
function basisFunctions(i,u,p,U)
    N = ones(p+1)
    left = zeros(p+1)
    right = zeros(p+1)
    for j=1:p
        left[j] = u-U[i+1-j]
        right[j] = U[i+j]-u
        saved = 0.0
        for r=0:j-1
            temp = N[r+1]/(right[r+1] + left[j-r])
            N[r+1] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        end
        N[j+1] = saved
    end
    return N
end

"""
    basisFunctionsDerivatives(i,u,p,U)

Calculate the non-vanishing basis function of the B-Spline of order p, defined by knots U at knot u.
"""
function basisFunctionsDerivatives(i,u,p,U)
    #Initialize
    n = length(U)-p-1
    ndu = ones(p+1,p+1)
    a = zeros(p+1,p+2)
    ders = zeros(p+1,p+1)
    left = zeros(p+1)
    right = zeros(p+1)

    #---Compute (and save) Basis Functions and Knot Differences
    for j=1:p
    left[j] = u-U[i+1-j]
    right[j] = U[i+j]-u
    saved = 0.0
        for r=0:j-1
            #upper triangle (basis functions)
            ndu[j+1,r+1] = right[r+1] + left[j-r]
            temp = ndu[r+1,j]/ndu[j+1,r+1]
            #lower triangle (knot differences)
            ndu[r+1,j+1] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        end #for r
    ndu[j+1,j+1] = saved
    end #for j

    #Load Basis Functions
    for j=1:p+1
        ders[1,j] = ndu[j,p+1]
    end #for j

    #---Compute Derivatives
    # println("\n\nDerivative Algorithm:\nFor knot span ", i)
    for r=0:p
        # println("For basis fuction: ", r)
        #Set row indices for coefficient array (swaps back and forth as they're computed/used)
        s1 = 0
        s2 = 1
        a[0+1,0+1] = 1.0
        for k=1:p
            # println("Computing Derivative ", k)
            # println("Coefficient Array Row Indices: s1 = ",s1,"\ts2 = ",s2)
            d = 0.0
            rk = r-k
            pk = p-k
            # println("rk = ", rk)
            # println("pk = ", pk)
            if r >= k
                # println("r >= k")
                a[s2+1,0+1] = a[s1+1,0+1]/ndu[pk+1+1,rk+1]
                d = a[s2+1,0+1]*ndu[rk+1,pk+1]
            end #if r>=k

            if rk >= -1
                # println("rk >= -1,\t j1 = 1")
                j1 = 1
            else
                # println("rk < -1,\t j1 = -rk")
                j1 = -rk
            end #if rk>=-1

            if r-1 <= pk
                # println("r-1 <= pk,\t j2 = k-1")
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
                    # println("ndu[pk+1,rk+j] = ", ndu[pk+1+1,rk+j+1])
                    # println("a[s1,j] = ", a[s1+1,j+1])
                    # println("a[s1,j-1] = ", a[s1+1,j-1+1])
                    a[s2+1,j+1] = (a[s1+1,j+1]-a[s1+1,j-1+1])/ndu[pk+1+1,rk+j+1]
                    # println("a[s2,j] = ", a[s2+1,j+1])
                    # println("ndu[rk+j,pk+1] = ", ndu[rk+j+1,pk+1])
                    d += a[s2+1,j+1]*ndu[rk+j+1,pk+1]
                    # println("d = ", d)
                end #for j
            end

            if r <= pk
                # println("r <= pk")
                # println("a[s1,k-1] = ", a[s1+1,k-1+1])
                # println("ndu[pk+1,r] = ", ndu[pk+1+1,r+1])
                a[s2+1,k+1] = -a[s1+1,k-1+1]/ndu[pk+1+1,r+1]
                # println("a[s2,k-1] = ", a[s2+1,k+1])
                # println("a[s2,k] = ", a[s2+1,k+1])
                # println("ndu[r,pk+1] = ", ndu[r+1,pk+1])
                d += a[s2+1,k+1]*ndu[r+1,pk+1]
                # println("d = ", d)
            end #if r
            ders[k+1,r+1] = d
            # println("derivatives: N",r+2,p,"^",k)
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
    for k=1:p
        for j=0:p
            ders[k+1,j+1] *= factorial(p)/factorial(p-k)
        end
    end

    # #---Compute derivatives the long way
    # display(ndu)
    # println()
    # ders[1,:] = ndu[:,p+1]
    # for r = 2:p+1
    #     println("r = ", r)
    #     for k = 1:p
    #         println("k = ", k)
    #         a[1,1] = 1.0
    #         for j=1:k
    #             println("j = ", j)
    #             println("a[$k,$j]")
    #             if j==1 && k>1
    #                 println("j==1 && k!=1")
    #                 a[k,j] = a[k-1,1]/(U[i+p-k+1] - U[i])
    #             elseif k>1 && j>1 && j!=k
    #                 println("k>1 && j>1 && j!=k")
    #                 a[k,j] = (a[k-1,j] - a[k-1,j-1])/(U[i+p+j-k+1] - U[i+j])
    #             elseif k>1 && j==k
    #                 println("k!=1 && j==k")
    #                 a[k,j] = -a[k-1,k-1]/(U[i+p+1] - U[i+k])
    #             end
    #             if a[k,j] == Inf || a[k,j] == -Inf || isnan(a[k,j])
    #                 println("a[k,j] == inf or nan")
    #                 a[k,j] = 0.0
    #             end
    #             println("a[k,j] = ", a[k,j])
    #             if i+j > i+1
    #                 N = 0
    #             else
    #                 N = ndu[i+j,p-k]
    #             end
    #             ders[r,k] += a[k,j]*N
    #         end
    #         ders[r,k] *= factorial(p)/factorial(p-k)
    #     end
    # end


    return ders

end #function

end #module BSpline