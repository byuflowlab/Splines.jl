"""
    binomialcoeff(n, i)

Calculate the Binomial Coefficient defined as:

```math
\\binom{n}{i} = \\frac{n!}{i!(n-1)!}
```
"""
function binomialcoeff(n, i)
    return factorial(n)./(factorial(i).*factorial(n-i))
end

"""
    bernsteincoeff(u, n, i)

Calculate Bernstein Coefficient (Bezier Basis Function) defined as:

```math
B_{i,n}(u) = \\binom{n}{i} u^i (1-u)^{n-1}
```

at parametric point, \$u\$, where \$ 0 \\leq u \\leq 1 \$. \$u\$ may either be a single value or an array.

(see NURBS, eqn 1.8)
"""
function bernsteincoeff(u, n, i)
    return binomialcoeff(n,i) .* u.^i .* (1-u).^(n-i)
end

"""
    simple_bezier1D(P, u)

Calculate a point along a Bezier curve at the parametric point, \$u\$, based on the control points  \$\\mathbf{P}\$ where the Bezier curve, \$\\mathbf{C}(u)\$ is defined as:

```math
\\mathbf{C}(u) = \\sum_{i=0}^n B_{i,n}(u) \\mathbf{P}_i~,~~~~ 0 \\leq u \\leq 1
```
where \$ B \$ is the basis (Bernstein Coefficient) at parametric point, \$u\$, as calculated from ```bernsteincoeff```, and n is the number of control points in vector \$\\mathbf{P}\$. Again, \$u\$ may either be a single value or an array.

(see NURBS eqn 1.7)
"""
function simple_bezier1D(P, u)
    n = length(P[:,1])
    bc = zeros(length(u),n)

    for i=1:n
        bc[:,i] = bernsteincoeff(u, n, i)
    end

    C = bc*P

    return C
end


"""
    decasteljau_bezier1D(P, u)

Calculate a point along a Bezier curve at the parametric point, \$u_0\$, based on the control points \$\\mathbf{P}\$ using deCasteljau's Algorithm:

```math
\\mathbf{P}_{k,i}(u_0) = (1- u_0)\\mathbf{P}_{k-1,i}(u_0) + u_0 \\mathbf{P}_{k-1,i+1}(u_0)
```
for \$ k = 1,...,n\$ and \$i = 0,...,n-k\$, where \$u_0\$ is any single value from 0 to 1.

(see NURBS, eqn 1.12 and A1.5)
"""
function decasteljau_bezier1D(P, u)
    n = length(P[:,1])
    Q = copy(P)
    for k=1:n
        for i=1:n-k
            Q[i,:] = (1-u)*Q[i,:] + u*Q[i+1,:]
        end
    end
    C = Q[1,:]

end