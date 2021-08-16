# Common How-to's

<!--
## Bezier Functions
```@docs
Splines.binomialcoeff
```

**Example**
```@example
import Splines # hide
n = 6
i = 2
nchoosei = Splines.binomialcoeff(n, i)
```



```@docs
Splines.bernsteincoeff
```

**Examples**

*example 1: single u value*
```@example
import Splines # hide
u = 0.5
n = 6
i = 2
b = Splines.bernsteincoeff(u, n, i)
```
*example 2: u as an array*
```@example
import Splines #hide
u = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
n = 6
i = 2
b = Splines.bernsteincoeff(u, n, i)
```



```@docs
Splines.simple_bezier1D
```

```@setup simplebez
import Splines
using Plots


P = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition
u = collect(0:0.05:1.0) #parametric points

bezierCurve = Splines.simple_bezier1D(P, u)

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24, linewidths = 5)
plot!(bezierCurve[:, 1], bezierCurve[:, 2], linewidths=10, aspectratio=:equal, grid=:off, label="Bezier")
plot!(P[:, 1], P[:, 2], markersizes=10, markershapes=:rect, linewidths=10, label="Control Points")
savefig("simplebezier.svg")
```
**Example**
```@example simplebez
P = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition
u = collect(0:0.05:1.0) #parametric points

bezierCurve = Splines.simple_bezier1D(P, u)
```

![](simplebezier.svg)






## B-Spline Functions
```@docs
Splines.getspanindex(n, p, u, U)
```

```@setup findspan
import Splines
U = [0,0,0,1,2,3,4,4,5,5,5]
p = 2
n = length(U)-p-1
```

**Example**
```@example findspan
U = [0,0,0,1,2,3,4,4,5,5,5]
p = 2
n = length(U)-p-1

u = 5/2
Splines.getspanindex(n,p,u,U)
```

```@example findspan
u = 5
Splines.getspanindex(n,p,u,U)
```

```@example findspan
u = 0
Splines.getspanindex(n,p,u,U)
```


```@docs
Splines.basisfunctions(i, u, p, U)
```

**Example**
```@example
import Splines # hide
U = [0,0,0,1,2,3,4,4,5,5,5]
u = 5/2
p = 2
i = 5
bases = Splines.basisfunctions(i,u,p,U)
```


```@docs
Splines.basisfunctionsderivatives(i, u, p, n, U)
```

**Example**
```@example
import Splines # hide
U = [0,0,0,1,2,3,4,4,5,5,5]
u = 5/2
p = 2
i = 5
n = p
derivatives = Splines.basisfunctionsderivatives(i,u,p,n,U)
```

```@docs
Splines.curvederivatives1(n, p, U, P, u, d)
```

**Example**
```@example
import Splines # hide
U = [0,0,0,1,2,3,4,4,5,5,5]
u = 5/2
p = 2
n = length(U)-p-1
P = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]
d = 1
curveDerivatives = Splines.curvederivatives1(n, p, U, P, u, d)
```

```@docs
Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)
```

**Example**
```@example
import Splines # hide
U = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]
u = 1/2
i = 4
d = 1
p = 3
P = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]
n = length(P[:,1])-1
r1 = 0
r2 = n
cprime = Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)
```

```@docs
Splines.globalcurveinterpolation(n,Q,r,p; knotplacement)
```
**Example**
```@example
import Splines # hide
Q = [0 0; 3 4; -1 4; -4 0; -4 -3]
r = 2
n = 4
p = 3
U = [0 0 0 0 28/51 1 1 1 1]
m, U, P = Splines.globalcurveinterpolation(n,Q,r,p;knotplacement="chordlength")
```





## NURBS Functions

```@docs
Splines.nurbsbasis
```
*Example*
```@example
import Splines #hide
U = [0,0,0,1,2,3,4,4,5,5,5] #knot vector
w = [1,1,1,1,1,1,1] #control point weights
u = 5/2 #parametric point of interest
p = 2 #curve degree
n = 1 #number of derivatives
R, dR = Splines.nurbsbasis(u,p,n,U,w) #rational bases and first derivatives
```

```@docs
Splines.curvepoint
```
**Examples**

*example 1: single point*
```@example
import Splines # hide
U = [0, 0, 0, 1, 1, 1] #knot vector
u = 0 #parametric point of interest
p = 2 #curve order
P = [1 0; 1 1; 0 1] #unweighted points
w = [1 1 2] #weights
Pw = [1 0 1; 1 1 1; 0 2 2] #weighted points
n = length(P[:, 1])-1
Cw = Splines.curvepoint(n, p, U, Pw, u)
```

```@setup nurbs_circle
import Splines
U = [0, 0, 0, 1, 1, 1] #knot vector
u = collect(0:0.05:1.0) #parametric points
p = 2 #curve order
P = [1 0; 1 1; 0 1] #unweighted points
w = [1 1 2] #weights
Pw = [1 0 1; 1 1 1; 0 2 2] #weighted points
n = length(P[:, 1])-1
Cw = zeros(length(u), length(Pw[1, :]))
for i = 1:length(u)
  Cw[i, :] = Splines.curvepoint(n, p, U, Pw, u[i])
end


using Plots

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)
plot!(Cw[:, 1], Cw[:, 2], linewidths=10, aspectratio=:equal, grid=:off, label="NURBS")
plot!(P[:, 1], P[:, 2], markersizes=10, markershapes=:square, linewidths=10, label="Control Points")
savefig("nurbscircle.svg")
```

*example 2: array of points*
```@example nurbs_circle
U = [0, 0, 0, 1, 1, 1] #knot vector
u = collect(0:0.05:1.0) #parametric points
p = 2 #curve order
P = [1 0; 1 1; 0 1] #unweighted points
w = [1 1 2] #weights
Pw = [1 0 1; 1 1 1; 0 2 2] #weighted points
n = length(P[:, 1])-1
Cw = zeros(length(u), length(Pw[1, :]))
for i = 1:length(u)
  Cw[i, :] = Splines.curvepoint(n, p, U, Pw, u[i])
end
```

![](nurbscircle.svg)



```@docs
Splines.rationalcurvederivatives
```
**Example**

```@example
import Splines # hide
U = [0, 0, 0, 1, 1, 1] #knot vector
u = 0 #parametric point of interest
p = 2 #curve order
P = [1 0; 1 1; 0 1] #unweighted points
w = [1 1 2] #weights
Pw = [1 0 1; 1 1 1; 0 2 2] #weighted points
n = length(P[:, 1])-1
d = 2 #max derivative level (2nd derivative)
#Calculate Cw(u) derivatives
ders = Splines.curvederivatives1(n, p, U, Pw, u, d)
#Separate derivatives
Aders = ders[:, 1:end-1]
wders = ders[:, end]
#Calculate NURBS derivatives
CK = Splines.rationalcurvederivatives(Aders, wders, d)
```



```@docs
Splines.curveknotinsertion
```

**Examples**

*example 1: Unique Knot Insertion*
```@setup unique_knot_insertion
import Splines
using Plots


UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
u = 5/2
p = 3
P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]
w = [1 1 1 1 1 1 1 1]
Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]
np = length(P[:, 1])-1
k = 5
s = 0
r = 1

nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)

curvepoints = collect(0:0.1:5)
Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw1[i, :] = Splines.curvepoint(np, p, UP, Pw, curvepoints[i])
end

curvepoints = collect(0:0.1:5)
Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw2[i, :] = Splines.curvepoint(nq, p, UQ, Qw, curvepoints[i])
end

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)
plot!(linewidths=10, aspectratio=:equal, grid=:off)
plot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label="Original Spline")
plot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label="New Spline")
plot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label="Original CP's")
plot!(Qw[4:6, 1], Qw[4:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label="New Control Points")
savefig("uniqueknotinsert.svg")
```

```@example unique_knot_insertion
UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
u = 5/2
p = 3
P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]
w = [1 1 1 1 1 1 1 1]
Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]
np = length(P[:, 1])-1
k = 5
s = 0
r = 1

nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)
```

![](uniqueknotinsert.svg)

*example 2: Repeated Knot Insertion*
```@setup unique_knot_insertion
import Splines
using Plots


UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
u = 2
p = 3
P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]
w = [1 1 1 1 1 1 1 1]
Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]
np = length(P[:, 1])-1
k = 5
s = 1
r = 1

nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)

curvepoints = collect(0.0:0.1:5.0)
Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw1[i, :] = Splines.curvepoint(np, p, UP, Pw, curvepoints[i])
end

curvepoints = collect(0.0:0.1:5.0)
Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw2[i, :] = Splines.curvepoint(nq, p, UQ, Qw, curvepoints[i])
end

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)
plot!(linewidths=10, aspectratio=:equal, grid=:off)
plot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label="Original Spline")
plot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label="New Spline")
plot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label="Original CP's")
plot!(Qw[4:6, 1], Qw[4:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label="New Control Points")
savefig("repeatknotinsert.svg")
```

```@example unique_knot_insertion
UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
u = 2
p = 3
P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]
w = [1 1 1 1 1 1 1 1]
Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]
np = length(P[:, 1])-1
k = 5
s = 1
r = 1

nq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)
```

![](repeatknotinsert.svg)


```@docs
Splines.refineknotvectorcurve
```

**Example**

```@setup multi_knot_insertion
import Splines
using Plots


U = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
X = [1.5, 2.5]
p = 3
P = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]
w = [1 1 1 1 1 1 1 1]
Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]
n = length(P[:,1])-1
r = length(X)-1

Ubar, Qwcalcd = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)

nq, UQ, Qw = Splines.curveknotinsertion(n, p, U, Pw, X[1], 4, 0, 1)
nq, UQ, Qw = Splines.curveknotinsertion(nq, p, UQ, Qw, X[2], 6, 0, 1)

curvepoints = collect(0:0.1:5)
Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
    Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])
end

curvepoints = collect(0:0.1:5)
Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
    Cw2[i, :] = Splines.curvepoint(nq, p, Ubar, Qwcalcd, curvepoints[i])
end

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)
plot!(linewidths=10, aspectratio=:equal, grid=:off)
plot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label="Original Spline")
plot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label="New Spline")
plot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label="Original CP's")
plot!(Qw[2:6, 1], Qw[2:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label="New Control Points")
savefig("multiknotinsert.svg")
```

```@example multi_knot_insertion
U = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5] #original knot vector
X = [1.5, 2.5] #knots to be added (ascending order)
p = 3 #curve order
Pw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1] #weighted control points
n = length(P[:,1])-1 #largest zero-based index in control point vector
r = length(X)-1 #largest zero-based index in X (vector of knots to be added)

Ubar, Qw = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)
```

![](multiknotinsert.svg)


```@docs
Splines.degreeelevatecurve
```

**Examples**

*example 1: Single Degree Elevation*
```@setup 1-degree-elevation
import Splines
using Plots


U = [0,0,0,0,3/10,7/10,1,1,1,1]
p = 3
P = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]
w = [1 1 1 1 1 1]
Pw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]
n = length(P[:,1])-1
t = 1

nh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)

curvepoints = collect(0:0.01:1)
Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])
end

Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])
end

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)
plot!(linewidths=10, aspectratio=:equal, grid=:off)
plot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label="Original Spline")
plot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label="New Spline")
plot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label="Original CP's")
plot!(Qw[2:end-1, 1], Qw[2:end-1, 2], markersizes=10, markershapes=:circle, linestyle=:dash,linewidths=7,  label="New Control Points")
savefig("1degreeelevate.svg")
```

```@example 1-degree-elevation
U = [0,0,0,0,3/10,7/10,1,1,1,1]
p = 3
P = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]
w = [1 1 1 1 1 1]
Pw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]
n = length(P[:,1])-1
t = 1

nh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)
```

![](1degreeelevate.svg)

*example 2: 2 Degree Elevation*
```@setup 2-degree-elevation
import Splines
using Plots


U = [0,0,0,0,3/10,7/10,1,1,1,1]
p = 3
P = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]
w = [1 1 1 1 1 1]
Pw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]
n = length(P[:,1])-1
t = 2

nh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)

curvepoints = collect(0:0.01:1)
Cw1 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])
end

Cw2 = zeros(length(curvepoints), length(Pw[1, :]))
for i = 1:length(curvepoints)
  Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])
end

plot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)
plot!(linewidths=10, aspectratio=:equal, grid=:off)
plot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label="Original Spline")
plot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label="New Spline")
plot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label="Original CP's")
plot!(Qw[2:end-1, 1], Qw[2:end-1, 2], markersizes=10, markershapes=:circle, linestyle=:dash,linewidths=7,  label="New Control Points")
savefig("2degreeelevation.svg")
```

```@example 2-degree-elevation
U = [0,0,0,0,3/10,7/10,1,1,1,1]
p = 3
P = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]
w = [1 1 1 1 1 1]
Pw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]
n = length(P[:,1])-1
t = 2

nh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)
```

![](2degreeelevation.svg) -->
