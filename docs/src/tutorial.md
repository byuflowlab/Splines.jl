# Quick Start Guide

## B-Spline Objects

```@docs
Splines.BSpline
```

To begin, we need to choose a degree, knot vector, and control points to create a B-Spline.

!!! note "Matching Vector Lengths"
    For a given degree, p, you must repeat the first and last knots in the knot vector p+1 times.

    For a given number of control points, cp, there must be cp + 2p - 1 elements in the knot vector, including the repeated ones at the endpoints and any other repeated knots throughout.

    If you've matched the lengths of the vectors like this, then you should end up with a spline. If not, you'll end up with an error sooner or later.

```@example
import Splines # hide

degree = 2 # degree of B-Spline.

knots = [0.0,0.0,0.0,0.5,1.0,1.0,1.0] # knot vector, with endpoint knots repeated degree+1 times

ctrlpts = [[0.0, 0.0], [0.0, 0.125], [0.5, 0.125], [1.0, 0.0]] # control point vector, with knots - 2degree + 1 points

bspline = Splines.BSpline(degree,knots,ctrlpts) # BSpline object
```







## NURBS Objects
```@docs
Splines.BSpline
```

A NURBS Object is nearly identical to a BSpline object with one exception: NURBS objects include weights for the control points, one for each.

```@example
import Splines # hide

degree = 2 # degree of B-Spline.

knots = [0.0,0.0,0.0,0.5,1.0,1.0,1.0] # knot vector, with endpoint knots repeated degree+1 times

ctrlpts = [[0.0, 0.0], [0.0, 0.125], [0.5, 0.125], [1.0, 0.0]] # control point vector, with knots - 2degree + 1 points

weights = [1; 1; 1; 1] # control point weights (all 1's results in the equivalent of a BSPline)

nurbs = Splines.NURBS(degree,knots,weights,ctrlpts) # NURBS object
```