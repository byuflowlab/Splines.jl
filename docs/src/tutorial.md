# Quick Start Guide

## B-Spline Objects

```@docs
Splines.BSpline
```

To begin, we need to choose a degree, knot vector, and control points to create a B-Spline.

```@example
import Splines # hide

degree = 2 # degree of B-Spline.

knots = [0.0,0.0,0.0,0.5,1.0,1.0,1.0]

ctrlpts = [[0.0, 0.0], [0.5, 0.125], [1.0, 0.0]]

bspline = Splines.BSpline(degree,knots,ctrlpts)
```







## NURBS Objects