# Quick Start Guide

## B-Spline Objects

```@docs
Splines.BSpline
```

To begin, we need to choose a degree, knot vector, and control points to create a B-Spline.

```@example
import Splines # hide

degree = 3 # degree of B-Spline.

knots = [0,0,0,1/3,2/3,1,1,1]

ctrlpts = [[0.0, 0.0], [0.0, 0.25], [0.75, 0.125], [1.0, 0.0]]

b = Splines.BSpline(degree,knots,ctrlpts)
```









## NURBS Objects