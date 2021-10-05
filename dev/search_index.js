var documenterSearchIndex = {"docs":
[{"location":"howto/#Common-How-to's","page":"Guided Examples","title":"Common How-to's","text":"","category":"section"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"<!–","category":"page"},{"location":"howto/#Bezier-Functions","page":"Guided Examples","title":"Bezier Functions","text":"","category":"section"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.binomialcoeff","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nn = 6\ni = 2\nnchoosei = Splines.binomialcoeff(n, i)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.bernsteincoeff","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Examples","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 1: single u value","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nu = 0.5\nn = 6\ni = 2\nb = Splines.bernsteincoeff(u, n, i)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 2: u as an array","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines #hide\nu = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]\nn = 6\ni = 2\nb = Splines.bernsteincoeff(u, n, i)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.simple_bezier1D","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nusing Plots\n\n\nP = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition\nu = collect(0:0.05:1.0) #parametric points\n\nbezierCurve = Splines.simple_bezier1D(P, u)\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24, linewidths = 5)\nplot!(bezierCurve[:, 1], bezierCurve[:, 2], linewidths=10, aspectratio=:equal, grid=:off, label=\"Bezier\")\nplot!(P[:, 1], P[:, 2], markersizes=10, markershapes=:rect, linewidths=10, label=\"Control Points\")\nsavefig(\"simplebezier.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"P = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition\nu = collect(0:0.05:1.0) #parametric points\n\nbezierCurve = Splines.simple_bezier1D(P, u)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: )","category":"page"},{"location":"howto/#B-Spline-Functions","page":"Guided Examples","title":"B-Spline Functions","text":"","category":"section"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.getspanindex(n, p, u, U)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nU = [0,0,0,1,2,3,4,4,5,5,5]\np = 2\nn = length(U)-p-1","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"U = [0,0,0,1,2,3,4,4,5,5,5]\np = 2\nn = length(U)-p-1\n\nu = 5/2\nSplines.getspanindex(n,p,u,U)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"u = 5\nSplines.getspanindex(n,p,u,U)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"u = 0\nSplines.getspanindex(n,p,u,U)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.basisfunctions(i, u, p, U)","category":"page"},{"location":"howto/#Splines.basisfunctions-NTuple{4,Any}","page":"Guided Examples","title":"Splines.basisfunctions","text":"basisfunctions(span, deg, knots, u)\n\n(private function) Compute nonvanishing basis functions (NURBS A2.2)\n\nArguments\n\ndeg::Integer: degree\nknots::Vector{Float64}:: a knot vector (u0, ... un+1)\nu::Float64`: nondimensional location we are searching for\nspan::Integer: corresponding index i for u between knotsi and knotsi+1 (computed from getspanindex)\n\nReturns\n\nN::Vector{Float64}: vector of length N0 ... Ndeg\n\n\n\n\n\n","category":"method"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nU = [0,0,0,1,2,3,4,4,5,5,5]\nu = 5/2\np = 2\ni = 5\nbases = Splines.basisfunctions(i,u,p,U)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.basisfunctionsderivatives(i, u, p, n, U)","category":"page"},{"location":"howto/#Splines.basisfunctionsderivatives-NTuple{5,Any}","page":"Guided Examples","title":"Splines.basisfunctionsderivatives","text":"basisfunctionsderivatives(span, deg, knots, u, n)\n\n(private function) Calculate the non-vanishing basis functions and derivatives of the B-Spline of order p, defined by knots U at parametric point,u`. (NURBS A 2.3)\n\nArguments\n\nspan::Integer: knot span containing u\ndeg::Integer: the curve order\nknots::Vector{Float}: the knot vector\nu::Float: parametric point of interest\nn::Integer : the max derivative order (n ≦ p)\n\nReturns\n\nders::Matrix{Float}: [0..n, 0..p]  ders[0, :] function values, ders[1: :], first derivatives, etc.\n\n\n\n\n\n","category":"method"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nU = [0,0,0,1,2,3,4,4,5,5,5]\nu = 5/2\np = 2\ni = 5\nn = p\nderivatives = Splines.basisfunctionsderivatives(i,u,p,n,U)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.curvederivatives1(n, p, U, P, u, d)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nU = [0,0,0,1,2,3,4,4,5,5,5]\nu = 5/2\np = 2\nn = length(U)-p-1\nP = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]\nd = 1\ncurveDerivatives = Splines.curvederivatives1(n, p, U, P, u, d)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nU = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]\nu = 1/2\ni = 4\nd = 1\np = 3\nP = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]\nn = length(P[:,1])-1\nr1 = 0\nr2 = n\ncprime = Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.globalcurveinterpolation(n,Q,r,p; knotplacement)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nQ = [0 0; 3 4; -1 4; -4 0; -4 -3]\nr = 2\nn = 4\np = 3\nU = [0 0 0 0 28/51 1 1 1 1]\nm, U, P = Splines.globalcurveinterpolation(n,Q,r,p;knotplacement=\"chordlength\")","category":"page"},{"location":"howto/#NURBS-Functions","page":"Guided Examples","title":"NURBS Functions","text":"","category":"section"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.nurbsbasis","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines #hide\nU = [0,0,0,1,2,3,4,4,5,5,5] #knot vector\nw = [1,1,1,1,1,1,1] #control point weights\nu = 5/2 #parametric point of interest\np = 2 #curve degree\nn = 1 #number of derivatives\nR, dR = Splines.nurbsbasis(u,p,n,U,w) #rational bases and first derivatives","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.curvepoint","category":"page"},{"location":"howto/#Splines.curvepoint","page":"Guided Examples","title":"Splines.curvepoint","text":"Evaluate point on B-spline curve (NURBS, A3.1)\n\nArguments\n\nbspline::BSpline: bspline object\nu::Float: point on spline to evaluate at\n\nReturns\n\nC::Vector{Float}: point in ND space\n\n\n\n\n\nEvaluate point on rational b-spline curve (NURBS, A4.1)\n\nArguments\n\nnurbs::NURBS: NURBS object\nu::Float: point on spline to evaluate at\n\nReturns\n\nC::Vector{Float}: point in ND space\n\n\n\n\n\n","category":"function"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Examples","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 1: single point","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nU = [0, 0, 0, 1, 1, 1] #knot vector\nu = 0 #parametric point of interest\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nCw = Splines.curvepoint(n, p, U, Pw, u)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nU = [0, 0, 0, 1, 1, 1] #knot vector\nu = collect(0:0.05:1.0) #parametric points\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nCw = zeros(length(u), length(Pw[1, :]))\nfor i = 1:length(u)\n  Cw[i, :] = Splines.curvepoint(n, p, U, Pw, u[i])\nend\n\n\nusing Plots\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(Cw[:, 1], Cw[:, 2], linewidths=10, aspectratio=:equal, grid=:off, label=\"NURBS\")\nplot!(P[:, 1], P[:, 2], markersizes=10, markershapes=:square, linewidths=10, label=\"Control Points\")\nsavefig(\"nurbscircle.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 2: array of points","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"U = [0, 0, 0, 1, 1, 1] #knot vector\nu = collect(0:0.05:1.0) #parametric points\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nCw = zeros(length(u), length(Pw[1, :]))\nfor i = 1:length(u)\n  Cw[i, :] = Splines.curvepoint(n, p, U, Pw, u[i])\nend","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: )","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.rationalcurvederivatives","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines # hide\nU = [0, 0, 0, 1, 1, 1] #knot vector\nu = 0 #parametric point of interest\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nd = 2 #max derivative level (2nd derivative)\n#Calculate Cw(u) derivatives\nders = Splines.curvederivatives1(n, p, U, Pw, u, d)\n#Separate derivatives\nAders = ders[:, 1:end-1]\nwders = ders[:, end]\n#Calculate NURBS derivatives\nCK = Splines.rationalcurvederivatives(Aders, wders, d)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.curveknotinsertion","category":"page"},{"location":"howto/#Splines.curveknotinsertion","page":"Guided Examples","title":"Splines.curveknotinsertion","text":"curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\nCompute a new curve from knot insertion. Using the formula:\n\nmathbfQ_i r^w = alpha_i r mathbfQ_i r-1^w + (1-alpha_i r) mathbfQ_i-1 r-1^w\n\nwhere\n\nalpha_i r =\nbegincases\n     1  i leq k-p+r-1 \n     fracbaru - u_iu_i+p-r+1 - baru_i  k-p+r leq ileq k-s \n     0  i geq k-s+1\nendcases\n\n(see NURBS eqn 5.15 and A5.1)\n\nInputs:\n\nnp : the number of control points minus 1 (the index of the last control point) before insertion\np : the curve order\nUP : the knot vector before insertion\nPw : the set of weighted control points and weights before insertion\nu : the knot to be added\nk : the span index at which the knot is to be inserted.\ns : numer of instances of the new knot alrady present in the knot vector, UP\nr : number of times the new knot is inserted (it is assumed that r+s leq p )\n\nOutputs:\n\nnq : the number of control points minus 1 (the index of the last control point) after insertion\nUQ : the knot vector after insertion\nQw : the set of weighted control points and weights after insertion\n\n\n\n\n\n","category":"function"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Examples","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 1: Unique Knot Insertion","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nusing Plots\n\n\nUP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 5/2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 0\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\ncurvepoints = collect(0:0.1:5)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(np, p, UP, Pw, curvepoints[i])\nend\n\ncurvepoints = collect(0:0.1:5)\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nq, p, UQ, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP's\")\nplot!(Qw[4:6, 1], Qw[4:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label=\"New Control Points\")\nsavefig(\"uniqueknotinsert.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 5/2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 0\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: )","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 2: Repeated Knot Insertion","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nusing Plots\n\n\nUP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 1\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\ncurvepoints = collect(0.0:0.1:5.0)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(np, p, UP, Pw, curvepoints[i])\nend\n\ncurvepoints = collect(0.0:0.1:5.0)\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nq, p, UQ, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP's\")\nplot!(Qw[4:6, 1], Qw[4:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label=\"New Control Points\")\nsavefig(\"repeatknotinsert.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 1\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: )","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.refineknotvectorcurve","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Example","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nusing Plots\n\n\nU = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nX = [1.5, 2.5]\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nn = length(P[:,1])-1\nr = length(X)-1\n\nUbar, Qwcalcd = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)\n\nnq, UQ, Qw = Splines.curveknotinsertion(n, p, U, Pw, X[1], 4, 0, 1)\nnq, UQ, Qw = Splines.curveknotinsertion(nq, p, UQ, Qw, X[2], 6, 0, 1)\n\ncurvepoints = collect(0:0.1:5)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n    Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])\nend\n\ncurvepoints = collect(0:0.1:5)\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n    Cw2[i, :] = Splines.curvepoint(nq, p, Ubar, Qwcalcd, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP's\")\nplot!(Qw[2:6, 1], Qw[2:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label=\"New Control Points\")\nsavefig(\"multiknotinsert.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"U = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5] #original knot vector\nX = [1.5, 2.5] #knots to be added (ascending order)\np = 3 #curve order\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1] #weighted control points\nn = length(P[:,1])-1 #largest zero-based index in control point vector\nr = length(X)-1 #largest zero-based index in X (vector of knots to be added)\n\nUbar, Qw = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: )","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Splines.degreeelevatecurve","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"Examples","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 1: Single Degree Elevation","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nusing Plots\n\n\nU = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 1\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)\n\ncurvepoints = collect(0:0.01:1)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])\nend\n\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP's\")\nplot!(Qw[2:end-1, 1], Qw[2:end-1, 2], markersizes=10, markershapes=:circle, linestyle=:dash,linewidths=7,  label=\"New Control Points\")\nsavefig(\"1degreeelevate.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"U = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 1\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: )","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"example 2: 2 Degree Elevation","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"import Splines\nusing Plots\n\n\nU = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 2\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)\n\ncurvepoints = collect(0:0.01:1)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])\nend\n\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP's\")\nplot!(Qw[2:end-1, 1], Qw[2:end-1, 2], markersizes=10, markershapes=:circle, linestyle=:dash,linewidths=7,  label=\"New Control Points\")\nsavefig(\"2degreeelevation.svg\")","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"U = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 2\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)","category":"page"},{"location":"howto/","page":"Guided Examples","title":"Guided Examples","text":"(Image: ) –>","category":"page"},{"location":"license/#License","page":"About","title":"License","text":"","category":"section"},{"location":"license/","page":"About","title":"About","text":"The Splines.jl package is licensed under the MIT \"Expat\" License:","category":"page"},{"location":"license/","page":"About","title":"About","text":"Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:","category":"page"},{"location":"license/","page":"About","title":"About","text":"The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.","category":"page"},{"location":"license/","page":"About","title":"About","text":"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"reference/","page":"API Reference","title":"API Reference","text":"Modules = [Splines]\nOrder   = [:type, :function]","category":"page"},{"location":"reference/#Splines.BSpline","page":"API Reference","title":"Splines.BSpline","text":"BSpline(degree, knots, ctrlpts)\n\nConstruct a B-Spline object\n\nArguments\n\ndeg::Integer: degree\nknots::Vector{Float64}:: a knot vector (u0, ... un+1)\nctrlpts::Vector{Vector{Float64}}:: control points.  outer index is number of control points, inner index of dimensionality of point.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Splines.NURBS","page":"API Reference","title":"Splines.NURBS","text":"NURBS(degree, knots, weights, ctrlpts)\n\nConstruct a NURBS object\n\nArguments\n\ndeg::Integer: degree\nknots::Vector{Float64}:: a knot vector (u0, ... un+1)\nweights::Vector{Float64}:: a corresponding vector of weights\nctrlpts::Vector{Vector{Float64}}:: control points.  outer index is number of control points, inner index of dimensionality of point.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Splines.curvederivativecontrolpoints-NTuple{4,Any}","page":"API Reference","title":"Splines.curvederivativecontrolpoints","text":"curvederivativecontrolpoints(n, p, U, P, d, r1, r2)\n\nCompute control points of curve derivatives:\n\nmathbfC^(k)(u) = sum_i=0^n-kN_ip-k(u) mathbfP_i^(k)\n\nwith\n\nmathbfP_i^(k) =\nbegincases\n    mathbfP_i  k=0 \n    fracp-k+1u_i+p+1-u_i+kleft(mathbfP_i+1^(k) - mathbfP_i^(k) right)  k  0\nendcases\n\n(see NURBS, eqn 3.8 and A3.3)\n\nArguments\n\nbspline::BSpline: bspline object\n`r1::Integer : first control point index to take derivatives at\n`r2::Integer : last control point index to take derivatives at\nd::Float: derivative order (0 ≤ k ≦ d)\n\nReturns\n\nPK::Matrix{Vector{Float}}: PK[i, j] is the ith derivative of the jth control point\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.curvederivatives-Tuple{BSpline,Any,Any}","page":"API Reference","title":"Splines.curvederivatives","text":"curvederivatives(bspline, u, d)\n\nCompute a curve point and its derivatives up do the dth derivative at parametric point, u. (NURBS, A3.2)\n\nArguments\n\nbspline::BSpline: bspline object\nu::Float: point on spline to evaluate at\nd::Float: derivative order (0 ≤ k ≤ d)\n\nReturns\n\nCK::Vector{Vector{Float}}.  where CK[0] is the point, CK[1] the first derivative, and so on.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.curvederivatives-Tuple{NURBS,Any,Any}","page":"API Reference","title":"Splines.curvederivatives","text":"curvederivatives(nurbs, u, d)\n\nCompute a curve point and its derivatives up do the dth derivative at parametric point, u. (NURBS, A3.2)\n\nArguments\n\nnurbs::NURBS: bspline object\nu::Float: point on spline to evaluate at\nd::Float: derivative order (0 ≤ k ≤ d)\n\nReturns\n\nCK::Vector{Vector{Float}} where CK[1] is the point, CK[2] the first derivative, and so on.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.curveknotinsertion-Tuple{Any,Any,Any}","page":"API Reference","title":"Splines.curveknotinsertion","text":"curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\nCompute a new curve from knot insertion. Using the formula:\n\nmathbfQ_i r^w = alpha_i r mathbfQ_i r-1^w + (1-alpha_i r) mathbfQ_i-1 r-1^w\n\nwhere\n\nalpha_i r =\nbegincases\n     1  i leq k-p+r-1 \n     fracbaru - u_iu_i+p-r+1 - baru_i  k-p+r leq ileq k-s \n     0  i geq k-s+1\nendcases\n\n(see NURBS eqn 5.15 and A5.1)\n\nInputs:\n\nnp : the number of control points minus 1 (the index of the last control point) before insertion\np : the curve order\nUP : the knot vector before insertion\nPw : the set of weighted control points and weights before insertion\nu : the knot to be added\nk : the span index at which the knot is to be inserted.\ns : numer of instances of the new knot alrady present in the knot vector, UP\nr : number of times the new knot is inserted (it is assumed that r+s leq p )\n\nOutputs:\n\nnq : the number of control points minus 1 (the index of the last control point) after insertion\nUQ : the knot vector after insertion\nQw : the set of weighted control points and weights after insertion\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.curvepoint-Tuple{BSpline,Any}","page":"API Reference","title":"Splines.curvepoint","text":"Evaluate point on B-spline curve (NURBS, A3.1)\n\nArguments\n\nbspline::BSpline: bspline object\nu::Float: point on spline to evaluate at\n\nReturns\n\nC::Vector{Float}: point in ND space\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.curvepoint-Tuple{NURBS,Any}","page":"API Reference","title":"Splines.curvepoint","text":"Evaluate point on rational b-spline curve (NURBS, A4.1)\n\nArguments\n\nnurbs::NURBS: NURBS object\nu::Float: point on spline to evaluate at\n\nReturns\n\nC::Vector{Float}: point in ND space\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.get_ctrlpts-Tuple{Any}","page":"API Reference","title":"Splines.get_ctrlpts","text":"get_ctrlpts\n\nreturn control points of spline\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.get_degree-Tuple{Any}","page":"API Reference","title":"Splines.get_degree","text":"get_degree\n\nreturn degree of spline\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.get_knots-Tuple{Any}","page":"API Reference","title":"Splines.get_knots","text":"get_knots\n\nreturn knot vector of spline\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.get_weights-Tuple{Any}","page":"API Reference","title":"Splines.get_weights","text":"get_weights\n\nreturn weights of spline\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.getspanindex-Tuple{Any,Any,Any}","page":"API Reference","title":"Splines.getspanindex","text":"getspanindex(deg, knots, u)\n\n(private function) binary search to find span index of vector, knots, in which the parametric point, u, lies. (NURBS A2.1)\n\nArguments\n\ndeg::Integer: degree\nknots::Vector{Float64}:: a knot vector (u0, ... un+1)\nu::Float64`: nondimensional location we are searching for\n\nReturns\n\nspan::Integer: corresponding index i for u between knotsi and knotsi+1\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.globalcurveinterpolation-Tuple{Any,Any}","page":"API Reference","title":"Splines.globalcurveinterpolation","text":"globalcurveinterpolation(pts, deg)\n\nInterpolate ctrl points pts, with a B-Spline of degree deg. (NURBS A9.1)\n\nArguments\n\npts::Vector{Vector{Float}}: outer vector of length npts, inner vector of length dimension of space\ndeg::Integer: degree of B-spline\n\nOutputs:\n\nbspline::BSpline: bspline object\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.leastsquarescurve-Tuple{Any,Any,Any}","page":"API Reference","title":"Splines.leastsquarescurve","text":"leastsquarescurve(pts, nctrl, deg)\n\nLeast squares fit to provided points.  (NURBS section 9.4.1) TODO: Currently hard-coded to 2D data.\n\nArguments\n\npts::Vector{Vector{Float}}: data points\nnctrl::Integer: number of control points to use in fit\ndeg::Integer: degree of bspline used in fit\n\nReturns\n\nbspline::BSpline: a BSpline object\n\n\n\n\n\n","category":"method"},{"location":"reference/#Splines.singlebasisfunction-NTuple{4,Any}","page":"API Reference","title":"Splines.singlebasisfunction","text":"singlebasisfunction(i, deg, knots, u)\n\n(private function) Compute single basis function N_i^p. (NURBS A2.4)\n\nArguments\n\ni::Integer : the index of the basis (the i in N_i)\ndeg::Integer : the basis degree up to the the curve order\nknots::Vector{Float} : the knot vector\nu::Float : parametric point of interest\n\nReturns\n\nNip::Float: the N_i^p basis function value\n\n\n\n\n\n","category":"method"},{"location":"#Splines.jl","page":"Intro","title":"Splines.jl","text":"","category":"section"},{"location":"","page":"Intro","title":"Intro","text":"CurrentModule = Splines","category":"page"},{"location":"","page":"Intro","title":"Intro","text":"Splines.jl is a work-in-progress splines package focusing on Bézier curves, Basis Splines (B-Splines), and Non-uniform Rational B-Splines (NURBS). As of now, the available methods are transcribed from algorithms, or created from equations found in The NURBS Book by Les Piegl and Wayne Tiller, hereafter referred to as \"NURBS\" when referencing equations and/or algorithms.","category":"page"},{"location":"","page":"Intro","title":"Intro","text":"The algorithms outlined in the text are for C/C++ code. The C language, and the theory as presented is zero indexed. In an effort to preserve the algorithms as displayed in the text, and to keep the native Julia 1-indexing, indices are typically left as-is from the text with the addition of a '+1'.","category":"page"},{"location":"","page":"Intro","title":"Intro","text":"For FLOW Lab students, a quick intro to pertinent spline theory can be found in the FLOW Lab Notebook, but for in depth explanation The NURBS Book is the recommended resource.","category":"page"},{"location":"","page":"Intro","title":"Intro","text":"Note that this package has been created primarily for research purposes, so only basic tools required for that research have been implemented at this time. As the research progresses, more methods will be added as required.","category":"page"},{"location":"","page":"Intro","title":"Intro","text":"<!– ## Links to Function Descriptions Descriptions of available methods along with example implementations can be found on the functions page.","category":"page"},{"location":"","page":"Intro","title":"Intro","text":"Pages = [\"Functions.md\"]\nDepth = 2 -->","category":"page"},{"location":"#About","page":"Intro","title":"About","text":"","category":"section"},{"location":"","page":"Intro","title":"Intro","text":"License","category":"page"},{"location":"tutorial/#Quick-Start-Guide","page":"Quick Start","title":"Quick Start Guide","text":"","category":"section"},{"location":"tutorial/#B-Spline-Objects","page":"Quick Start","title":"B-Spline Objects","text":"","category":"section"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"Splines.BSpline","category":"page"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"To begin, we need to choose a degree, knot vector, and control points to create a B-Spline.","category":"page"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"note: Matching Vector Lengths\nFor a given degree, p, you must repeat the first and last knots in the knot vector p+1 times.For a given number of control points, cp, there must be cp + 2p - 1 elements in the knot vector, including the repeated ones at the endpoints and any other repeated knots throughout.If you've matched the lengths of the vectors like this, then you should end up with a spline. If not, you'll end up with an error sooner or later.","category":"page"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"import Splines # hide\n\ndegree = 2 # degree of B-Spline.\n\nknots = [0.0,0.0,0.0,0.5,1.0,1.0,1.0] # knot vector, with endpoint knots repeated degree+1 times\n\nctrlpts = [[0.0, 0.0], [0.0, 0.125], [0.5, 0.125], [1.0, 0.0]] # control point vector, with knots - 2degree + 1 points\n\nbspline = Splines.BSpline(degree,knots,ctrlpts) # BSpline object","category":"page"},{"location":"tutorial/#NURBS-Objects","page":"Quick Start","title":"NURBS Objects","text":"","category":"section"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"Splines.BSpline","category":"page"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"A NURBS Object is nearly identical to a BSpline object with one exception: NURBS objects include weights for the control points, one for each.","category":"page"},{"location":"tutorial/","page":"Quick Start","title":"Quick Start","text":"import Splines # hide\n\ndegree = 2 # degree of B-Spline.\n\nknots = [0.0,0.0,0.0,0.5,1.0,1.0,1.0] # knot vector, with endpoint knots repeated degree+1 times\n\nctrlpts = [[0.0, 0.0], [0.0, 0.125], [0.5, 0.125], [1.0, 0.0]] # control point vector, with knots - 2degree + 1 points\n\nweights = [1; 1; 1; 1] # control point weights (all 1's results in the equivalent of a BSPline)\n\nnurbs = Splines.NURBS(degree,knots,weights,ctrlpts) # NURBS object","category":"page"}]
}
