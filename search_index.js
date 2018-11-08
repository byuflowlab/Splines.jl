var documenterSearchIndex = {"docs": [

{
    "location": "Functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "Functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "The majority of the functions listed below come from The NURBS Book by Les Piegl and Wayne Tiller, hereafter referred to as \"NURBS\" when referencing equations and/or algorithms."
},

{
    "location": "Functions/#Splines.binomialcoeff",
    "page": "Functions",
    "title": "Splines.binomialcoeff",
    "category": "function",
    "text": "binomialcoeff(n, i)\n\nCalculate the Binomial Coefficient defined as:\n\nbinomni = fracni(n-1)\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.bernsteincoeff",
    "page": "Functions",
    "title": "Splines.bernsteincoeff",
    "category": "function",
    "text": "bernsteincoeff(u, n, i)\n\nCalculate Bernstein Coefficient (Bezier Basis Function) defined as:\n\nB_i n(u) = binomni u^i (1-u)^n-1\n\nat parametric point, u, where 0leq uleq1. u may either be a single value or an array.\n\n(see NURBS, eqn 1.8)\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.simple_bezier1D",
    "page": "Functions",
    "title": "Splines.simple_bezier1D",
    "category": "function",
    "text": "simple_bezier1D(P, u)\n\nCalculate a point along a Bezier curve at the parametric point, u, based on the control points mathbfP where the Bezier curve, mathbfC(u) is defined as:\n\nmathbfC(u) = sum_i=0^n B_i n(u) mathbfP_i  0 leq u leq 1\n\nwhere B is the basis (Bernstein Coefficient) at parametric point, u, as calculated frombernsteincoeff, and n is the number of control points in vector mathbfP. Again, u may either be a single value or an array.\n\n(see NURBS eqn 1.7)\n\n\n\n\n\n"
},

{
    "location": "Functions/#Bezier-Functions-1",
    "page": "Functions",
    "title": "Bezier Functions",
    "category": "section",
    "text": "Splines.binomialcoefflet n = 6 and i = 2. We calculate n choose i (which in this case is 15) by calling the binomialcoeff function:Exampleimport Splines # hide\nn = 6\ni = 2\nnchoosei = Splines.binomialcoeff(n, i)Splines.bernsteincoeffExamplesexample 1: single u valueimport Splines # hide\nu = 0.5\nn = 6\ni = 2\nb = Splines.bernsteincoeff(u, n, i)example 2: u as an arrayimport Splines #hide\nu = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]\nn = 6\ni = 2\nb = Splines.bernsteincoeff(u, n, i)Splines.simple_bezier1Dimport Splines\nusing Plots\n\n\nP = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition\nu = collect(0:0.05:1.0) #parametric points\n\nbezierCurve = Splines.simple_bezier1D(P, u)\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24, linewidths = 5)\nplot!(bezierCurve[:, 1], bezierCurve[:, 2], linewidths=10, aspectratio=:equal, grid=:off, label=\"Bezier\")\nplot!(P[:, 1], P[:, 2], markersizes=10, markershapes=:rect, linewidths=10, label=\"Control Points\")\nsavefig(\"simplebezier.svg\")ExampleP = [0.0 0.0; 0.0 0.1; 0.3 0.25; 1.0 0.0] #control point definition\nu = collect(0:0.05:1.0) #parametric points\n\nbezierCurve = Splines.simple_bezier1D(P, u)(Image: )"
},

{
    "location": "Functions/#Splines.getspanindex-NTuple{4,Any}",
    "page": "Functions",
    "title": "Splines.getspanindex",
    "category": "method",
    "text": "getspanindex(n, p, u, U)\n\nComplete binary search to find span index of vector, U, in which knot, u, lies. (NURBS A2.1)\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.basisfunctions-NTuple{4,Any}",
    "page": "Functions",
    "title": "Splines.basisfunctions",
    "category": "method",
    "text": "basisFunctions(i, u, p, U)\n\nCalculate the non-vanishing basis functions of the B-Spline of order p, defined by knots U at knot u.\n\nThe formula for the basis functions is:\n\nN_i0(u) =\nbegincases\n      1  textrmif  u_i leq u leq u_i+1 \n      0  textrmotherwise\nendcases\n\nN_ip(u) = fracu-u_iu_i+p - u_i N_ip-1(u) - fracu_i+p+1 - uu_i+p+1 - u_i+1 N_i+1p-1(u)\n\nNote that the algorithm used in basisFunctions removes redunant calculation and potential division by zero (see NURBS, eqn 2.5 and A2.2).\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.basisfunctionsderivatives-NTuple{5,Any}",
    "page": "Functions",
    "title": "Splines.basisfunctionsderivatives",
    "category": "method",
    "text": "basisfunctionsderivatives(i, u, p, n, U)\n\nCalculate the non-vanishing basis functions and derivatives of the B-Spline of order p, defined by knots U at parametric location u.\n\nThe basis function derivative is given by\n\nN_ip^ = fracpu_i+p - u_i N_ip-1(u) - fracpu_i+p+1 - u_i+1 N_i+1p-1(u)\n\n(see NURBS, eqn 2.7 and A2.3)\n\nInputs:\n\ni : knot span containing u\nu : parametric point of interest\np : the curve order\nn : the max derivative order (n ≦ p)\nU : the knot vector\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.curvederivatives1-NTuple{6,Any}",
    "page": "Functions",
    "title": "Splines.curvederivatives1",
    "category": "method",
    "text": "curvederivatives1(n, p, U, P, u, d)\n\nCompute a curve point and its derivatives up do the dth derivative at parametric point u. (NURBS, A3.2)\n\nInputs\n\nn : the number of control points is n+1\np : the degree of the curve\nU : the knot vector\nP : the control points\nu : the parametric point of interest\nd : derivative order (0 ≤ k ≦ d)\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.curvederivativecontrolpoints-NTuple{7,Any}",
    "page": "Functions",
    "title": "Splines.curvederivativecontrolpoints",
    "category": "method",
    "text": "curvederivativecontrolpoints(n, p, U, P, d, r1, r2)\n\nCompute control points of curve derivatives:\n\nmathbfC^(k)(u) = sum_i=0^n-kN_ip-k(u) mathbfP_i^(k)\n\nwith\n\nmathbfP_i^(k) =\nbegincases\n    mathbfP_i  k=0 \n    fracp-k+1u_i+p+1-u_i+kleft(mathbfP_i+1^(k) - mathbfP_i^(k) right)  k  0\nendcases\n\n(see NURBS, eqn 3.8 and A3.3)\n\nInputs\n\nn : the number of control points is n+1\np : the degree of the curve\nU : the knot vector\nP : the control points\nu : the parametric point of interest\nd : derivative order (0 ≤ k ≦ d)\nr1 : first control point index\nr2 : last control point index\n\n\n\n\n\n"
},

{
    "location": "Functions/#B-Spline-Functions-1",
    "page": "Functions",
    "title": "B-Spline Functions",
    "category": "section",
    "text": "Splines.getspanindex(n, p, u, U)import Splines\nU = [0,0,0,1,2,3,4,4,5,5,5]\np = 2\nn = length(U)-p-1ExampleU = [0,0,0,1,2,3,4,4,5,5,5]\np = 2\nn = length(U)-p-1\n\nu = 5/2\nSplines.getspanindex(n,p,u,U)u = 5\nSplines.getspanindex(n,p,u,U)u = 0\nSplines.getspanindex(n,p,u,U)Splines.basisfunctions(i, u, p, U)Exampleimport Splines # hide\nU = [0,0,0,1,2,3,4,4,5,5,5]\nu = 5/2\np = 2\ni = 5\nbases = Splines.basisfunctions(i,u,p,U)Splines.basisfunctionsderivatives(i, u, p, n, U)Exampleimport Splines # hide\nU = [0,0,0,1,2,3,4,4,5,5,5]\nu = 5/2\np = 2\ni = 5\nn = p\nderivatives = Splines.basisfunctionsderivatives(i,u,p,n,U)Splines.curvederivatives1(n, p, U, P, u, d)Exampleimport Splines # hide\nU = [0,0,0,1,2,3,4,4,5,5,5]\nu = 5/2\np = 2\nn = length(U)-p-1\nP = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]\nd = 1\ncurveDerivatives = Splines.curvederivatives1(n, p, U, P, u, d)Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)Exampleimport Splines # hide\nU = [0,0,0,0,2/5,3/5,3/5,1,1,1,1]\nu = 1/2\ni = 4\nd = 1\np = 3\nP = [0 0; 1/2 1/2; 1 0; 3/2 1/2; 2 0; 5/2 1/2; 3 0]\nn = length(P[:,1])-1\nr1 = 0\nr2 = n\ncprime = Splines.curvederivativecontrolpoints(n, p, U, P, d, r1, r2)"
},

{
    "location": "Functions/#Splines.nurbsbasis",
    "page": "Functions",
    "title": "Splines.nurbsbasis",
    "category": "function",
    "text": "nurbsbasis(i,p,u,U,w)\n\nGet rational basis functions and derivatives. see eqn 4.2\n\nR_ip(u) = fracN_ip(u)w_isum_j=0^n N_jp(u)w_j\n\nwhere N_ip(u ) are B-Spline Basis Functions and w_i are weights associated with the NURBS control points.\n\nInputs:\n\nu : parametric point of interest\np : the curve order\nn : the max derivative order (n ≦ p)\nU : the knot vector\nweights : control point weights\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.curvepoint",
    "page": "Functions",
    "title": "Splines.curvepoint",
    "category": "function",
    "text": "curvepoint(n, p, U, Pw, u)\n\nCompute point on rational B-Spline curve defined as:\n\nmathbfC^w(u) = sum_i=0^n N_i p(u) mathbfP_i^w\n\nwhere mathbfP_i^w are the set of weighted control points and weights such that mathbfP_i^w = (w_ix_i w_iy_i w_iz_i w_i).\n\n(see NURBS, eqn 4.5 and A4.1)\n\nInputs:\n\nn : the number of control points minus 1 (the index of the last control point)\np : the curve order\nU : the knot vector\nPw : the set of weighted control points and weights\nu : the parametric point of interest\n\nTODO: if u value outside of U vector range is given, function hangs, but doesn\'t throw error. Need to add a check/error.\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.rationalcurvederivatives",
    "page": "Functions",
    "title": "Splines.rationalcurvederivatives",
    "category": "function",
    "text": "rationalcurvederivatives(Aders, wders, d)\n\nCompute the point mathbfC(u) and the derivatives mathbfC^(k)(u) for 1 leq k leq d where:\n\nmathbfC^(k)(u) = frac mathbfA^(k)(u) - sum_i=1^k binomki w^(i)(u) mathbfC^(k-1)(u) w(u)\n\nwhere mathbfA^(k)(u) and w^(i)(u) are precomputed using preweighted control points for some parametric point, 0 leq u leq 1, from curvederivatives1 and are inputs Aders and wders, respectively.\n\n(see NURBS eqn 4.8 and A4.2)\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.curveknotinsertion",
    "page": "Functions",
    "title": "Splines.curveknotinsertion",
    "category": "function",
    "text": "curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\nCompute a new curve from knot insertion. Using the formula:\n\nmathbfQ_i r^w = alpha_i r mathbfQ_i r-1^w + (1-alpha_i r) mathbfQ_i-1 r-1^w\n\nwhere\n\nalpha_i r =\nbegincases\n     1  i leq k-p+r-1 \n     fracbaru - u_iu_i+p-r+1 - baru_i  k-p+r leq ileq k-s \n     0  i geq k-s+1\nendcases\n\n(see NURBS eqn 5.15 and A5.1)\n\nInputs:\n\nnp : the number of control points minus 1 (the index of the last control point) before insertion\np : the curve order\nUP : the knot vector before insertion\nPw : the set of weighted control points and weights before insertion\nu : the parametric point of interest\nk : the span index at which the knot is to be inserted.\ns : numer of instances of the new knot alrady present in the knot vector, UP\nr : number of times the new knot is inserted (it is assumed that r+s leq p )\n\nOutputs:\n\nnq : the number of control points minus 1 (the index of the last control point) after insertion\nUQ : the knot vector after insertion\nQw : the set of weighted control points and weights after insertion\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.refineknotvectorcurve",
    "page": "Functions",
    "title": "Splines.refineknotvectorcurve",
    "category": "function",
    "text": "refineknotvectorcurve(n, p, U, Pw, X, r)\n\nRefine curve knot vector using NURBS A5.4.\n\nThis algorithm is simply a knot insertion algorithm that allows for multiple knots to be added simulataneously, i.e., a knot refinement procedure.\n\nInputs:\n\nn : the number of control points minus 1 (the index of the last control point) before insertion\np : the curve order\nU : the knot vector before insertion\nPw : the set of weighted control points and weights before insertion\nX : elements, in ascending order, to be inserted into U (elements should be repeated according to their multiplicities, e.g., if x and y have multiplicites 2 and 3, X = [x,x,y,y,y])\nr : length of X vector\n\nOutputs:\n\nUbar : the knot vector after insertion\nQw : the set of weighted control points and weights after insertion\n\n\n\n\n\n"
},

{
    "location": "Functions/#Splines.degreeelevatecurve",
    "page": "Functions",
    "title": "Splines.degreeelevatecurve",
    "category": "function",
    "text": "degreeelevatecurve(n,p,U,Pw,t)\n\nRaise degree of spline from p to p +t, t geq 1 by computing the new control point vector and knot vector.\n\nKnots are inserted to divide the spline into equivalent Bezier Curves. These curves are then degree elevated using the following equation.\n\nmathbfP^t_i = sum^textrmmin(pi)_j=textrmmax(0i-t) fracbinompj binomti-j mathbfP_jbinomp+tii=0p+t\n\nwhere mathbfP^t_i are the degree elevated control points after t -degree elevations\n\nFinally, the excess knots are removed and the degree elevated spline is returned.\n\n(see NURBS eqn 5.36, A5.9)\n\nInputs:\n\nn : the number of control points minus 1 (the index of the last control point) before degree elevation\np : the curve order\nU : the knot vector before degree elevation\nPw : the set of weighted control points and weights before degree elevation\nt : the number of degrees to elevate, i.e. the new curve degree is p+t\n\nOutputs:\n\nnh : the number of control points minus 1 (the index of the last control point) after degree elevation\nUh : the knot vector after degree elevation\nQw : the set of weighted control points and weights after degree elevation\n\n\n\n\n\n"
},

{
    "location": "Functions/#NURBS-Functions-1",
    "page": "Functions",
    "title": "NURBS Functions",
    "category": "section",
    "text": "Splines.nurbsbasisExampleimport Splines #hide\nU = [0,0,0,1,2,3,4,4,5,5,5] #knot vector\nw = [1,1,1,1,1,1,1] #control point weights\nu = 5/2 #parametric point of interest\np = 2 #curve degree\nn = 1 #number of derivatives\nR, dR = Splines.nurbsbasis(u,p,n,U,w) #rational bases and first derivativesSplines.curvepointExamplesexample 1: single pointimport Splines # hide\nU = [0, 0, 0, 1, 1, 1] #knot vector\nu = 0 #parametric point of interest\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nCw = Splines.curvepoint(n, p, U, Pw, u)import Splines\nU = [0, 0, 0, 1, 1, 1] #knot vector\nu = collect(0:0.05:1.0) #parametric points\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nCw = zeros(length(u), length(Pw[1, :]))\nfor i = 1:length(u)\n  Cw[i, :] = Splines.curvepoint(n, p, U, Pw, u[i])\nend\n\n\nusing Plots\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(Cw[:, 1], Cw[:, 2], linewidths=10, aspectratio=:equal, grid=:off, label=\"NURBS\")\nplot!(P[:, 1], P[:, 2], markersizes=10, markershapes=:square, linewidths=10, label=\"Control Points\")\nsavefig(\"nurbscircle.svg\")example 2: array of pointsU = [0, 0, 0, 1, 1, 1] #knot vector\nu = collect(0:0.05:1.0) #parametric points\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nCw = zeros(length(u), length(Pw[1, :]))\nfor i = 1:length(u)\n  Cw[i, :] = Splines.curvepoint(n, p, U, Pw, u[i])\nend(Image: )Splines.rationalcurvederivativesExampleimport Splines # hide\nU = [0, 0, 0, 1, 1, 1] #knot vector\nu = 0 #parametric point of interest\np = 2 #curve order\nP = [1 0; 1 1; 0 1] #unweighted points\nw = [1 1 2] #weights\nPw = [1 0 1; 1 1 1; 0 2 2] #weighted points\nn = length(P[:, 1])-1\nd = 2 #max derivative level (2nd derivative)\n#Calculate Cw(u) derivatives\nders = Splines.curvederivatives1(n, p, U, Pw, u, d)\n#Separate derivatives\nAders = ders[:, 1:end-1]\nwders = ders[:, end]\n#Calculate NURBS derivatives\nCK = Splines.rationalcurvederivatives(Aders, wders, d)Splines.curveknotinsertionExamplesexample 1: Unique Knot Insertionimport Splines\nusing Plots\n\n\nUP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 5/2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 0\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\ncurvepoints = collect(0:0.1:5)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(np, p, UP, Pw, curvepoints[i])\nend\n\ncurvepoints = collect(0:0.1:5)\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nq, p, UQ, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP\'s\")\nplot!(Qw[4:6, 1], Qw[4:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label=\"New Control Points\")\nsavefig(\"uniqueknotinsert.svg\")UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 5/2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 0\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)(Image: )example 2: Repeated Knot Insertionimport Splines\nusing Plots\n\n\nUP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 1\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)\n\ncurvepoints = collect(0.0:0.1:5.0)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(np, p, UP, Pw, curvepoints[i])\nend\n\ncurvepoints = collect(0.0:0.1:5.0)\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nq, p, UQ, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP\'s\")\nplot!(Qw[4:6, 1], Qw[4:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label=\"New Control Points\")\nsavefig(\"repeatknotinsert.svg\")UP = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nu = 2\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nnp = length(P[:, 1])-1\nk = 5\ns = 1\nr = 1\n\nnq, UQ, Qw = Splines.curveknotinsertion(np, p, UP, Pw, u, k, s, r)(Image: )Splines.refineknotvectorcurveExampleimport Splines\nusing Plots\n\n\nU = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]\nX = [1.5, 2.5]\np = 3\nP = [0 0; 1 1; 2 0; 3 0; 4 1; 3 2; 2 2; 1.5 1]\nw = [1 1 1 1 1 1 1 1]\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1.5 1 1]\nn = length(P[:,1])-1\nr = length(X)-1\n\nUbar, Qwcalcd = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)\n\nnq, UQ, Qw = Splines.curveknotinsertion(n, p, U, Pw, X[1], 4, 0, 1)\nnq, UQ, Qw = Splines.curveknotinsertion(nq, p, UQ, Qw, X[2], 6, 0, 1)\n\ncurvepoints = collect(0:0.1:5)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n    Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])\nend\n\ncurvepoints = collect(0:0.1:5)\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n    Cw2[i, :] = Splines.curvepoint(nq, p, Ubar, Qwcalcd, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP\'s\")\nplot!(Qw[2:6, 1], Qw[2:6, 2], markersizes=10, markershapes=:circle, linestyle=:dash, linewidths=7, label=\"New Control Points\")\nsavefig(\"multiknotinsert.svg\")U = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5] #original knot vector\nX = [1.5, 2.5] #knots to be added (ascending order)\np = 3 #curve order\nPw = [0 0 1; 1 1 1; 2 0 1; 3 0 1; 4 1 1; 3 2 1; 2 2 1; 1 1 1] #weighted control points\nn = length(P[:,1])-1 #largest zero-based index in control point vector\nr = length(X)-1 #largest zero-based index in X (vector of knots to be added)\n\nUbar, Qw = Splines.refineknotvectorcurve(n, p, U, Pw, X, r)(Image: )Splines.degreeelevatecurveExamplesexample 1: Single Degree Elevationimport Splines\nusing Plots\n\n\nU = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 1\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)\n\ncurvepoints = collect(0:0.01:1)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])\nend\n\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP\'s\")\nplot!(Qw[2:end-1, 1], Qw[2:end-1, 2], markersizes=10, markershapes=:circle, linestyle=:dash,linewidths=7,  label=\"New Control Points\")\nsavefig(\"1degreeelevate.svg\")U = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 1\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)(Image: )example 2: 2 Degree Elevationimport Splines\nusing Plots\n\n\nU = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 2\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)\n\ncurvepoints = collect(0:0.01:1)\nCw1 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw1[i, :] = Splines.curvepoint(n, p, U, Pw, curvepoints[i])\nend\n\nCw2 = zeros(length(curvepoints), length(Pw[1, :]))\nfor i = 1:length(curvepoints)\n  Cw2[i, :] = Splines.curvepoint(nh, p+t, Uh, Qw, curvepoints[i])\nend\n\nplot(size=(2400,1600), titlefontsize=24, legendfontsize=24, tickfontsize=24, guidfontsize=24)\nplot!(linewidths=10, aspectratio=:equal, grid=:off)\nplot!(Cw1[:, 1], Cw1[:, 2], linewidths=10, label=\"Original Spline\")\nplot!(Cw2[:, 1], Cw2[:, 2], linestyle=:dot, linewidths=10, label=\"New Spline\")\nplot!(Pw[:, 1], Pw[:, 2], markersizes=10, markershapes=:square, linewidths=7, label=\"Original CP\'s\")\nplot!(Qw[2:end-1, 1], Qw[2:end-1, 2], markersizes=10, markershapes=:circle, linestyle=:dash,linewidths=7,  label=\"New Control Points\")\nsavefig(\"2degreeelevation.svg\")U = [0,0,0,0,3/10,7/10,1,1,1,1]\np = 3\nP = [-1 0; -1.5 1; -0.5 2; 0.5 2; 1.5 1; 1 0]\nw = [1 1 1 1 1 1]\nPw = [-1 0 1; -1.5 1 1; -0.5 2 1; 0.5 2 1; 1.5 1 1; 1 0 1]\nn = length(P[:,1])-1\nt = 2\n\nnh, Uh, Qw = Splines.degreeelevatecurve(n,p,U,Pw,t)(Image: )"
},

{
    "location": "#",
    "page": "Splines.jl Documentation",
    "title": "Splines.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#Splines.jl-Documentation-1",
    "page": "Splines.jl Documentation",
    "title": "Splines.jl Documentation",
    "category": "section",
    "text": "CurrentModule = SplinesPages = [\"Functions.md\"]\nDepth = 2"
},

]}
