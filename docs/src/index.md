# Splines.jl
```@meta
using Splines
```

Splines.jl is a work-in-progress splines package focusing on Bézier curves, Basis Splines (B-Splines), and Non-uniform Rational B-Splines (NURBS). As of now, the available methods are transcribed from algorithms, or created from equations found in *The NURBS Book* by Les Piegl and Wayne Tiller, hereafter referred to as "NURBS" when referencing equations and/or algorithms.

The algorithms outlined in the text are for C/C++ code. The C language, and the theory as presented is zero indexed. In an effort to preserve the algorithms as displayed in the text, and to keep the native Julia 1-indexing, indices are typically left as-is from the text with the addition of a '+1'.

For FLOW Lab students, a quick intro to pertinent spline theory can be found in the [FLOW Lab Notebook](https://github.com/byuflowlab/flowlab-notebook/blob/master/theory/splines/splines.pdf), but for in depth explanation *The NURBS Book* is the recommended resource.

Note that this package has been created primarily for research purposes, so only basic tools required for that research have been implemented at this time. As the research progresses, more methods will be added as required.



<!-- ## Links to Function Descriptions
Descriptions of available methods along with example implementations can be found on the functions page.

```@contents
Pages = ["Functions.md"]
Depth = 2 -->
```

## About
[License](license.md)
