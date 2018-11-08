# Splines.jl

This package is an implementation of various algorithms outlined in "The NURBS Book" by Les Piegl and Wayne Tiller.

The algorithms outlined in the text are for C/C++ code. The C language, and the theory as presented is zero indexed. In an effort to preserve the algorithms as displayed in the text, and to keep the native Julia 1-indexing, indices are typically left as-is from the text with the addition of a +1.

Note that this package is still in early, active development for research purposes, so only tools required for that research have been implemented at this time.


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juddmehr.github.io/Splines.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juddmehr.github.io/Splines.jl/latest)