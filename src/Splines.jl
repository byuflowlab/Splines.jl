module Splines

using LinearAlgebra: norm
using OffsetArrays

# convenience functions
onebased = OffsetArrays.no_offset_view
zerobased(arr) = OffsetArray(arr, 0:length(arr)-1)

# include("bezier.jl")

export BSpline
export curvepoint, curvederivatives, curvederivativecontrolpoints
export globalcurveinterpolation, leastsquarescurve
include("bspline.jl")

export NURBS  # also defines (overloaded) curvepoint, curvederivatives
include("nurbs.jl")

end # module
