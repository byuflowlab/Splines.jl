module Splines

using LinearAlgebra: norm
using OffsetArrays
onebased = OffsetArrays.no_offset_view

include("bezier.jl")

export BSpline, curvepoint, curvederivatives, curvederivativecontrolpoints
export globalcurveinterpolation, leastsquarescurve
include("bspline.jl")

include("nurbs.jl")

end # module
