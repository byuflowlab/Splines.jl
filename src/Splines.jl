module Splines

using LinearAlgebra: norm
using OffsetArrays

import Roots

# convenience functions
onebased = OffsetArrays.no_offset_view
zerobased(arr) = OffsetArray(arr, 0:length(arr)-1)

include("Bezier.jl")

export BSpline
export curvepoint, curvederivatives, curvederivativecontrolpoints
export globalcurveinterpolation, leastsquarescurve
include("BSpline.jl")

export NURBS  # also defines (overloaded) curvepoint, curvederivatives
include("NURBS.jl")

end # module
