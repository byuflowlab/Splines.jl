using Splines
using LinearAlgebra
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("beziertests.jl")
include("bsplinetests.jl")
include("nurbstests.jl")