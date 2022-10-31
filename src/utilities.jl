
"""
    curvepointparameterized(spline, paramdim=1) -> Function

While [`curvepoint`](@Ref) returns the points corresponding to a given spline
parameter `u`, this function generates a function `f_curve_u(x)` that reverses
that relation to return the `u` that corresponds `x`, where x is the coordinate
in dimension `paramdim`.

# Arguments
- `spline: any spline object (e.g., Bspline, NURBS, etc)`
- `paramdim::Int`: dimension to parameterized

# Returns
- `f_curve_u::Function`: Inverted spline function returning `u` for a target
                            position `x`.
"""
function curvepointparameterized(spline, paramdim=1)

    """
        Return the spline parameter u for a given x (where x is the coordinate
        in dimension `paramdim`)
    """
    function f_curve_u(x)

        "Function to zero out"
        function f(u)
            # Simultaneous evaluation of x(u) and dxdu
            CK = Splines.curvederivatives(spline, u, 1)

            # Return f and f/f' (providing derivatives speeds up the Newton solve)
            return (x - CK[1][paramdim], -(x - CK[1][paramdim])/CK[2][paramdim])
        end

        # Netwon solve
        return Roots.newton(f, x)
    end

    return f_curve_u
end
