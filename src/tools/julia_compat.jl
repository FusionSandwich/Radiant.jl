"""
Julia-version compatibility helpers used by Radiant.

Radiant v1.1.67 declares Julia 1.6.7 compatibility, but the upstream source uses the
`range(first, last)` integer shorthand introduced after Julia 1.6 in more than one thousand loop
sites. Rewriting all transport and cross-section kernels in one unqualified mechanical change
would be high risk. On Julia versions older than 1.7, provide only the missing two-integer method
with the inclusive semantics used by those loops. Newer Julia versions retain Base's native
implementation.

This is deliberately restricted to `Integer` endpoints. It must not be generalized to floating
point endpoints because their spacing is ambiguous without an explicit `step` or `length`.
"""
@static if VERSION < v"1.7.0"
    function Base.range(first::Integer,last::Integer)
        return first:last
    end
end
