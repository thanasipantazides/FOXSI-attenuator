module Rotations
using LinearAlgebra

export so3
"""
    so3(u)

Get skew-symmetric matrix for crossing `u` with another 3-vector. I.e.: so3(u)*v == cross(u,v)
"""
function so3(u)
    return [0 -u[3] u[2]; u[3] 0 -u[1]; -u[2] u[1] 0]
end

export dcm1
"""
    dcm1(angle)

Returns the direction cosine matrix corresponding to a rotation about the 1-axis by the given `angle`.
"""
function dcm1(angle)
    return [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)]
end

export dcm2
"""
    dcm2(angle)

Returns the direction cosine matrix corresponding to a rotation about the 2-axis by the given `angle`.
"""
function dcm2(angle)
    return [cos(angle) 0 -sin(angle); 0 1 0; sin(angle) 0 cos(angle)]
end

export dcm3
"""
    dcm3(angle)

Returns the direction cosine matrix corresponding to a rotation about the 3-axis by the given `angle`.
"""
function dcm3(angle)
    return [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1]
end

export axisangle2dcm
"""
    axisangle2dcm(axis, angle)

Returns the direction cosine matrix corresponding to the attitude transformation described by the `axis` and `angle` supplied.
"""
function axisangle2dcm(axis, angle)
    return I*cos(angle) + (1 - cos(angle))*(axis*axis') + so3(axis)*sin(angle)
end

# export dcm2axisangle


end