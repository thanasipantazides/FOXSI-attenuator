module Attenuator3D
using Distributed
using SharedArrays
using LinearAlgebra
using Plots

import SpecialFunctions

export Particle
struct Particle
    r0  # initial position
    v   # propagation direction (unit)
    E   # energy
end

export Cylinder
struct Cylinder
    R   # radius
    h   # height
    c   # center point
    a   # axis (unit) (NOTE: c + h*a yields a point on the opposite cap => this is an inward normal axis)
end

export PixelatedAttenuator
struct PixelatedAttenuator
    holes
    density     # in kg/m^3
    toppoint    # a point (x,y,z) on the one side of the attenuator (in m)
    bottompoint # a point on the other side of the attenuator
    normal      # normal vector to attenuator surface
    largeradius # largest hole radius in attenuator 
    massattenuation # table of mass attenuation coefficients per energy
end

export so3
"""
    so3(u)

Get skew-symmetric matrix for crossing `u` with another 3-vector. I.e.: so3(u)*v == cross(u,v)
"""
function so3(u)
    return [0 -u[3] u[2]; u[3] 0 -u[1]; -u[2] u[1] 0]
end

export interpolateattenuation
"""
    interpolateattenuation(E, massattenuation)

Interpolate in array `massattenuation` by energy E (in eV). First column of `massattenuation` is energy in MeV, second column is mass attenuation coefficient μ/ρ in m^2/kg. Returns mass attenuation coefficient corresponding to input energy E.
"""
function interpolateattenuation(E, massattenuation)
    energies = massattenuation[:,1]
    coefficients = massattenuation[:,2]

    diff = E .- energies

    if all(diff .> 0) || all(diff .< 0)
        error("interpolant outside range")
    end

    upper = findfirst(diff .≤ 0)
    lower = findlast(diff .≥ 0)

    if lower == upper
        return coefficients[lower]
    else
        return coefficients[lower] + (E - energies[lower])/(energies[upper] - energies[lower])*(coefficients[upper] - coefficients[lower])
    end
end

export cylinderentryexit
"""
    cylinderentryexit(particle, cylinder)

Compute times at which a ray passes through a cylinder. The ray has initial position `particle.r0` and velocity `particle.v`. The cylinder has radius `cylinder.R`, height `cylinder.h`, a point one one of its caps `cylinder.c`, and a unit-length axis pointing from `cylinder.c` towards its other cap `cylinder.a`. 

Two scalars are returned: the time at which the ray enters the cylinder and the time at which it leaves the cylinder. If the ray misses the cylinder, zeros are returned for both entry and exit time.
"""
function cylinderentryexit(particle, cylinder)
    r0 = particle.r0
    v = particle.v
    c = cylinder.c
    a = cylinder.a
    R = cylinder.R
    h = cylinder.h

    vdir = v/norm(v)

    ltop = norm(v)*(c'*a - r0'*a)/(v'*a)
    lbottom = norm(v)*(c'*a + h - r0'*a)/(v'*a)

    topcapradius = norm(r0 + vdir*ltop - c)
    bottomcapradius = norm(r0 + vdir*lbottom - c - a*h)

    qa = 1 - (v'*a/norm(v))^2
    qb = 2*(r0' - c' + (a'*c - a'*r0)*a')*vdir
    qc = norm(r0 - c)^2 - (r0'*a - c'*a)^2 - R^2

    # particle coaxial with cylinder:
    if qa == 0
        # check if it hits caps, if so return caps times
        if topcapradius ≤ R
            if bottomcapradius ≤ R
                # rectify to zero if particle starts within cylinder
                ltop = max(ltop, 0)
                lbottom = max(lbottom, 0)
                return (min(ltop, lbottom), max(ltop, lbottom))
            else
                error("particle supposedly parallel to cylinder axis")
            end
        else
            # outside cap radius, miss
            if bottomcapradius ≤ R
                error("particle supposedly parallel to cylinder axis")
            end

            return (0, 0)
        end
        
    # photon oblique to cylinder axis:
    else
        root1 = (-qb + sqrt(Complex(qb^2 - 4*qa*qc)))/(2*qa)
        root2 = (-qb - sqrt(Complex(qb^2 - 4*qa*qc)))/(2*qa)

        if imag(root1) == 0
            root1 = real(root1)
            root2 = real(root2)
            # hits infinite cylinder somewhere
            l1 = min(root1, root2)
            l2 = max(root1, root2)

            # list of all intersection lengths with cap planes and cylinder
            dl = [lbottom, ltop, l1, l2]
            # order intersection times
            I = sortperm(dl)

            # check if first two intersection lengths are with cylinder (t1, t2)
            if length(intersect(I[1:2], [3, 4])) == 2
                # miss
                return (0, 0)
            end
            # check if last two intersection lengths are with cylinder (t1, t2)
            if length(intersect(I[3:4], [3, 4])) == 2
                # miss
                return (0, 0)
            end

            # if you've made it this far, cylinder is hit somewhere
            lin = dl[I[2]]
            lout = dl[I[3]]

            lin = max(lin, 0)
            lout = max(lout, 0)
            return (lin, lout)

        # miss altogether (complex roots for intersection with cylinder) 
        else
            return (0, 0)
        end
    end
end

export lengthsinattenuator
"""
    lengthsinattenuator(photon, attenuator)

Get distances a `photon::Particle` travels within an `attenuator::PixelatedAttenuator` which has cylinders cut out of it. 
"""
function lengthsinattenuator(photon, attenuator)
    topcylinderslengths = Array{Float64}(undef, 0, 2)
    # filter cylinders along path
    hits = falses(length(attenuator.holes[:]))
    for c = 1:length(attenuator.holes[:])
        if norm(cross(photon.v,attenuator.holes[c].a)) == 0
            # if cylinder parallel to ray
            hits[c] = norm(attenuator.holes[c].c - photon.r0 - photon.v*((photon.v'*(attenuator.holes[c].c - photon.r0))./(photon.v'*photon.v))) ≤ attenuator.holes[c].R
        else
            # if cylinder oblique to ray
            hits[c] = cross(photon.v, attenuator.holes[c].a)'*(photon.r0 - attenuator.holes[c].c)/norm(cross(photon.v, attenuator.holes[c].a)) ≤ attenuator.holes[c].R
        end
        
        # if this photon flies within this cylinder's radius, check for collision
        if hits[c]
            (lin, lout) = cylinderentryexit(photon, attenuator.holes[c])

            # make sure it actually hits cylinder
            if lout - lin > 0
                topcylinderslengths = cat(topcylinderslengths, [lin lout], dims=1)
            end
        end
    end
    
    # sort cylinders visited in order of visit time (increasing distance from r0 along v)
    sortedcylinderslengths = zeros(size(topcylinderslengths))
    if length(topcylinderslengths[:,1]) > 1
        topcylinderorder = sortperm(topcylinderslengths[:,1])
        sortedcylinderslengths = topcylinderslengths[topcylinderorder,:]
    else
        sortedcylinderslengths = topcylinderslengths
    end
    
    # lengths at which photon hits attenuator top surface and bottom surface
    lslabbottom = (attenuator.bottompoint'*attenuator.normal - photon.r0'*attenuator.normal)/(photon.v'*attenuator.normal)*norm(photon.v)
    lslabtop = (attenuator.toppoint'*attenuator.normal - photon.r0'*attenuator.normal)/(photon.v'*attenuator.normal)*norm(photon.v)

    lslabin = min(lslabbottom, lslabtop)
    lslabout = max(lslabbottom, lslabtop)

    if length(sortedcylinderslengths) > 0
        # times spent within attenuator material (complement of time spent in cylinders)
        lengthdiffs = [sortedcylinderslengths[:,1]; lslabout] - [lslabin; sortedcylinderslengths[:,2]]
    else
        lengthdiffs = lslabout - lslabin
    end

    # path lengths in material
    lengths = lengthdiffs

    return lengths
end

export absorptionprobability
"""
    absorptionprobability(photon, attenuator)

Return the probability of a photon::Particle being absorbed in attenuator::PixelatedAttenuator.
"""
function absorptionprobability(photon, attenuator)
    lengths = lengthsinattenuator(photon, attenuator)
    if length(lengths) != 0
        transmissionlikelihoods = zeros(size(lengths))
        for i = 1:length(lengths)
            # interpolate attenuation coefficients from table data:
            attenuationcoeff = interpolateattenuation(photon.E, attenuator.massattenuation)
            
            # probability of passing through material:
            transmissionlikelihoods[i] = exp(-lengths[i]/attenuationcoeff)
        end

        # complement of total transmission likelihood is absorption likelihood
        absorptionlikelihood = 1 - prod(transmissionlikelihoods)
        
        return absorptionlikelihood
    else
        return 0.0
    end
end

export batchphotons
"""
    batchphotons(photons, attenuator)

Compute absorption likelihood for a set of `photons` passing through `attenuator`. If multiple processes are available, computation will be parallelized. To execute in parallel, use the `Distributed` module in the calling context, add processes using `addprocs(n)`, and include this source with the `@everywhere macro` before `using` the module:
    
    @everywhere include("path_to_this_file.jl")
    using .Attenuator3D
"""
function batchphotons(photons, attenuator)
    # check if multiple processes available:
    if nprocs() > 1
        absorblikelihood = SharedArray{Float64}(size(photons))
    
        # parallel loop over photons
        @inbounds @sync @distributed for i = 1:length(photons)
            absorblikelihood[i] = absorptionprobability(photons[i], attenuator)
        end
        return Array(absorblikelihood)
    else
        absorblikelihood = zeros(size(photons))

        @inbounds for i = 1:length(photons)
            absorblikelihood[i] = absorptionprobability(photons[i], attenuator)
        end
        return absorblikelihood
    end
end

export airyintensity
"""
    airyintensity(radius, angle, wavenumber)

Get intensity fraction for diffracted light of given `wavenumber` projected at `angle` through hole of `radius` (following Airy's law for large distances).
"""
function airyintensity(radius, angle, wavelength)
    x = 2*π*radius/wavelength*sin(angle)
    I = (2*SpecialFunctions.besselj1(x)/x)^2
    return I
end

export plotcyl
"""
    plotcyl(cylinder::Cylinder)

Plot a transparent cylinder described by in Cylinder struct
"""
function plotcyl(cylinder)
    R = cylinder.R
    h = cylinder.h
    c = cylinder.c
    a = cylinder.a
    a = a./norm(a)

    # angular resolution:
    n = 100

    t = range(0, stop=2π, length=n)
    
    X = R.*[cos.(t) cos.(t)]
    Y = R.*[sin.(t) sin.(t)]
    Z = [zeros(n,1) h*ones(n,1)]

    transform = zeros(3,3)
    if a'*[0;0;1] == 1
        # if cylinder axis parallel to z-axis:
        transform = I
    elseif a'*[0;0;1] == -1
        # if cylinder axis antiparallel to z-axis:
        transform = [1 0 0; 0 -1 0; 0 0 -1]
    else
        # cylinder axis oblique to z-axis:
        ax = cross([0;0;1], a)
        ax = ax./norm(ax)
        cosang = a'*[0;0;1]
        # Rodrigues's formula/axis-angle representation:
        transform = I*cosang + (1 - cosang)*(ax*ax') + [0 -ax[3] ax[2]; ax[3] 0 -ax[1]; -ax[2] ax[1] 0]*sqrt(1 - cosang^2)
    end

    # transform all cylinder points:
    for i = 1:n
        posl = [X[i,1]; Y[i,1]; Z[i,1]]
        posr = [X[i,2]; Y[i,2]; Z[i,2]]

        newposl = transform*posl
        newposr = transform*posr

        X[i,:] = [newposl[1] newposr[1]] + [c[1] c[1]]
        Y[i,:] = [newposl[2] newposr[2]] + [c[2] c[2]]
        Z[i,:] = [newposl[3] newposr[3]] + [c[3] c[3]]
    end

    # draw surface:
    s = surface(X,Y,Z,
            color=:green,
            alpha=0.75,
            legend=false)
    return s
end

export plotcyl!
"""
    plotcyl!(cylinder::Cylinder)

Plot a transparent cylinder described by in Cylinder struct
"""
function plotcyl!(cylinder)
    R = cylinder.R
    h = cylinder.h
    c = cylinder.c
    a = cylinder.a
    a = a./norm(a)

    # angular resolution:
    n = 100

    t = range(0, stop=2π, length=n)

    X = R.*[cos.(t) cos.(t)]
    Y = R.*[sin.(t) sin.(t)]
    Z = [zeros(n,1) h*ones(n,1)]

    transform = zeros(3,3)
    if a'*[0;0;1] == 1
        # if cylinder axis parallel to z-axis:
        transform = I
    elseif a'*[0;0;1] == -1
        # if cylinder axis antiparallel to z-axis:
        transform = [1 0 0; 0 -1 0; 0 0 -1]
    else
        # cylinder axis oblique to z-axis:
        ax = cross([0;0;1], a)
        ax = ax./norm(ax)
        cosang = a'*[0;0;1]
        # Rodrigues's formula/axis-angle representation:
        transform = I*cosang + (1 - cosang)*(ax*ax') + so3(ax)*sqrt(1 - cosang^2)
    end

    # transform all cylinder points:
    for i = 1:n
        posl = [X[i,1]; Y[i,1]; Z[i,1]]
        posr = [X[i,2]; Y[i,2]; Z[i,2]]

        newposl = transform*posl
        newposr = transform*posr

        X[i,:] = [newposl[1] newposr[1]] + [c[1] c[1]]
        Y[i,:] = [newposl[2] newposr[2]] + [c[2] c[2]]
        Z[i,:] = [newposl[3] newposr[3]] + [c[3] c[3]]
    end

    # draw surface:
    s = surface!(X,Y,Z,
            color=:green,
            alpha=0.75,
            legend=false)
    return s
end

export plotparticles
"""
    plotparticles(particles::Array{Particle}, vscale, energycode)

Plot particles as rays with option (set by `energycode` bit) to colorcode. Option `vscale` scales length of velocity.
"""
function plotparticles(particles, vscale, energycode)
    
    n = length(particles)

    r0 = zeros(3, n)
    v = zeros(3, n)

    for i = 1:n
        r0[:, i] = particles[i].r0
        v[:, i] = particles[i].r0 .+ vscale*particles[i].v
    end

    data = cat(r0, v, dims=3)
    data = cat(data, fill(NaN,3,n), dims=3)
    x = vec(data[1,:,:]')
    y = vec(data[2,:,:]')
    z = vec(data[3,:,:]')

    p = plot(x, y, z, color=:yellow, legend=false)
    return p
end

export plotparticles!
"""
    plotparticles!(particles::Array{Particle}, energycode)

Plot particles as rays with option (set by `energycode` bit) to colorcode. Option `vscale` scales length of velocity.
"""
function plotparticles!(particles, vscale, energycode)
    
    n = length(particles)

    r0 = zeros(3, n)
    v = zeros(3, n)

    for i = 1:n
        r0[:, i] = particles[i].r0
        v[:, i] = particles[i].r0 .+ vscale*particles[i].v
    end

    data = cat(r0, v, dims=3)
    data = cat(data, fill(NaN,3,n), dims=3)
    x = vec(data[1,:,:]')
    y = vec(data[2,:,:]')
    z = vec(data[3,:,:]')

    p = plot!(x, y, z, color=:yellow, legend=false)
    return p
end

export plotattenuator
"""
    plotattenuator(attenuator::PixelatedAttenuator)

Plot bounding box of attenuator.
"""
function plotattenuator(attenuator)
    cornerstop = [attenuator.holes[1,1,1] attenuator.holes[1,end,1]; 
                  attenuator.holes[end,1,1] attenuator.holes[end,end,1]]
    cornersbottom = [attenuator.holes[1,1,end] attenuator.holes[1,end,end]; 
                    attenuator.holes[end,1,end] attenuator.holes[end,end,end]]
    corners = cat(cornerstop, cornersbottom, dims=3)

    a = corners[1,1,1].c'*attenuator.normal
    b = corners[1,1,2].c'*attenuator.normal

    extremes = zeros(2,2,2,3)

    if a > b
        # front = a
        if corners[1,1,1].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                    extremes[i,j,2,:] = corners[i,j,2].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c
                    extremes[i,j,2,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                end
            end
        end
    else
        # front = b
        if corners[1,1,2].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                    extremes[i,j,2,:] = corners[i,j,1].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c
                    extremes[i,j,2,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                end
            end
        end
    end

    s = surface(0,0,0,legend=false)
    for i = 1:3
        for j = 1:2
            if i == 1
                surf = extremes[j,:,:,:]
            elseif i == 2
                surf = extremes[:,j,:,:]
            elseif i == 3
                surf = extremes[:,:,j,:]
            end 
            surface!(surf[:,:,1], surf[:,:,2], surf[:,:,3], color=:purple, alpha=0.5, legend=false)
        end
    end
    
    current()
    return s
end

export plotattenuator
"""
    plotattenuator!(attenuator::PixelatedAttenuator)

Plot bounding box of attenuator.
"""
function plotattenuator!(attenuator)
    cornerstop = [attenuator.holes[1,1,1] attenuator.holes[1,end,1]; 
                  attenuator.holes[end,1,1] attenuator.holes[end,end,1]]
    cornersbottom = [attenuator.holes[1,1,end] attenuator.holes[1,end,end]; 
                    attenuator.holes[end,1,end] attenuator.holes[end,end,end]]
    corners = cat(cornerstop, cornersbottom, dims=3)

    a = corners[1,1,1].c'*attenuator.normal
    b = corners[1,1,2].c'*attenuator.normal

    extremes = zeros(2,2,2,3)

    if a > b
        # front = a
        if corners[1,1,1].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                    extremes[i,j,2,:] = corners[i,j,2].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c
                    extremes[i,j,2,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                end
            end
        end
    else
        # front = b
        if corners[1,1,2].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                    extremes[i,j,2,:] = corners[i,j,1].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c
                    extremes[i,j,2,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                end
            end
        end
    end

    for i = 1:3
        for j = 1:2
            if i == 1
                surf = extremes[j,:,:,:]
            elseif i == 2
                surf = extremes[:,j,:,:]
            elseif i == 3
                surf = extremes[:,:,j,:]
            end 
            s = surface!(surf[:,:,1], surf[:,:,2], surf[:,:,3], color=:purple, alpha=0.5, legend=false)
            current()
            return s
        end
    end
end

end # /module