module Gravityfield_unitless
# ] activate .
"Newton's gravity constant"
const γ=6.67428E-11 #m^3*kg^-1*s^-2
const m_earth = 5.97E24 # kg
const r_earth = 6300000. # m
const r_sun = 695510000. #m
const d_sun_earth = 1.496e11  # m
const d_moon_earth = 384400000. # m
const m_sun = 1.989E30 # kg
const m_moon = 7.34767309e22
const yr = 31556926 # s
const ρ_sun = 1410. # kg/m^3

"Volume of a sphere"
V(r) = 4π / 3 * r^3
"Density of a sphere"
ρ(;m=1., r=1.) = m / V(r)
const ρ_earth = ρ(m = m_earth, r = r_earth)
"Position vector"
const Pos = Vector{T} where T<:Number
Pos(p) = Vector(p)
import Base.show
function show(io::IO, p::Pos) # short form
    typ = typeof(p)
    ioc = IOContext(io, :typeinfo => typ)
    show(ioc, ustrip(p))
    printstyled(ioc, unit(eltype(typ)); color=:cyan)
end
function show(io::IO, ::MIME"text/plain", p::Pos{T}) where T# long form
    typ = typeof(p)
    ioc = IOContext(io, :typeinfo => typ)
    print(ioc, "Pos{", T, "}(")
    show(ioc, ustrip(p))
    printstyled(ioc, unit(eltype(typ)); color=:cyan)
    print(ioc, ")")
end
"Nested positions. Each element is a position vector; the container can be a vector or matrix"
const Npos = Union{Array{S, 2}, Array{S,1}} where S<:Pos
Npos(p::Npos) = Array(p)
function show(io::IO, np::Npos) # short form
    np_unitless = map(p->ustrip(p), np)
    typ = typeof(np_unitless)
    ioc = IOContext(io, :typeinfo => typ)
    show(ioc, np_unitless)
    printstyled(ioc, unit(eltype(np[1])); color=:cyan)
    nothing
end
function show(io::IO, ::MIME"text/plain", np::Npos) # long form
    etyp = eltype(np[1])
    np_unitless = map(p->ustrip(p), np)
    typul = eltype(np_unitless[1])
    if ndims(np) == 1
        r= size(np)[1]
        println(io, r, "-element Npos{Pos{", typul, "}}:")
    else
        r, c = size(np)
        println(io, r, "x", c, " Npos{Pos{", typul, "}}:")
    end
    typ = typeof(np_unitless)
    ioc = IOContext(io, :typeinfo => typ)
    show(ioc, np_unitless)
    printstyled(io, unit(eltype(etyp)); color=:cyan)
end
const AA = AbstractArray
function nest(x::AA, y::AA, z::Number)
    Npos(map(x, y) do x, y
            [x, y, z]
        end)
end
function nest(x::AA, y::Number, z::Number)
    Npos(map(x) do x
            [x, y, z]
        end)
end
function nest(x::Number, y::AA, z::Number)
    Npos(map(y) do y
            [x, y, z]
        end)
end
function nest(x::Number, y::Number, z::AA)
    Npos(map(z) do z
            [x, y, z]
        end)
end

function nest(x::AA, y::Number)
    Npos(map(x) do x
            [x, y]
        end)
end
function nest(x::Number, y::AA)
    Npos(map(y) do y
            [x, y]
        end)
end
function nest(x::AA, y::AA, z::AA, α::AA, β::AA, γ::AA)
    Npos(map(x, y, z, α, β, γ) do x, y, z, α, β, γ
            [x, y, z, α, β, γ]
        end)
end
function nest(x::AA, y::AA, z::AA, α::AA, β::AA)
    Npos(map(x, y, z, α, β) do x, y, z, α, β
            [x, y, z, α, β]
        end)
end
function nest(x::AA, y::AA, z::AA, α::AA)
    Npos(map(x, y, z, α, β) do x, y, z, α
            [x, y, z, α]
        end)
end
function nest(x::AA, y::AA, z::AA)
    Npos(map(x, y, z) do x, y, z
            [x, y, z]
        end)
end
function nest(x::AA, y::AA)
    Npos(map(x, y) do x, y
            [x, y]
        end)
end
function nest(x::AA)
    Npos(map(x) do x
            [x]
        end)
end
"Slice out 1st component from a nested vector or matrix"
slice_1(np::Npos) = map(p-> p[1], np)
"Slice out 2nd component from a nested vector or matrix"
slice_2(np::Npos) = map(p-> p[2], np)
"Slice out 3rd component from a nested vector or matrix"
slice_3(np::Npos) = map(p-> p[3], np)
"Slice out 4th component from a nested vector or matrix"
slice_4(np::Npos) = map(p-> p[4], np)
"Slice out 5th component from a nested vector or matrix"
slice_5(np::Npos) = map(p-> p[5], np)
"Slice out sixth component from a nested vector or matrix"
slice_6(np::Npos) = map(p-> p[6], np)






abstract type Body end
"""
Constructor defaults to a spherical earth. Alternatively, use 
```
    earthlikeSphere(p = my_p, m = my_p)
    earthlikeSphere(p = my_p, r = my_r)
```
"""
Base.@kwdef struct Sphere <: Body
    p::Pos = [0., 0., 0.]
    m::Number = m_earth
    r::Number = r_earth
    Sphere(p, m, r) = m > zero(m) && r > zero(r) ? new(p, m, r) : error("no negative m or r")
end

show(io::IO, b::Sphere) = print(io, "Sphere(p = ", b.p, ", m = ", b.m, ", r = ", b.r, ")")

"Alternative constructors of Sphere for when mass or radius is unknown"
earthlikeSphere_m_known(;p = [0, 0, 0], m = m_earth) = Sphere(p, m , (3m / (4π * ρ_earth))^(1/3))
earthlikeSphere_r_known(;p = [0, 0, 0], r = r_earth) = Sphere(p, ρ_earth * V(r), r)

import Base.show
function show(io::IO, b::Sphere)
    ioc = IOContext(io, :typeinfo => typeof(b.p))
    print(ioc, "Sphere(p = ", b.p, ", m = ", b.m, ", r = ", b.r, ")")
end


const Bodies = Vector{T} where T <:Body
Bodies(bds::Bodies) = Vector(bds)
function bodies_3d()
    earth = Sphere()
    moon = earthlikeSphere_m_known(p = [d_moon_earth, 0m, 0m], m = m_moon)
    sun = Sphere([-d_sun_earth, 0, 0], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    Bodies([sun, earth, moon])
end
function bodies_2d()
    earth = Sphere(p = [0,0]m)
    moon = earthlikeSphere_m_known(p = [d_moon_earth, 0m], m = m_moon)
    sun = Sphere([-d_sun_earth, 0], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    Bodies([sun, earth, moon])
end
function bodies_1d()
    earth = Sphere(p = [0]m)
    moon = earthlikeSphere_m_known(p = [d_moon_earth], m = m_moon)
    sun = Sphere([-d_sun_earth], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    Bodies([sun, earth, moon])
end

import Base.show_circular
import Base._truncate_at_width_or_chars
function show(io::IO, ::MIME"text/plain", bds::Bodies) # long form.
    # show more descriptively, with one line per body. Based on Dict.
    recur_io = IOContext(io, :SHOWN_SET => bds)
    limit::Bool = get(io, :limit, false)
    if !haskey(io, :compact)
        recur_io = IOContext(recur_io, :compact => true)
    end
    summary(io, bds)
    isempty(bds) && return
    print(io, ":")
    show_circular(io, bds) && return
    if limit
        sz = displaysize(io)
        rows, cols = sz[1] - 3, sz[2]
        rows < 2   && (print(io, " …"); return)
        cols < 12  && (cols = 12) # Minimum widths of 2 for key, 4 for value
        cols -= 6 # Subtract the widths of prefix "  " separator " => "
        rows -= 1 # Subtract the summary

        # determine max key width to align the output, caching the strings
        ps = Vector{AbstractString}(undef, min(rows, length(bds)))
        ms = Vector{AbstractString}(undef, min(rows, length(bds)))
        rs = Vector{AbstractString}(undef, min(rows, length(bds)))
        plen = 0
        mlen = 0
        rlen = 0
        for (i, b) in enumerate(bds)
            i > rows && break
            st = sprint(show, b.p, context=recur_io, sizehint=0)
            ps[i] = "p = " * st 
            st = sprint(show,  b.m, context=recur_io, sizehint=0)
            ms[i] = " m = " * st
            st = sprint(show, b.r, context=recur_io, sizehint=0)
            rs[i] = " r = " * st
            plen = clamp(length(ps[i]), plen, cols)
            mlen = clamp(length(ms[i]), mlen, cols)
            rlen = clamp(length(rs[i]), rlen, cols)
        end
        if plen > max(div(cols, 2), cols - mlen - rlen)
            plen = max(cld(cols, 3), cols - mlen - rlen)
        end
    else
        rows = cols = typemax(Int)
    end

    for (i, bd) in enumerate(bds)
        print(io, "\n  ")
        i == rows < length(bds) && (print(io, rpad("⋮", plen), ", ⋮"); break)
        if limit
            pos = rpad(_truncate_at_width_or_chars(ps[i], plen, "\r\n"), plen)
        else
            st = sprint(show, bd.p, context=recur_io, sizehint=0)
            pos = "p = " * st
        end
        print(recur_io, pos)
        if limit
            val = rpad(_truncate_at_width_or_chars(ms[i], cols - plen, "\r\n"), mlen)
            print(io, val)
            val = _truncate_at_width_or_chars(rs[i], cols - plen - mlen, "\r\n")
            print(io, val)
        else
            st = " m = " * sprint(show, bd.m, context=recur_io, sizehint=0) *
                " r = " * sprint(show, bd.r, context=recur_io, sizehint=0)
            print(recur_io, st)
        end
    end
end



"Vector from position (of) a to position (of) b."
rvec(a::T, b::T) where T<: Body = b.p - a.p
rvec(a::Pos, b::Body) = b.p - a
rvec(a::Body, b::Pos) = b - a.p
rvec(a::T, b::T) where T<: Pos = b - a
"Unity vector of vector"
uvec(p::Pos) = len(p) > zero(eltype(p)) ? p / len(p) : zero(p)
"Unity vector from position (of) a to position (of) b."
uvec(a::T, b::T) where T<: Body = uvec(b.p - a.p)
uvec(a::Pos, b::Body) = uvec(b.p - a)
uvec(a::Body, b::Pos) = uvec(b - a.p)
uvec(a::T, b::T) where T<: Pos = uvec(b - a)


"Euclidean vector length"
len(p) =  sqrt(sum(p.*p))
"Euclidean length between positions (of) a and b"
len(a, b) = len(rvec(a, b))

"Gravitational acceleration magnitude"
g(m::Number, r::Number) = r > zero(r) ? γ * m / r^2 : zero(γ) / one(r)^2  # Test thoroughly for overflow

"Gravity vector on thisbody or position due to otherbody or bodies"
gvec(p::Pos, body::Body) = uvec(p, body) * g(body.m, len(p, body))
gvec(thisbody::Body, otherbody::Body) = gvec(thisbody.p, otherbody)
function gvec(p::Pos, sphere::Sphere)
    l = len(p, sphere)
    if l > sphere.r || l == 0.
        return uvec(p, sphere) * g(sphere.m, l)
    else
        # Inside sphere
        return uvec(p, sphere) * g(sphere.m, sphere.r) * l / sphere.r
    end
end
function gvec(p::Pos, bodies::Bodies) 
    gvec_one = b -> gvec(p, b)  
    sum(gvec_one.(bodies))
end
gvec(body::Body, bodies::Bodies) = gvec(body.p, bodies)
gvec(vp::Npos, bodies::Bodies) = map(p -> gvec(p, bodies), vp)


"Gravity field strength [length / time^2], nested position vectors"
gfieldvals(vp::Npos, bodies::Bodies) = map(p -> len(gvec(p, bodies)), vp)

"Return matrices mat_x, mat_y, size n * n"
function xygrids(;xmin = -3., xmax = 3., ymin, ymax, n = 100)
    mat_x = repeat(range(xmin, xmax, length = n)', n, 1)
    mat_y = repeat(range(ymin, ymax, length = n), 1, n)
    mat_x, mat_y
end
"100x100 grids around origo"
sqgrids(mi, ma) = xygrids( xmin = mi,xmax = ma, ymin = -(ma-mi)/2, ymax = (ma-mi)/2, n = 100)

nothing
end