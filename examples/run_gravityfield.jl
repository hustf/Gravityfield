# ] activate .
# in vscode, display(gcf()) or just gcf() after each plot.
# VSCode use won't be interactive.
using Gravityfield
include(joinpath(@__DIR__, "src/plotfield.jl"))

"Plot a sphere for phun"
let
    n = 100
    u = range(0, 2π, length = n)
    v = range(0, π, length = n)
    x = cos.(u) * sin.(v)'m
    y = sin.(u) * sin.(v)'m
    z = ones(n) * cos.(v)'m
    np = nest(x, y, z)
    wrap_surface_and_countour(np, title = "Just a sphere")
end


# 3d space, plot gravity at a plane offset by a moon radius,
# crop to moon and earth
let 
    bodies = bodies_3d()
    x, y = xygrids( xmin = -d_moon_earth/5, 
                            xmax = d_moon_earth * 6/5, 
                            ymin = -d_moon_earth*(6/5 + 1/5)/2,
                            ymax = d_moon_earth*(6/5 + 1/5)/2,
                            n = 300)
    nxyz = nest(x, y, 1.45456e6m)
    np = nest(x, y, gfieldvals(nxyz, bodies)s²)
    wrap_surface_and_countour(np, title = "1 Gravity acceleration")
end

# 3d space, plot gravity in the planet and moon plane, 
# cap the plot roughly at moon's gravity.
let 
    bodies = bodies_3d()
    x,y = xygrids( xmin = -d_moon_earth/5, 
                            xmax = d_moon_earth * 6/5, 
                            ymin = -d_moon_earth*(6/5 + 1/5)/2,
                            ymax = d_moon_earth*(6/5 + 1/5)/2,
                            n = 300)
    nxyz = nest(x, y, 0m)
    np = nest(x, y, clamp.(gfieldvals(nxyz, bodies)s², 0m, 2m))
    wrap_surface_and_countour(np, title = "2 Clamped gravitational acceleration")
end

# 3d space, plot gravity in the planet and moon plane, 
# cap the plot at 0.5 m/s²
let 
    bodies = bodies_3d()
    x,y = xygrids( xmin = -d_moon_earth/5, 
                            xmax = d_moon_earth * 6/5, 
                            ymin = -d_moon_earth*(6/5 + 1/5)/2,
                            ymax = d_moon_earth*(6/5 + 1/5)/2,
                            n = 300)
    nxyz = nest(x, y, 0m)
    np = nest(x, y, clamp.(gfieldvals(nxyz, bodies)*s², 0m, 0.5m))
    wrap_surface_and_countour(np, title = "3 Clamped gravitational acceleration")
end

# 3d space, plot gravity in the planet and moon plane, 
# cap the plot at 0.5 m/s²
let 
    bodies = bodies_3d()
    x,y = xygrids( xmin = -d_moon_earth/5, 
                            xmax = d_moon_earth * 6/5, 
                            ymin = -d_moon_earth*(6/5 + 1/5)/2,
                            ymax = d_moon_earth*(6/5 + 1/5)/2,
                            n = 300)
    nxyz = nest(x, y, 0m)
    np = nest(x, y, clamp.(gfieldvals(nxyz, bodies)s², 0m, 0.5m))
    wrap_surface_and_countour(np, title = "4 Clamped gravitational acceleration")
end

# 3d space, plot gravity in the planet and moon plane, 
# cap the plot at 0.05 m/s², fewer levels.
let 
    bodies = bodies_3d()
    x,y = xygrids( xmin = -d_moon_earth/5, 
                            xmax = d_moon_earth * 6/5, 
                            ymin = -d_moon_earth*(6/5 + 1/5)/2,
                            ymax = d_moon_earth*(6/5 + 1/5)/2,
                            n = 300)
    nxyz = nest(x, y, 0m)
    np = nest(x, y, clamp.(gfieldvals(nxyz, bodies)s², 0m, 0.05m))
    wrap_surface_and_countour(np, 
                title = "5 Clamped gravitational acceleration",
                levels = 20)
end
# 3d space, plot gravity in the planet and moon plane, 
# cap the plot at 0.025 m/s², fewer levels.
# Zoom far out, both earth and moon. Zero gravity points visible.
let 
    bodies = bodies_3d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 300)
    nxyz = nest(x, y, 0m)
    np = nest(x, y, clamp.(gfieldvals(nxyz, bodies)s², 0m, 0.025m))
    wrap_surface_and_countour(np, 
                            title = "6 Clamped gravity",
                            levels = 20)
end

# 1d space, magnitude of g in one dimension, higher resolution.
# Zoom far out, both earth and moon. Zero gravity points visible,
# but hard to see.
let 
    bodies = bodies_1d()
    n = 10000
    nx = nest(range(-d_moon_earth * 1.5, d_moon_earth * 1.5, length = n))
    gmagn = gfieldvals(nx, bodies)
    clf()
    plot(slice_1(nx)/m, gmagn*s²/m)
    grid()
    title("7 1d gravity, moon, earth and zero gravity points")
end

# 1d space, magnitude of g in one dimension, higher resolution.
# Zoom far out, both earth and moon. Zero gravity points visible,
# capped gravity at 0.025 m/s²
let 
    bodies = bodies_1d()
    n = 1000
    nx = nest(range(-d_moon_earth * 1.5, d_moon_earth * 1.5, length = n))
    gmagn = clamp.(gfieldvals(nx, bodies)s², 0m, 0.025m)
    clf()
    plot(slice_1(nx)/m, gmagn/m)
    grid()
    title("8 1d gravity, moon, earth and zero gravity points")
    # Add a leader text
    dx = maximum(slice_1(nx)) - minimum(slice_1(nx))
    dy = maximum(gmagn) - minimum(gmagn)
    subsetx = slice_1(nx[1:floor(Int, n/2)])
    subsetg = gmagn[1:floor(Int, n/2)]
    y, xind = findmin(subsetg)
    x = subsetx[xind][1]
    annotate("0 m/s² at x = $(round(x/r_earth)) · r_earth",
        xy = [x; y]/m,
        xytext = [x + 0.1dx; y + 0.2dy]/m,
        xycoords = "data", arrowprops=Dict("facecolor"=>"blue"))
    # And one at the other minimum
    subsetx = slice_1(nx[floor(Int, n/2):end])
    subsetg = gmagn[floor(Int, n/2):end]
    y, xind = findmin(subsetg)
    x = subsetx[xind][1]
    annotate("0 m/s² at x = $(round(x/r_earth)) · r_earth",
        xy = [x; y]/m,
        xytext = [x - 0.3dx; y + 0.1dy]/m,
        xycoords = "data", arrowprops=Dict("facecolor"=>"blue"))
end
# 2d space, quiver plot gravity vectors in the planet and moon plane, 
# axes are vectors not grids.
let
    bodies = bodies_2d()
    n = 20
    x = collect(range(-r_earth*1.5, d_moon_earth+ 1.45456e6m * 1.5, length = n))
    dy = maximum(x) -minimum(x)
    y = collect(range(-dy/2, dy/2, length = n))
    mg= fill([0.,0.]m/s², n, n)
    for i = 1:n, j = 1:n
        mg[i, j] = gvec([x[i], y[j]], bodies)
    end
    wrap_quiver(x, y, mg*s², title = "9 Sun, moon, earth")
end

# 2d space, quiver plot gravity vectors in the planet and moon plane, 
# axes are vectors not grids. Cap vector length > 0.02.
let
    bodies = bodies_2d()
    n = 50
    x = collect(range(-d_moon_earth * 1.5, d_moon_earth * 1.5, length = n))
    dy = maximum(x) -minimum(x)
    y = collect(range(-dy/2, dy/2, length = n))
    mg = fill([0.,0.]m/s², n, n)
    for i = 1:n, j = 1:n
        g = gvec([x[i], y[j]], bodies)
        if len(g) < 0.02m/s²
            mg[i, j] = g
        else
            mg[i,j] = [NaN, NaN]m/s²
        end 
    end
    wrap_quiver(x, y, mg*s², title = "10 Sun, moon, earth gravity")
end

# 2d space, stream plot gravity vectors in the planet and moon plane.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 50)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    wrap_stream_a(x, y, np, title = "11 Sun, moon, earth gravity")
end

# 2d space, stream plot gravity vectors in the planet and moon plane.
# Colour lines based on field strength.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 50)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    wrap_stream_b(x, y, np, title = "12 Sun, moon, earth gravity")
end
# 2d space, stream plot gravity vectors in the planet and moon plane.
# Line widths based on field strength.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 50)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    wrap_stream_c(x, y, np, title = "13 Sun, moon, earth gravity")
end

# 2d space, stream plot gravity vectors in the planet and moon plane.
# Line widths based on field strength.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 50)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    wrap_stream_d(x, y, np, title = "14 Sun, moon, earth gravity")
end

# 2d space, stream plot gravity vectors in the planet and moon plane.
# Cap vector length > 0.02. Much denser grid.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 150)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    np = map(np) do g
        if len(g) < 0.02m
            g
        else
            [NaN, NaN]m
        end 
    end
    wrap_stream_a(x, y, np, title = "15 Sun, moon, earth gravity", density = 6.)
end

# 2d space, stream plot gravity vectors in the planet and moon plane.
# Cap vector length > 0.02. Much denser grid.
# Colour lines based on field strength.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 150)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    np = map(np) do g
        if len(g) < 0.02m
            g
        else
            [NaN, NaN]m
        end 
    end
    wrap_stream_b(x, y, np, title = "16 Sun, moon, earth gravity", density = 6.)
end
# 2d space, stream plot gravity vectors in the planet and moon plane.
# Cap vector length > 0.02. Much denser grid.
# Line widths based on field strength.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 150)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    np = map(np) do g
        if len(g) < 0.02m
            g
        else
            [NaN, NaN]m
        end 
    end
    wrap_stream_c(x, y, np, title = "17 Sun, moon, earth gravity", density = 3.)
end
# 2d space, stream plot gravity vectors in the planet and moon plane.
# Cap vector length > 0.02. Much denser grid.
# Line widths based on field strength.
# Colour lines based on field strength.
let
    bodies = bodies_2d()
    x, y = xygrids( xmin = -d_moon_earth * 1.5, 
                    xmax =  d_moon_earth * 1.5, 
                    ymin = -d_moon_earth * 1.5,
                    ymax = d_moon_earth * 1.5,
                    n = 150)
    nxy = nest(x, y)
    np = gvec(nxy, bodies)s²
    np = map(np) do g
        if len(g) < 0.02m
            g
        else
            [NaN, NaN]m
        end 
    end
    wrap_stream_d(x, y, np, title = "18 Sun, moon, earth gravity", density = 3.)
end

#using ForwardDif
