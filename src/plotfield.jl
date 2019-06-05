using PyPlot
using Gravityfield
"Rounding for legends"
round_sig(x; sig = 3) = round(x, digits = sig - Int(floor(log10(abs(x))))-1)


"Pyplot 3d surface and countour plot"
function wrap_surface_and_countour(np::Npos; title = "3d surface", levels = 100)
    wrap_surface_and_countour(slice_1(np)/m, slice_2(np)/m, slice_3(np)/m; title = title, levels = levels)
end
function wrap_surface_and_countour(mat_x, mat_y, mat_z; title = "3d surface", levels = 100)
    clf()
    using3D()
    fig = figure("pyplot_surfaceplot", figsize=(10, 10)) # inches
    ax = fig.add_subplot(2,1,1, projection="3d")
    nx, ny = size(mat_x)
    plot_surface(mat_x, mat_y, mat_z, rcount=nx, ccount = ny, edgecolors="k", 
                 cmap=ColorMap("Blues"), alpha=0.8, linewidth=0.25)
    xlabel("x")
    ylabel("y")
    PyPlot.title(title)
    #subplot(212)
    ax2 = fig.add_subplot(2,1,2, aspect = "equal")
    ax2.contourf(mat_x, mat_y, mat_z, levels = levels, cmap="Reds")
    cp = ax2.contour(mat_x, mat_y, mat_z, levels = levels, colors="k")
    ax2.clabel(cp, inline=1, fontsize=10, colors = "black")
    xlabel("x")
    ylabel("y")
    PyPlot.title(title * " contour")
    #tight_layout()
end
function wrap_quiver(x::Pos, y::Pos, np::Npos; title = "Quiver")
    if (size(x)[1], size(y)[1]) != size(np)
        error("incorrect dims")
    end 
    clf()
    fig = figure("pyplot_quiverplot", figsize = (10, 10))
    q = quiver(x'/m, y'/m, slice_1(np)'/m, slice_2(np)'/m)
    maxlen = maximum(len.(replace(np/m, [NaN, NaN]=>[0.0]))) |> round_sig
    ax = gca()
    ax.quiverkey(q, X=0.8, Y = 0.05, U = maxlen, coordinates="figure", label="Quiver key, length = $maxlen",labelpos = "E")
    PyPlot.title(title)
end
function wrap_stream_a(x::Matrix, y::Matrix, np::Npos; density = 1., title = "Stream a")
    clf()
    fig = figure("pyplot_streamplot", figsize = (10,10))
    streamplot(x = x/m, y = y/m, u = slice_1(np)/m, v = slice_2(np)/m,
            density = density)
    PyPlot.title(title)
end
function wrap_stream_b(x::Matrix, y::Matrix, np::Npos; density = 1., title = "Stream b")
    clf()
    fig = figure("pyplot_streamplot", figsize = (10,10))
    vl = len.(np)/m
    streamplot(x = x/m, y = y/m, u = slice_1(np)/m, v = slice_2(np)/m,
             color = vl,
             linewidth = 2,
             cmap = PyPlot.cm.autumn,
             density = density)
    PyPlot.title(title)
end
function wrap_stream_c(x::Matrix, y::Matrix, np::Npos; density = 1., title = "Stream c")
    clf()
    fig = figure("pyplot_streamplot", figsize = (10,10))
    vl = len.(np)
    maxlen = maximum(len.(replace(np, [NaN, NaN]m=>[0.0]m)))
    lw = 5 .* vl ./ maxlen
    streamplot(x = x/m, y = y/m, u = slice_1(np)/m, v = slice_2(np)/m,
        color = "k",
        density = density, 
        linewidth = lw)
    PyPlot.title(title)
end
function wrap_stream_d(x::Matrix, y::Matrix, np::Npos; density = 1., title = "Stream d")
    clf()
    fig = figure("pyplot_streamplot", figsize = (10,10))
    vl = len.(np)
    maxlen = maximum(len.(replace(np, [NaN, NaN]m=>[0.0]m)))
    lw = 5 .* vl ./ maxlen
    streamplot(x = x/m, y = y/m, u = slice_1(np)/m, v = slice_2(np)/m,
        color = vl/m,
        density = density,
        cmap = PyPlot.cm.YlGnBu,
        linewidth = lw)
    PyPlot.title(title)
end
nothing