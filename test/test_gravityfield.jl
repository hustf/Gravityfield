using Test
using Gravityfield

let
    @test V(1.0m) == (4π / 3)m³
    @test ρ(m = (4π / 3)kg, r = 1m) == 1kg/m³
    @test ρ(r = 1m, m = (4π / 3)kg) == 1kg/m³
    @test ρ_earth |> ustrip |> typeof == Float64
end

@testset "Gravity acceleration" begin
    @test g(1, 0) == 0.0m^3*kg^-1*s^-2
    @test g(1., 1.) == γ
    @test g(m_earth, r_earth) ≈ 10.0391664m/s^2
    ω_earth_sun = 2π / yr
    a_earth_sun = ω_earth_sun^2 * d_sun_earth
    @test g(m_sun, d_sun_earth) / a_earth_sun ≈ 1.000 atol = 0.001
end

@testset "Force vector on 1 due to 2" begin
    p1 = [0.,0.,0.]m
    p2 = [0.,1.,0.]m
    m1 = 316227kg
    m2 = 129099kg
    f1_2 = m1* uvec(p1, p2) * g(m2, len(p2-p1))
    @test f1_2 |> length == 3
    @test ustrip.(f1_2) |> typeof == Vector{Float64}
    @test f1_2 ≈ [0.0, 2.7247474, 0.]N
end
@testset "Euclidean vector length" begin
    @test len(-1m) == 1m
    @test len([1]m) == 1m
    @test len([2,0]m) == 2m
    @test len([1m,1m]) == sqrt(2)m
end
@testset "Vector between bodies or points" begin
    p1 = [0.,0.,0.]m
    p2 = [0.,1.,0.]m
    s = Sphere(p1, m_earth, r_earth)
    @test rvec(s,s) == [0., 0., 0.]m
    @test len(s,s) == 0m
    @test uvec(s,s) == [0., 0., 0.]
    @test rvec(p2,s) == [0., -1.0, 0.]m
    @test len(p2,s) == 1m
    @test uvec(p2,s) == [0., -1.0, 0.]

    s2 = Sphere(p2, m_earth, r_earth)
    @test rvec(s,s2) == [0., 1.0, 0.]m
    @test len(s,s2) == 1m
    @test uvec(s,s2) == [0., 1.0, 0.]

    s3 = Sphere([-d_sun_earth], 1, 1)
    s4 = Sphere([0.0m], 1, 1)
    @test rvec(s3,s4) == [d_sun_earth]
    @test len(s4,s3) == d_sun_earth
    @test uvec(s3,s4) == [1]
    @test uvec(s4,s3) == [-1]
end


@testset "Gravity acceleration vector btw. bodies or point and body" begin
    p1 = [0.,0.,0.]m
    s1 = Sphere(p = p1, m = 1kg, r = 1m)
    @test gvec(s1,s1) == [0., 0., 0.]m/s²
    @test gvec([-1, 0., 0.]m, s1) == [γ*kg/m², 0m/s², 0m/s²]
    s2 = Sphere([0.0]m, m_earth, r_earth)
    @test gvec([r_earth], s2) == [-10.039166439909296]m/s²
    @test gvec([r_earth/2], s2) == [-10.039166439909296]m/s²/2
end

@testset "Volume" begin
    @test V(1m) == (4π/3)m³
    @test V(r_sun) * ρ_sun / m_sun ≈ 0.999 atol = 0.001
    @test V(r_sun) isa Quantity
end

@testset "Sphere constructors" begin
    s1 = Sphere(p = [0]m, r = r_earth)
    s2 = Sphere([0]m, m_earth, r_earth)
    @test s1.m == s2.m
    @test s1.r == s2.r
    s3 = Sphere()
    @test s3.m ≈ s2.m
    @test s3.r ≈ s2.r
    s4 = earthlikeSphere_r_known(r = 10.E6m)
    s5 = earthlikeSphere_m_known(m = ρ_earth * V(10.E6m))
    @test s4.m ≈ s5.m
    @test s4.r ≈ s5.r
end

@testset "Position vectors" begin
    @test [1, 2, 3]m isa Pos
    @test !isa([1 2; 3 4]m, Pos)
    @test !isa([1m, 2, 3], Pos)
    @test !isa([1, sin, 3], Pos)
    @test !isa([[1, 2, 3, 4,5,6]], Pos)
    @test [1, 2, 3, 4,5,6]m isa Pos
    @test !isa([[1, 2, 3, 4,5,6]], Pos)
    @test !isa([1, 2, 3, 4,5,6]kg, Pos)
end


@testset "Gravity acceleration vector in system of bodies" begin
    earth = Sphere()
    moon = earthlikeSphere_m_known(p = [d_moon_earth, 0m, 0m], m = m_moon)
    sun = Sphere([-d_sun_earth, 0m, 0m], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    bodies = [sun, earth, moon]
    @test gvec(earth.p, bodies)[1] < 0m/s² # Sun stronger than moon in centre of earth
    @test gvec([r_earth/2, 0m, 0m], bodies)[1] < -5m/s² # Halfway to centre of earth, earth dominates 
    xgrid, ygrid = xygrids( xmin = -d_moon_earth/2., 
        xmax = d_moon_earth * 3/2, 
        ymin =-d_moon_earth/2, 
        ymax = d_moon_earth * 3/2, n = 100)
    gmagn = gfieldvals(nest(xgrid, ygrid, moon.r),  bodies)
    @test foldr(+, gmagn)  ≈  222.5649149449058m/s²
    @test gvec(earth, bodies) == gvec(earth.p, bodies)
end

@testset "Nested position vectors" begin
    earth = Sphere()
    moon = earthlikeSphere_m_known(p = [d_moon_earth, 0m, 0m], m = m_moon)
    sun = Sphere([-d_sun_earth, 0m, 0m], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    bodies = [sun, earth, moon]
    p1 = Pos([1., 2.]m)
    p2 = Pos([-1., -2.]m)
    vp = Npos([p1, p2])
    @test vp == [[1., 2.], [-1., -2.]]m
    mp = hcat(vp, vp)
    @test mp isa Npos
    pu1 = [1., 2.]
    pu2 = [-1., -2.]
    mup = hcat([pu1, pu2], [pu1, pu2])
    @test mp == mup*m
    @test nest(1m:1m:3m, 2m, 3m) == [[1, 2, 3], [2, 2, 3], [3, 2, 3]]m
    @test nest(1m,2m:1m:3m, 4m) == [[1, 2, 4], [1, 3, 4]]m
    @test nest(1m,2m, 3m:1m:4m) == [[1, 2, 3], [1, 2, 4]]m
    @test nest(1m:1m:2m, 1m:1m:2m) == [[1,1], [2,2]]m
    @test nest(1m,2m:1m:3m) == [[1,2], [1,3]]m
    @test nest(1m:1m:2m, 3m) == [[1,3], [2,3]]m
    @test nest(range(1.5m, 3.5m, length = 2)) == [[1.5], [3.5]]m
    @test gfieldvals(nest(1m,2m:1m:3m,4m), bodies)[1] ≈ 0.0059000661748m/s²
    @test nest([1m, 2m, 3m], [4, 5, 6]m, [7, 8, 9]m) ==  [[1, 4, 7], [2, 5, 8], [3, 6, 9]]m
    np = nest([1, 10]m, [2, 20]m, [3, 30]m, [4, 40]m, [5, 50]m, [6, 60]m)
    @test np ==  [[1, 2, 3, 4, 5, 6]m, [10, 20, 30, 40, 50, 60]m]
    @test slice_1(np) == [1, 10]m
    @test slice_2(np) == [2, 20]m
    @test slice_3(np) == [3, 30]m
    @test slice_4(np) == [4, 40]m
    @test slice_5(np) == [5, 50]m
    @test slice_6(np) == [6, 60]m
    @test gvec([[0m,0m,0m]], bodies)[1][1] ≈ -0.005898468352763m/s²
end

@testset "Bodies constructors" begin
    @test eltype(bodies_3d()) == Sphere
    @test eltype(bodies_2d()) == Sphere
    @test eltype(bodies_1d()) == Sphere
end

@testset "Pretty-printing" begin
    function sout(p)
        iob=IOBuffer()
        show(IOContext(iob, :color=>true), p)
        take!(iob)|>String
    end
    p = [1,2]m
    @test sout(p) == "[1.0, 2.0]\e[36mm\e[39m"
    function slongout(p)
        iob=IOBuffer()
        show(IOContext(iob, :color=>true), :"text/plain", p)
        take!(iob)|>String
    end
    @test slongout(p) == "Pos{Float64}([1.0, 2.0]\e[36mm\e[39m)"
    vp = [p, p]
    @test sout(vp) == "[[1.0, 2.0], [1.0, 2.0]]\e[36mm\e[39m"
    @test slongout(vp) == "2-element Npos{Pos{Float64}}:\n[[1.0, 2.0], [1.0, 2.0]]\e[36mm\e[39m"
    p1 = [1,1]m
    p2 = [2,2]m
    p3 = [3,3]m
    p4 = [4,4]m
    mp = hcat([p1, p2], [p3, p4])
    @test mp[2,1] == p2
    @test sout(mp) == "[[1.0, 1.0] [3.0, 3.0]; [2.0, 2.0] [4.0, 4.0]]\e[36mm\e[39m"
    @test slongout(mp) == "2x2 Npos{Pos{Float64}}:\n[[1.0, 1.0] [3.0, 3.0]; [2.0, 2.0] [4.0, 4.0]]\e[36mm\e[39m"
    @test sout(bodies_3d()) == "Sphere[Sphere(p = [-1.496e11, 0.0, 0.0]\e[36mm\e[39m, m = 1.989e30 kg, r = 6.95733e8 m), Sphere(p = [0.0, 0.0, 0.0]\e[36mm\e[39m, m = 5.97e24 kg, r = 6.3e6 m), Sphere(p = [3.844e8, 0.0, 0.0]\e[36mm\e[39m, m = 7.34767e22 kg, r = 1.45456e6 m)]"
    @test sout(bodies_2d()) == "Sphere[Sphere(p = [-1.496e11, 0.0]\e[36mm\e[39m, m = 1.989e30 kg, r = 6.95733e8 m), Sphere(p = [0.0, 0.0]\e[36mm\e[39m, m = 5.97e24 kg, r = 6.3e6 m), Sphere(p = [3.844e8, 0.0]\e[36mm\e[39m, m = 7.34767e22 kg, r = 1.45456e6 m)]"
    @test sout(bodies_1d()) == "Sphere[Sphere(p = [-1.496e11]\e[36mm\e[39m, m = 1.989e30 kg, r = 6.95733e8 m), Sphere(p = [0.0]\e[36mm\e[39m, m = 5.97e24 kg, r = 6.3e6 m), Sphere(p = [3.844e8]\e[36mm\e[39m, m = 7.34767e22 kg, r = 1.45456e6 m)]"
    @test slongout(bodies_3d()) == "3-element Array{Sphere,1}:\n  p = [-1.496e11, 0.0, 0.0]\e[36mm\e[39m m = 1.989e30 kg r = 6.95733e8 m\n  p = [0.0, 0.0, 0.0]\e[36mm\e[39m m = 5.97e24 kg r = 6.3e6 m\n  p = [3.844e8, 0.0, 0.0]\e[36mm\e[39m m = 7.34767e22 kg r = 1.45456e6 m"
    @test slongout(bodies_2d()) == "3-element Array{Sphere,1}:\n  p = [-1.496e11, 0.0]\e[36mm\e[39m m = 1.989e30 kg r = 6.95733e8 m\n  p = [0.0, 0.0]\e[36mm\e[39m m = 5.97e24 kg r = 6.3e6 m\n  p = [3.844e8, 0.0]\e[36mm\e[39m m = 7.34767e22 kg r = 1.45456e6 m"
    @test slongout(bodies_1d()) == "3-element Array{Sphere,1}:\n  p = [-1.496e11]\e[36mm\e[39m m = 1.989e30 kg r = 6.95733e8 m\n  p = [0.0]\e[36mm\e[39m m = 5.97e24 kg r = 6.3e6 m\n  p = [3.844e8]\e[36mm\e[39m m = 7.34767e22 kg r = 1.45456e6 m" 
end


nothing
