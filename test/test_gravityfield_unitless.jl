using Test
using Gravityfield_unitless

let
    @test V(1.) == 4π / 3
    @test ρ(m = 4π / 3, r = 1) == 1.
    @test ρ(r = 1, m = 4π / 3) == 1.
    @test ρ_earth |> typeof == Float64
end

let
    "Example position"
    p1 = [0.,0.,0.]
    "Example position"
    p2 = [0.,1.,0.]
    "Example mass"
    m1 = 316227. #kg
    "Example mass"
    m2 = 129099. #kg
end
@testset "Gravity acceleration" begin
    @test g(1, 0) == 0.
    @test g(1., 1.) == γ
    @test g(m_earth, r_earth) ≈ 10.0391664
    ω_earth_sun = 2π / yr
    a_earth_sun = ω_earth_sun^2 * d_sun_earth
    @test round(g(m_sun, d_sun_earth), digits = 4) == round(a_earth_sun, digits = 4)
end

@testset "Force vector on 1 due to 2" begin
    p1 = [0.,0.,0.] #m
    p2 = [0.,1.,0.] #m
    m1 = 316227. #kg
    m2 = 129099. #kg
    f1_2 = m1* uvec(p1, p2) * g(m2, len(p2-p1))
    @test f1_2 |> length == 3
    @test f1_2 |> typeof == Vector{Float64}
    @test f1_2 ≈ [0.0, 2.7247474, 0.]
end
@testset "Euclidean vector length" begin
    @test len(1) == 1.
    @test len([1]) == 1.
    @test len([2,0]) == 2.
    @test len([1,1]) == sqrt(2)
end
@testset "Vector between bodies or points" begin
    p1 = [0.,0.,0.] #m
    p2 = [0.,1.,0.] #m
    s = Sphere(p1, m_earth, r_earth)
    @test rvec(s,s) == [0., 0., 0.]
    @test len(s,s) == 0.
    @test uvec(s,s) == [0., 0., 0.]
    @test rvec(p2,s) == [0., -1.0, 0.]
    @test len(p2,s) == 1.
    @test uvec(p2,s) == [0., -1.0, 0.]
    s2 = Sphere(p2, m_earth, r_earth)
    @test rvec(s,s2) == [0., 1.0, 0.]
    @test len(s,s2) == 1
    @test uvec(s,s2) == [0., 1.0, 0.]
    s3 = Sphere([-d_sun_earth], 1, 1)
    s4 = Sphere([0.], 1, 1)
    @test rvec(s3,s4) == [d_sun_earth]
    @test len(s4,s3) == d_sun_earth
    @test uvec(s3,s4) == [1]
    @test uvec(s4,s3) == [-1]
end

@testset "Gravity acceleration vector btw. bodies or point and body" begin
    p1 = [0.,0.,0.] #m
    s = Sphere(p1, 1., 1.)
    @test gvec(s,s) == [0., 0., 0.]
    @test gvec([-1, 0., 0.], s) == [γ, 0., 0.]
    s2 = Sphere([0.0], m_earth, r_earth)
    @test gvec([r_earth], s2) == [-10.039166439909296]
    @test gvec([r_earth/2], s2) == [-10.039166439909296]/2
end

@testset "Volume" begin
    @test V(1) == 4π/3
    @test round(V(r_sun) * ρ_sun /m_sun, digits=3) == 0.999
    @test typeof(V(r_sun)) == Float64
end

@testset "Sphere constructors" begin
    s1 = Sphere(p = [0], r = r_earth)
    s2 = Sphere([0], m_earth, r_earth)
    @test s1.m == s2.m
    @test s1.r == s2.r
    s3 = Sphere()
    @test s3.m ≈ s2.m
    @test s3.r ≈ s2.r
    s4 = earthlikeSphere_r_known(r = 10.E6)
    s5 = earthlikeSphere_m_known(m = ρ_earth * V(10.E6))
    @test s4.m ≈ s5.m
    @test s4.r ≈ s5.r
end

@testset "Position vectors" begin
    @test [1, 2, 3] isa Pos
    @test !isa([1 2; 3 4], Pos)
    @test !isa([1, "2", 3], Pos)
    @test !isa([1, sin, 3], Pos)
    @test !isa([[1, 2, 3, 4,5,6]], Pos)
    @test Pos([1, 2, 3, 4,5,6]) isa Pos
    @test !isa([[1, 2, 3, 4,5,6]], Pos)
end


@testset "Gravity acceleration vector in system of bodies" begin
    earth = Sphere()
    moon = earthlikeSphere_m_known(p = [d_moon_earth, 0., 0.], m = m_moon)
    sun = Sphere([-d_sun_earth, 0., 0.], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    bodies = [sun, earth, moon]
    @test gvec(earth.p, bodies)[1] < 0 # Sun stronger than moon in centre of earth
    @test gvec([r_earth/2, 0., 0.], bodies)[1] < -5 # Halfway to centre of earth, earth dominates 
    xgrid, ygrid = xygrids( xmin = -d_moon_earth/2., 
    xmax = d_moon_earth * 3/2, 
    ymin =-d_moon_earth/2, 
    ymax = d_moon_earth * 3/2, n = 100)
    gmagn = gfieldvals(nest(xgrid, ygrid, moon.r),  bodies)
    @test foldr(+, gmagn)  ≈  222.5649149449058
    @test gvec(earth, bodies) == gvec(earth.p, bodies)
end

@testset "Nested position vectors" begin
    earth = Sphere()
    moon = earthlikeSphere_m_known(p = [d_moon_earth, 0., 0.], m = m_moon)
    sun = Sphere([-d_sun_earth, 0., 0.], m_sun, (3m_sun / (4π * ρ_sun))^(1/3))
    bodies = [sun, earth, moon]
    p1 = Pos([1., 2.])
    p2 = Pos([-1., -2.])
    @test Npos([p1, p2]) == [[1., 2.], [-1., -2.]]
    @test nest(1:3, 2, 3) == [[1, 2, 3], [2, 2, 3], [3, 2, 3]]
    @test nest(1,2:3, 4) == [[1, 2, 4], [1, 3, 4]]
    @test nest(1,2, 3:4) == [[1, 2, 3], [1, 2, 4]]
    @test nest(1:2, 1:2) == [[1,1], [2,2]]
    @test nest(1,2:3) == [[1,2], [1,3]]
    @test nest(1:2,3) == [[1,3], [2,3]] 
    @test nest(range(1.5, 3.5, length = 2)) == [[1.5], [3.5]]
    @test gfieldvals(nest(1,2:3,4), bodies)[1] ≈ 0.0059000661748
    @test nest([1, 2, 3], [4, 5, 6], [7, 8, 9]) ==  [[1, 4, 7], [2, 5, 8], [3, 6, 9]]
    np = nest([1, 10], [2, 20], [3, 30], [4, 40], [5, 50], [6, 60])
    @test np ==  [[1, 2, 3, 4, 5, 6], [10, 20, 30, 40, 50, 60]]
    @test slice_1(np) == [1, 10]
    @test slice_2(np) == [2, 20]
    @test slice_3(np) == [3, 30]
    @test slice_4(np) == [4, 40]
    @test slice_5(np) == [5, 50]
    @test slice_6(np) == [6, 60]
    @test gvec([[0.,0.,0]], bodies)[1][1] ≈ -0.005898468352763
end
nothing
