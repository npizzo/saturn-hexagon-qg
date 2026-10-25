using Test
if !isdefined(@__MODULE__, :PolarDiagnostics)
    include(joinpath(@__DIR__, "..", "src", "PolarDiagnostics.jl"))
end
using .PolarDiagnostics

function synthetic_jet(; amplitude=0.12, phase=0.47, n=256, radial_speed=0.0)
    x = collect(range(-4.0; step=8.0/n, length=n))
    y = copy(x)
    u, v = zeros(n,n), zeros(n,n)
    for j in eachindex(y), i in eachindex(x)
        r, θ = hypot(x[i],y[j]), atan(y[j],x[i])
        peak = 2.0 + amplitude*cos(6θ+phase)
        speed = 10exp(-0.5*((r-peak)/0.18)^2)
        if r > 0
            u[i,j] = -sin(θ)*speed + cos(θ)*radial_speed
            v[i,j] =  cos(θ)*speed + sin(θ)*radial_speed
        end
    end
    return u, v, x, y
end

@testset "Polar jet shape, amplitude, orientation, and velocities" begin
    phase, amplitude = 0.47, 0.12
    u,v,x,y = synthetic_jet(; amplitude, phase, radial_speed=0.7)
    z = zeros(size(u))
    d = polar_diagnostics(u,v,z,z,z,x,y; jet_radius_m=2.,jet_halfwidth_m=.6,nr=120)
    @test d.jet_detected
    @test d.coverage == 1
    @test d.meanjet ≈ 2.0 atol=0.002
    @test d.mode_amplitudes[6] ≈ amplitude atol=0.004
    @test angle(d.mode_coefficients[6]) ≈ phase atol=0.02
    @test maximum(d.mode_amplitudes[[1,2,3,4,5,7,8,9,10,11,12]]) < 0.004
    @test maximum(abs.(d.mean_ur .- .7)) < 0.002
    @test maximum(d.mean_utheta) > 7
    @test d.mean_peak_speed > 9.8
    @test all(d.radial_coverage .== 1)
    @test maximum(abs.(d.mean_pv)) == 0

    # Reversing angular phase reverses the complex coefficient's imaginary part.
    um,vm,_,_ = synthetic_jet(; amplitude, phase=-phase)
    dm = polar_diagnostics(um,vm,z,z,z,x,y; jet_radius_m=2.,jet_halfwidth_m=.6,nr=120)
    @test angle(dm.mode_coefficients[6]) ≈ -phase atol=0.02
    # Wrong velocity sign must not be presented as an eastward Saturn jet.
    dw = polar_diagnostics(-u,-v,z,z,z,x,y; jet_radius_m=2.,jet_halfwidth_m=.6,nr=120)
    @test !dw.jet_detected
end

@testset "Axisymmetric and degenerate snapshots" begin
    u,v,x,y = synthetic_jet(; amplitude=0.0)
    z = zeros(size(u))
    d = polar_diagnostics(u,v,z,z,z,x,y; jet_radius_m=2.,jet_halfwidth_m=.6,nr=120)
    @test d.jet_detected
    @test d.mode_amplitudes[6] < 1e-10
    @test maximum(d.mode_amplitudes) < 0.002
    empty = polar_diagnostics(z,z,z,z,z,x,y; jet_radius_m=2.,jet_halfwidth_m=.6)
    @test !empty.jet_detected
    @test empty.coverage == 0
    @test isnan(empty.meanjet)
    @test all(isnan, empty.mode_amplitudes)
    @test all(isnan, empty.jet_radius_theta)
    # A constant azimuthal wind has no distinguished jet radius. Bilinear
    # interpolation of its Cartesian components must not manufacture one.
    flat_u=[-10sin(atan(yy,xx)) for xx in x, yy in y]
    flat_v=[10cos(atan(yy,xx)) for xx in x, yy in y]
    flat=polar_diagnostics(flat_u,flat_v,z,z,z,x,y;jet_radius_m=2.,jet_halfwidth_m=.6)
    @test !flat.jet_detected
    @test flat.coverage == 0
    @test all(isnan, flat.mode_amplitudes)
    umasked = copy(u)
    umasked[180:220,120:136] .= NaN
    masked = polar_diagnostics(umasked,v,z,z,z,x,y; jet_radius_m=2.,jet_halfwidth_m=.6)
    @test 0 < masked.coverage < 1
    @test all(isnan, masked.mode_amplitudes)
    @test_throws ArgumentError polar_diagnostics(u,v,z,z,z,x,y; jet_radius_m=3.9,jet_halfwidth_m=.6)
end

@testset "Radial anomaly and periodic bilinear interpolation" begin
    n=64
    x=collect(range(-2.;step=4/n,length=n)); y=copy(x)
    constant=fill(5.,n,n)
    @test radial_anomaly(constant,x,y) == zeros(n,n)
    radius=[hypot(xx,yy) for xx in x, yy in y]
    radial=radius.^2
    residual=radial_anomaly(radial,x,y;nbins=64)
    @test maximum(abs.(residual[radius .< 2])) < .002
    # Constant extension of the last populated bin has O(dr) edge error.
    @test maximum(abs.(residual)) < 2maximum(radius)^2/64
    # Every point's residual is computed against the same radial interpolant.
    @test residual[20,30] ≈ residual[30,20] atol=1e-12
    field=[sin(2π*xx/4)+cos(2π*yy/4) for xx in x, yy in y]
    grid=PolarDiagnostics._grid(x,y,size(field))
    @test PolarDiagnostics._sample(field,2.01,.2,grid) ≈ PolarDiagnostics._sample(field,-1.99,.2,grid) atol=1e-13
    @test PolarDiagnostics._sample(field,.1,2.02,grid) ≈ PolarDiagnostics._sample(field,.1,-1.98,grid) atol=1e-13
    @test_throws DimensionMismatch radial_anomaly(zeros(2,3),x,y)
end

@testset "Local cubic PV anomaly and spatial convergence" begin
    # A nonuniform axisymmetric Gaussian ring should have zero anomaly.
    # Keep the physical field and sample locations fixed while refining dx.
    # Quadratic background is reproduced by the cubic interpolant; this test
    # therefore exposes its residual on the smooth non-polynomial PV ring.
    errors = Float64[]
    for n in (64,128,256)
        x=collect(range(-4.;step=8/n,length=n))
        field=[3+.2*(xx^2+yy^2)+exp(-.5*((hypot(xx,yy)-2)/.35)^2)
               for xx in x, yy in x]
        theta=(collect(0:47).+.31).*(2π/48)
        residual=[local_pv_anomaly(field,x,x,r*cos(t),r*sin(t))
                  for r in range(1.3,2.7;length=13), t in theta]
        push!(errors,sqrt(sum(abs2,residual)/length(residual)))
    end
    @test errors[1] < 2e-4
    @test errors[3] < 1e-6
    @test all(log2.(errors[1:2]./errors[2:3]) .> 3.5)

    n=256; x=collect(range(-4.;step=8/n,length=n))
    phase=.37
    field=[3+.2*(xx^2+yy^2)+exp(-.5*((hypot(xx,yy)-2)/.35)^2)*
           (2+.2*cos(6atan(yy,xx)+phase)) for xx in x, yy in x]
    for theta in (.19,.53,1.28)
        @test local_pv_anomaly(field,x,x,2cos(theta),2sin(theta)) ≈ .2*cos(6theta+phase) atol=3e-6
        @test local_pv_anomaly(field,x,x,2cos(theta),2sin(theta);ntheta=2048) ≈ .2*cos(6theta+phase) atol=3e-6
    end
    # The mean and local value use precisely the same interpolant, including
    # constant backgrounds. Missing circle sectors are not treated as zeros.
    @test local_pv_anomaly(fill(5.,n,n),x,x,1.17,.49) == 0
    @test local_pv_anomaly(field,x,x,0.,0.) == 0
    masked=copy(field); masked[192:196,128:130].=NaN
    @test isnan(local_pv_anomaly(masked,x,x,0.,2.))
    @test_throws ArgumentError local_pv_anomaly(field,x,x,NaN,0.)
    @test_throws ArgumentError local_pv_anomaly(field,x,x,1.,1.;ntheta=7)
    @test_throws ArgumentError local_pv_anomaly(field,x,x,1.,1.;ntheta=64.0)
    @test_throws DimensionMismatch local_pv_anomaly(zeros(3,4),x,x,0.,0.)
    @test_throws ArgumentError local_pv_anomaly(zeros(3,3),[0.,1.,2.],[0.,1.,2.],0.,0.)
end

@testset "Periodic tensor cubic interpolant" begin
    errors=Float64[]
    for n in (16,32,64)
        x=collect(range(-2.;step=4/n,length=n))
        f=[sin(2π*xx/4)+cos(2π*yy/4) for xx in x, yy in x]
        g=PolarDiagnostics._grid(x,x,size(f))
        positions=[(1.91,.27),(-.73,1.94),(.33,-.48)]
        sampled=[PolarDiagnostics._sample_cubic(f,xx,yy,g) for (xx,yy) in positions]
        exact=[sin(2π*xx/4)+cos(2π*yy/4) for (xx,yy) in positions]
        push!(errors,maximum(abs.(sampled.-exact)))
        @test PolarDiagnostics._sample_cubic(f,2.01,.2,g) ≈ PolarDiagnostics._sample_cubic(f,-1.99,.2,g) atol=1e-13
        @test PolarDiagnostics._sample_cubic(f,.1,2.02,g) ≈ PolarDiagnostics._sample_cubic(f,.1,-1.98,g) atol=1e-13
        @test PolarDiagnostics._sample_cubic(fill(5.,n,n),.13,-.79,g) == 5.
        @test PolarDiagnostics._sample_cubic(f,x[4],x[7],g) == f[4,7]
    end
    @test errors[2] < errors[1]/8
    @test errors[3] < errors[2]/8
end

@testset "Explicit wind CSV contract" begin
    mktempdir() do dir
        path=joinpath(dir,"profile.csv")
        write(path,"# converted using documented geometry and frame\nradius_m,u_theta_m_s\n0,0\n1000,12\n2000,2\n")
        p=load_wind_profile(path;required_radius_range=(0.,1900.))
        @test p.radius_m == [0.,1000.,2000.]
        @test p.u_theta_m_s == [0.,12.,2.]
        @test p.radius_min_m == 0
        @test p.radius_max_m == 2000
        @test_throws ArgumentError load_wind_profile(path;required_radius_range=(0.,3000.))
        for content in (
            "latitude,u\n75,100\n76,90\n",
            "radius_m,u_theta_m_s\n0,0\n0,2\n",
            "radius_m,u_theta_m_s\n1000,1\n0,2\n",
            "radius_m,u_theta_m_s\n-1,1\n1000,2\n",
            "radius_m,u_theta_m_s\n0,1\n1000,NaN\n",
            "radius_m,u_theta_m_s\n0,1\n1000,2,3\n",
            "radius_m,u_theta_m_s\n0,1\n")
            write(path,content)
            @test_throws ArgumentError load_wind_profile(path)
        end
    end
end
