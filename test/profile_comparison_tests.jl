using Test
if !isdefined(@__MODULE__, :ProfileComparison)
    include(joinpath(@__DIR__,"..","scripts","compare_profile.jl"))
end
using .ProfileComparison

@testset "Profile comparison without extrapolation" begin
    mktempdir() do dir
        observation=joinpath(dir,"observed.csv")
        simulation=joinpath(dir,"simulated.csv")
        write(observation,"radius_m,u_theta_m_s\n0,2\n1,4\n2,6\n3,8\n")
        write(simulation,"radius_m,u_theta_m_s\n0,3\n2,7\n")
        result=compare_profiles(observation,simulation)
        @test result.n_observed == 4
        @test result.n_compared == 3
        @test result.coverage == .75
        @test result.radius_m == [0.,1.,2.]
        @test result.simulated_u_theta_m_s == [3.,5.,7.]
        @test result.rmse_m_s ≈ 1
        @test result.bias_m_s ≈ 1
        buffer=IOBuffer()
        write_summary(buffer,result)
        summary=String(take!(buffer))
        @test startswith(summary,"n_observed,n_compared,coverage,rmse_m_s,bias_m_s,")
        @test occursin("4,3,0.75,1.0,1.0,0.0,2.0",summary)
        output=joinpath(dir,"nested","summary.csv")
        @test ProfileComparison.main([observation,simulation,output]) == 0
        @test read(output,String) == summary
        @test_throws ArgumentError ProfileComparison.main([observation,simulation,observation])
        write(simulation,"radius_m,u_theta_m_s\n4,1\n5,2\n")
        @test_throws ArgumentError compare_profiles(observation,simulation)
        write(simulation,"radius_m,u_theta_m_s\n0,NaN\n5,2\n")
        @test_throws ArgumentError compare_profiles(observation,simulation)
    end
end
