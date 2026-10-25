using Test, TOML, SHA
if !isdefined(@__MODULE__, :SuiteGenerator)
    include(joinpath(@__DIR__,"..","scripts","make_suite.jl"))
end
using .SuiteGenerator

@testset "Paired suites preserve hypotheses and scheduler hashes" begin
    pilot=joinpath(@__DIR__,"..","configs","pilot.toml")
    mktempdir() do dir
        result=make_suite(pilot,joinpath(dir,"suite");seeds=[42],include_emergence=true)
        @test result.count == 24
        @test all(isfile,result.paths)
        @test result.manifest_sha256 == bytes2hex(sha256(read(result.manifest)))
        lines=readlines(result.manifest)
        @test length(lines) == 24
        for line in lines
            digest,path=split(line,'\t')
            @test isabspath(path)
            @test digest == bytes2hex(sha256(read(path)))
        end
        configs=TOML.parsefile.(result.paths)
        @test count(c -> c["experiment"]["purpose"] == "maintenance",configs) == 9
        @test count(c -> c["experiment"]["purpose"] == "selection",configs) == 9
        emergence=filter(c -> c["experiment"]["purpose"] == "emergence",configs)
        @test length(emergence) == 6
        @test all(c -> c["initial"]["jet_speed_m_s"] == 0 && c["initial"]["cyclone_speed_m_s"] == 0 && c["initial"]["perturbation"] == "none" && c["initial"]["displacement_m"] == 0 && c["convection"]["enabled"], emergence)
        @test all(c -> c["physics"]["deep_jet_speed_m_s"] == 0,filter(c -> isinf(c["physics"]["deformation_radius_m"]),configs))
        @test all(c -> c["initial"]["perturbation"] == "broadband",filter(c -> c["experiment"]["purpose"] == "selection",configs))
        @test all(c -> c["initial"]["mode"] == 6 && c["initial"]["perturbation"] == "mode",filter(c -> c["experiment"]["purpose"] == "maintenance",configs))
        @test_throws ArgumentError make_suite(pilot,result.directory)
        @test_throws ArgumentError make_suite(pilot,joinpath(dir,"bad");seeds=[42,42])
        @test !ispath(joinpath(dir,"bad"))
        default=make_suite(pilot,joinpath(dir,"default"))
        @test default.count == 54
        @test !any(c -> c["experiment"]["purpose"] == "emergence",TOML.parsefile.(default.paths))
        @test readlines(joinpath(result.directory,"config_paths.txt")) == result.paths
    end
end
