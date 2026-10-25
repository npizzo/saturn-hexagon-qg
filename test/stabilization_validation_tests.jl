using Test
include(joinpath(@__DIR__,"..","scripts","validate_stabilization.jl"))

@testset "Stabilization CPU/GPU validation checks fail closed" begin
    v = StabilizationValidation
    fields0 = (; s=ones(3,3),psi=fill(2.0,3,3),u=fill(3.0,3,3),v=fill(4.0,3,3),totalpv=fill(5.0,3,3))
    @test all(e["relative_l2"] == 0 for e in values(v.field_errors(fields0,deepcopy(fields0))))
    changed = deepcopy(fields0); changed.psi[1] += 0.1
    @test_throws ErrorException v.field_errors(fields0,changed)
    changed = deepcopy(fields0); changed.s[1] = NaN
    @test_throws ErrorException v.field_errors(fields0,changed)
    off = (; time_s=400.0,step=20,events=0,next_event_s=Inf,rng_probe=UInt64[1,2])
    @test v.check_forcing_states(off,deepcopy(off),false) === nothing
    @test_throws ErrorException v.check_forcing_states(off,off,true)
    changed_clock = merge(off,(;step=21))
    @test_throws ErrorException v.check_forcing_states(off,changed_clock,false)
    changed_rng = merge(off,(;rng_probe=UInt64[1,3]))
    @test_throws ErrorException v.check_forcing_states(off,changed_rng,false)
    on = merge(off,(;events=4,next_event_s=450.0))
    @test v.check_forcing_states(on,on,true) === nothing
    @test_throws ErrorException v.check_forcing_states(on,on,false)
    @test_throws ErrorException v.options(["config.toml"])
    @test_throws ErrorException v.options(["config.toml","--device","cpu"])
    c = Dict("physics"=>Dict{String,Any}("deformation_radius_m"=>3.2e6,
        "deep_profile"=>"matched_components","deep_jet_fraction"=>1.0,"deep_cyclone_fraction"=>1.0))
    @test v.full_match_applicable(c)
    c["physics"]["deep_cyclone_fraction"] = 0.0
    @test !v.full_match_applicable(c)
    c["physics"]["deep_cyclone_fraction"] = 1.0
    c["physics"]["deformation_radius_m"] = Inf
    @test !v.full_match_applicable(c)
    mktempdir() do dir
        write(joinpath(dir,"COMPLETED"),"400.0\n")
        write(joinpath(dir,"CPU_GPU_VALIDATION_PASSED"),"stale\n")
        v.mark_failed(dir,ErrorException("intentional comparison failure"))
        @test !isfile(joinpath(dir,"COMPLETED"))
        @test !isfile(joinpath(dir,"CPU_GPU_VALIDATION_PASSED"))
        @test occursin("intentional comparison failure",read(joinpath(dir,"FAILED.txt"),String))
    end
end
