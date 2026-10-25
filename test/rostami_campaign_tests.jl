using Test, TOML, SHA
include(joinpath(@__DIR__,"..","scripts","make_rostami.jl"))
const ROSTAMI = RostamiExperiments

@testset "2017 QG analogue preserves experiment controls and immutable provenance" begin
    mktempdir() do tmp
        pilot=joinpath(@__DIR__,"..","configs","pilot.toml")
        out=joinpath(tmp,"campaign")
        result=ROSTAMI.make_campaign(pilot,out)
        @test result.count==32
        @test result.stages==Dict("validation"=>3,"controls"=>3,"formation"=>6,"competition"=>3,
            "qg_regime"=>3,"deep"=>4,"deformation"=>2,"core"=>3,"numerics"=>4,"boundary"=>1)
        @test_throws ArgumentError ROSTAMI.make_campaign(pilot,out)
        configurations=Dict{String,Any}()
        for stage in keys(result.stages), row in readlines(joinpath(out,stage,"members.tsv"))
            digest,name=split(row,'\t'); path=joinpath(out,stage,name)
            @test bytes2hex(sha256(read(path)))==digest
            c=TOML.parsefile(path); configurations[stage*"/"*name]=c
            @test !c["convection"]["enabled"]
            @test c["initial"]["profile"]=="rostami2017"
            @test c["physics"]["deep_profile"]=="matched_components"
            @test c["physics"]["deep_jet_speed_m_s"]==0
            @test !haskey(c["physics"],"deep_match_fraction")
            @test c["boundary"]["sponge_rate_s"]==0
            @test ROSTAMI.damping_bound(c)<=2
            stage=="validation" || @test c["grid"]["nx"]>=1024
        end
        baseline=configurations["formation/reference_bb2to8_s42.toml"]
        @test baseline["planet"]==Dict("radius_m"=>55e6,"rotation_rate_s"=>1.6e-4)
        @test baseline["initial"]["jet_speed_m_s"]==81.92
        @test baseline["initial"]["jet_radius_m"]==10.784e6
        @test baseline["initial"]["jet_width_m"]≈3.2e6/sqrt(3)
        @test baseline["initial"]["cyclone_radius_m"]≈(.42/1.3)^(1/1.3)*3.2e6
        @test baseline["numerics"]["duration_s"]==3.125e6
        @test baseline["numerics"]["duration_s"]*2baseline["planet"]["rotation_rate_s"]≈1000
        @test baseline["initial"]["m_min"]==2
        @test baseline["initial"]["m_max"]==8
        @test baseline["output"]["m_max"]==32
        @test baseline["experiment"]["purpose"]=="selection"
        # Deep-flow and deformation comparisons preserve the physical upper IC.
        for (key,c) in configurations
            if startswith(key,"deep/") || startswith(key,"deformation/")
                @test c["initial"]==baseline["initial"]
                @test c["planet"]==baseline["planet"]
                @test c["experiment"]["seed"]==baseline["experiment"]["seed"]
                @test c["grid"]==baseline["grid"]
                @test c["numerics"]==baseline["numerics"]
            end
        end
        fractions=Set((c["physics"]["deep_jet_fraction"],c["physics"]["deep_cyclone_fraction"]) for (key,c) in configurations if startswith(key,"deep/"))
        @test fractions==Set(((.5,0.0),(1.0,0.0),(0.0,1.0),(1.0,1.0)))
        competition=[c for (key,c) in configurations if startswith(key,"competition/")]
        @test Set(c["experiment"]["seed"] for c in competition)==Set(42:44)
        @test all(c->c["initial"]["m_min"]==1 && c["initial"]["m_max"]==32,competition)
        weak=[c for (key,c) in configurations if startswith(key,"qg_regime/")]
        @test Set(c["initial"]["cyclone_speed_m_s"] for c in weak)==Set(153.6 .* [.1,.25,.5])
        for c in weak
            a=copy(c["initial"]);a["cyclone_speed_m_s"]=153.6
            @test a==baseline["initial"]
            @test c["physics"]==baseline["physics"]
            @test c["numerics"]==baseline["numerics"]
        end
        small=configurations["core/epsilon250km_nx2048.toml"]
        reference=configurations["core/epsilon500km_nx2048.toml"]
        @test small["grid"]==reference["grid"]
        @test small["numerics"]==reference["numerics"]
        a=copy(small["initial"]); a["core_regularization_m"]=500e3
        @test a==reference["initial"]
        nu(c)=1/(c["numerics"]["cutoff_damping_time_s"]*(2pi/c["grid"]["length_m"]*floor(c["grid"]["nx"]/3))^8)
        @test nu(reference)≈nu(baseline)/256
        fixed=configurations["numerics/nx2048_fixed_nu.toml"]
        @test nu(fixed)≈nu(baseline)
        @test fixed["numerics"]["dt_s"]==25/32
        boundary=configurations["boundary/box100Mm.toml"]
        @test nu(boundary)≈nu(baseline)
        @test boundary["grid"]["length_m"]/boundary["grid"]["nx"]==baseline["grid"]["length_m"]/baseline["grid"]["nx"]
        @test boundary["initial"]==baseline["initial"]
        @test configurations["validation/reference_full_deep_smoke.toml"]["experiment"]["purpose"]=="smoke"
        for row in readlines(joinpath(out,"provenance.tsv"))
            digest,path=split(row,'\t')
            @test digest==bytes2hex(sha256(read(joinpath(ROSTAMI.ROOT,path))))
        end
        meta=TOML.parsefile(joinpath(out,"campaign.toml"))
        @test meta["paper_doi"]=="10.1016/j.icarus.2017.06.006"
        @test occursin("not 1000 inertial periods",meta["interpretation"]["duration"])
        @test occursin("declarative",meta["interpretation"]["regime"])
        @test !ispath(joinpath(out,"maintenance"))
    end
end
