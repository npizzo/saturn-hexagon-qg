using Test, TOML, SHA
include(joinpath(@__DIR__,"..","scripts","make_experiments.jl"))
using .LowerLayerExperiments

@testset "Staged lower-layer experiment controls" begin
    pilot=joinpath(@__DIR__,"..","configs","unity_followup_01","maintenance_unforced.toml")
    input=TOML.parsefile(pilot)
    function stage_configs(root,stage)
        paths=String[]
        for line in readlines(joinpath(root,stage,"members.tsv"))
            digest,name=split(line,'\t')
            @test basename(name)==name
            path=joinpath(root,stage,name)
            @test bytes2hex(sha256(read(path)))==digest
            push!(paths,path)
        end
        TOML.parsefile.(paths)
    end
    mktempdir() do tmp
        out=joinpath(tmp,"campaign")
        result=LowerLayerExperiments.make_experiments(pilot,out)
        @test result.count==25
        @test result.counts==Dict("screen"=>6,"expanded"=>11,"numerics"=>8)
        screen=stage_configs(out,"screen")
        expanded=stage_configs(out,"expanded")
        numerical=stage_configs(out,"numerics")
        physical=vcat(screen,expanded)
        @test length(unique([(c["physics"]["deformation_radius_m"],c["physics"]["deep_jet_speed_m_s"]) for c in physical]))==17
        @test count(c->isinf(c["physics"]["deformation_radius_m"]),physical)==1
        @test all(c->c["initial"]==first(screen)["initial"],physical)
        @test all(c->c["initial"]["jet_speed_m_s"]==input["initial"]["jet_speed_m_s"],physical)
        @test all(c->c["initial"]["cyclone_speed_m_s"]==input["initial"]["cyclone_speed_m_s"],physical)
        @test all(c->c["planet"]==input["planet"],physical)
        @test all(c->!c["convection"]["enabled"] && c["initial"]["perturbation"]=="broadband",physical)
        @test all(c->c["initial"]["m_max"]==32 && c["output"]["m_max"]==32,physical)
        @test length(readlines(joinpath(out,"physics-audit.tsv")))==18
        matched=deepcopy(first(screen))
        matched["physics"]["deformation_radius_m"]=1e6
        matched["physics"]["deep_jet_speed_m_s"]=matched["initial"]["jet_speed_m_s"]
        matched["initial"]["cyclone_speed_m_s"]=0.0
        audit=LowerLayerExperiments.physical_regime(matched)
        @test audit.min_h_over_H==1.0 && audit.max_abs_hprime_over_H==0.0
        resting=deepcopy(matched)
        resting["physics"]["deep_jet_speed_m_s"]=0.0
        a1=LowerLayerExperiments.physical_regime(resting)
        resting["physics"]["deformation_radius_m"]=2e6
        a2=LowerLayerExperiments.physical_regime(resting)
        @test a1.max_abs_hprime_over_H≈4a2.max_abs_hprime_over_H
        @test (1-a1.min_h_over_H)≈4(1-a2.min_h_over_H)
        @test a2.rossby≈100/(0.00033*1e6)
        finer=LowerLayerExperiments.physical_regime(resting;samples=16385)
        @test abs(a2.max_abs_hprime_over_H-finer.max_abs_hprime_over_H)<1e-7
        @test LowerLayerExperiments.physical_regime(first(screen)).flag=="barotropic_limit_no_finite_layer_inference"
        small_ld=only(filter(c->c["physics"]["deformation_radius_m"]==500000 && c["physics"]["deep_jet_speed_m_s"]==0,physical))
        @test LowerLayerExperiments.physical_regime(small_ld).min_h_over_H<0
        @test count(c->c["physics"]["deep_jet_speed_m_s"]<0,physical)==4
        @test count(c->c["grid"]["nx"]==2048 && c["numerics"]["dt_s"]==25,numerical)==2
        @test count(c->c["initial"]["displacement_m"]==5000,numerical)==2
        @test count(c->c["numerics"]["cutoff_damping_time_s"]==1e6,numerical)==2
        meta=TOML.parsefile(joinpath(out,"campaign.toml"))
        @test !meta["exact_2024_family"]
        @test meta["frame"]["observed_angular_speed_model_s"]≈2pi/38363.01-input["planet"]["rotation_rate_s"]
        @test meta["frame"]["observed_speed_at_model_ring_m_s"]≈-18.264211862044
        @test_throws ArgumentError LowerLayerExperiments.make_experiments(pilot,out)
        @test_throws ArgumentError LowerLayerExperiments.make_experiments(pilot,joinpath(tmp,"bad");candidates=["bogus"])
        @test !ispath(joinpath(tmp,"bad"))
        @test_throws ArgumentError LowerLayerExperiments.make_experiments(pilot,joinpath(tmp,"bad2");seeds=[42,42])
        @test_throws ArgumentError LowerLayerExperiments.make_experiments(pilot,joinpath(tmp,"bad3");long_duration_s=1e6)
        promoted=joinpath(tmp,"promoted")
        r=LowerLayerExperiments.make_experiments(pilot,promoted;candidates=["ld4000_dp050"])
        @test r.count==35
        @test r.counts["maintenance"]==1 && r.counts["selection"]==3 && r.counts["convection"]==6
        maint=only(stage_configs(promoted,"maintenance"))
        selection=stage_configs(promoted,"selection")
        conv=stage_configs(promoted,"convection")
        @test all(c->c["physics"]["deformation_radius_m"]==4000000 && c["physics"]["deep_jet_speed_m_s"]==50,vcat([maint],selection,conv))
        @test maint["initial"]["perturbation"]=="mode" && maint["initial"]["displacement_m"]==100000
        @test Set(c["experiment"]["seed"] for c in selection)==Set(42:44)
        @test all(c->c["initial"]["perturbation"]=="broadband",selection)
        @test all(c->c["convection"]["enabled"],conv)
        @test !maint["convection"]["enabled"]
        @test count(c->c["convection"]["enabled"] && c["convection"]["bias_max"]==0,conv)==3
        @test count(c->c["convection"]["enabled"] && c["convection"]["bias_max"]==input["convection"]["bias_max"],conv)==3
        @test all(c->c["initial"]==maint["initial"],conv)
        @test all(c->c["numerics"]["duration_s"]==2e6,conv)
        @test all(c->c["output"]["diagnostic_interval_s"]==2000,vcat([maint],selection,conv))
        extra=LowerLayerExperiments.make_experiments(pilot,joinpath(tmp,"extra");candidates=["ld0500_dn050"])
        @test extra.counts["numerics"]==12
        for line in readlines(joinpath(out,"provenance.tsv"))
            digest,relative=split(line,'\t')
            @test !isabspath(relative)
            @test digest==bytes2hex(sha256(read(joinpath(@__DIR__,"..",relative))))
        end
        @test TOML.parsefile(pilot)==input
    end
end
