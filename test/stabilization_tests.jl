using TOML, SHA
include(joinpath(@__DIR__,"..","scripts","make_stabilization.jl"))
const STAB = StabilizationExperiments

@testset "Unforced whole-profile stabilization protocol" begin
    mktempdir() do tmp
        pilot=joinpath(@__DIR__,"..","configs","pilot.toml")
        out=joinpath(tmp,"campaign")
        result=STAB.make_campaign(pilot,out)
        @test result.count==39
        @test result.stages==Dict("validation"=>3,"controls"=>5,"screen"=>9,"transition"=>4,"width"=>6,"numerics"=>10,"boundary"=>2)
        @test_throws ArgumentError STAB.make_campaign(pilot,out)
        configs=Dict{String,Any}()
        for stage in keys(result.stages)
            for row in readlines(joinpath(out,stage,"members.tsv"))
                hash,name=split(row,'\t')
                path=joinpath(out,stage,name)
                @test bytes2hex(sha256(read(path)))==hash
                c=TOML.parsefile(path);configs[stage*"/"*name]=c
                @test !c["convection"]["enabled"]
                @test c["physics"]["deep_profile"]=="matched_upper"
                @test c["physics"]["deep_jet_speed_m_s"]==0
                @test STAB.damping_bound(c)<=2
                @test c["boundary"]["sponge_rate_s"]==0
            end
        end
        baseline=configs["screen/ld3200_a100_w1000__bb10km_s42.toml"]
        for key in ("screen/ld2000_a000_w1000__bb10km_s42.toml", "screen/ld4000_a050_w1000__bb10km_s42.toml")
            @test configs[key]["initial"]==baseline["initial"]
            @test configs[key]["planet"]==baseline["planet"]
        end
        fixed=configs["numerics/ld3200_a100_w1000__nx2048_fixed_nu.toml"]
        nu(c)=1/(c["numerics"]["cutoff_damping_time_s"]*(2pi/c["grid"]["length_m"]*floor(c["grid"]["nx"]/3))^8)
        @test nu(fixed)≈nu(baseline)
        boundary=configs["boundary/ld3200_a100_w1000__box75Mm.toml"]
        @test nu(boundary)≈nu(baseline)
        @test boundary["grid"]["length_m"]/boundary["grid"]["nx"]==baseline["grid"]["length_m"]/baseline["grid"]["nx"]
        matched=STAB.physical_audit(baseline)
        @test matched.base_max==0
        @test matched.base_min==1
        @test matched.perturb_bound>0
        rest=STAB.physical_audit(configs["screen/ld3200_a000_w1000__bb10km_s42.toml"])
        @test rest.base_max>0
        bt=STAB.physical_audit(configs["controls/bt_a000_w1000__axisymmetric.toml"])
        @test isnan(bt.base_max)
        for row in readlines(joinpath(out,"provenance.tsv"))
            hash,path=split(row,'\t')
            @test hash==bytes2hex(sha256(read(joinpath(STAB.ROOT,path))))
        end
        meta=TOML.parsefile(joinpath(out,"campaign.toml"))
        @test occursin("ambient PV",meta["published_states"]["pizzo_salmon_2024"])
        @test meta["frame"]["rotation_rate_s"]==baseline["planet"]["rotation_rate_s"]
        @test_throws ArgumentError STAB.make_campaign(pilot,joinpath(tmp,"bad");candidates=["unknown"])
        @test !ispath(joinpath(tmp,"bad"))
        candidate=STAB.make_campaign(pilot,joinpath(tmp,"candidate");candidates=["ld500_a100_w1000"],seeds=42:43)
        @test candidate.count==42
        m=TOML.parsefile(joinpath(tmp,"candidate","maintenance","ld500_a100_w1000__m6_100km.toml"))
        @test m["initial"]["mode"]==6
        @test !m["convection"]["enabled"]
        @test m["numerics"]["duration_s"]==5e6
        @test_throws ArgumentError STAB.make_campaign(pilot,joinpath(tmp,"too_short");long_duration_s=1e6)
    end
end
