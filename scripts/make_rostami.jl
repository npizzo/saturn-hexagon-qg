#!/usr/bin/env julia
"""Staged QG analogue of Rostami, Zeitlin & Spiga (2017). Never submits jobs."""
module RostamiExperiments
using TOML, SHA
include(joinpath(@__DIR__, "make_experiments.jl"))
const Legacy = LowerLayerExperiments
const Validator = Legacy.Validator
const ROOT = dirname(@__DIR__)
const PAPER_LENGTH = 3.2e6
const PAPER_F0 = 3.2e-4
const PAPER_DURATION = 3.125e6 # 1000/f0, written exactly to end on a timestep boundary.
const CORE_RADIUS = (.42/1.3)^(1/1.3)*PAPER_LENGTH

"Worst-mode explicit hyperviscosity bound used by the solver, without sponge."
function damping_bound(c)
    n=c["grid"]["nx"]; p=c["numerics"]["hyper_order"]
    c["numerics"]["dt_s"]/c["numerics"]["cutoff_damping_time_s"] *
        (2*(n/2)^2/floor(n/3)^2)^p
end

function base_config(pilot)
    c=Legacy.base_config(pilot)
    c["planet"]["radius_m"]=55e6
    c["planet"]["rotation_rate_s"]=PAPER_F0/2
    c["grid"]["length_m"]=80e6
    c["physics"]["deformation_radius_m"]=PAPER_LENGTH
    c["physics"]["deep_profile"]="matched_components"
    c["physics"]["deep_jet_speed_m_s"]=0.0
    c["physics"]["deep_jet_fraction"]=0.0
    c["physics"]["deep_cyclone_fraction"]=0.0
    delete!(c["physics"],"deep_match_fraction")
    c["initial"]["profile"]="rostami2017"
    c["initial"]["core_regularization_m"]=500e3
    c["initial"]["jet_speed_m_s"]=81.92
    c["initial"]["jet_radius_m"]=3.37PAPER_LENGTH
    c["initial"]["jet_width_m"]=PAPER_LENGTH/sqrt(3)
    c["initial"]["cyclone_speed_m_s"]=153.6
    c["initial"]["cyclone_radius_m"]=CORE_RADIUS
    c["initial"]["m_min"]=2
    c["initial"]["m_max"]=8
    # The diagnostic band remains wide even when the initial seeds are 2:8.
    c["output"]["m_max"]=32
    c["output"]["jet_search_halfwidth_m"]=4e6
    c["boundary"]["inner_radius_m"]=28e6
    c["boundary"]["outer_radius_m"]=36e6
    c["numerics"]["dt_s"]=25.0
    c["numerics"]["duration_s"]=PAPER_DURATION
    c
end

function make_campaign(pilot_path,new_directory)
    out=abspath(new_directory)
    (ispath(out)||islink(out)) && throw(ArgumentError("Choose a NEW campaign directory; existing campaigns are frozen"))
    any(ch->ch in ('\t','\r','\n'),out) && throw(ArgumentError("Invalid output path"))
    base=base_config(Validator.load_config(pilot_path))
    members=NamedTuple[]
    function add(stage,id;core=153.6,ld=PAPER_LENGTH,jet_fraction=0.0,core_fraction=0.0,
                 family="broadband",m_min=2,m_max=8,seed=42,displacement=1e4,
                 epsilon=500e3,nx=1024,dt=25.0,duration=PAPER_DURATION,tau=1e5,
                 box=80e6,inner=28e6,outer=36e6,role="")
        c=deepcopy(base)
        c["experiment"]["name"]=stage*"__"*id
        c["experiment"]["purpose"]=stage=="validation" ? "smoke" : "selection"
        c["experiment"]["seed"]=seed
        c["physics"]["deformation_radius_m"]=ld
        c["physics"]["deep_jet_fraction"]=jet_fraction
        c["physics"]["deep_cyclone_fraction"]=core_fraction
        c["initial"]["cyclone_speed_m_s"]=core
        c["initial"]["core_regularization_m"]=epsilon
        c["initial"]["perturbation"]=family
        c["initial"]["m_min"]=m_min
        c["initial"]["m_max"]=m_max
        c["initial"]["displacement_m"]=displacement
        c["grid"]["nx"]=nx
        c["grid"]["length_m"]=box
        c["boundary"]["inner_radius_m"]=inner
        c["boundary"]["outer_radius_m"]=outer
        c["numerics"]["dt_s"]=dt
        c["numerics"]["duration_s"]=duration
        c["numerics"]["cutoff_damping_time_s"]=tau
        c["output"]["diagnostic_interval_s"]=min(2000.0,duration/4)
        c["output"]["snapshot_interval_s"]=min(50_000.0,duration/4)
        c["output"]["checkpoint_interval_s"]=min(200_000.0,duration/2)
        Validator.validate_config(c)
        damping_bound(c)<=2 || throw(ArgumentError("Explicit damping bound exceeds 2 for $id"))
        push!(members,(;stage,id,role,cfg=c))
    end
    add("validation","jet_only_smoke";core=0.0,nx=512,dt=10.0,duration=400.0,
        role="Software CPU/GPU check, resting deep flow; not spatial convergence")
    add("validation","reference_half_deep_jet_smoke";jet_fraction=.5,nx=512,dt=10.0,duration=400.0,
        role="Software CPU/GPU check with jet-only deep matching; core epsilon spans 3.2 cells")
    add("validation","reference_full_deep_smoke";jet_fraction=1.0,core_fraction=1.0,nx=512,dt=10.0,duration=400.0,
        role="Software CPU/GPU check with both deep components matched")
    for (label,core) in (("jet_only",0.0),("reference",153.6),("strong_core",204.8))
        add("controls",label*"_axisymmetric";core,family="none",displacement=0.0,duration=1e5,
            role="Axisymmetric QG equilibrium and polar-core numerical drift check")
        add("formation",label*"_m6";core,family="mode",
            role="Seeded sixfold formation from circular flow; imposed symmetry is not mode selection")
        add("formation",label*"_bb2to8_s42";core,
            role="Competing modes 2:8; test the central cyclone effect on formation and persistence")
    end
    for seed in 42:44
        add("competition","reference_bb1to32_s$seed";m_min=1,m_max=32,seed,
            role="Wider mode band and independent random phases; inspect mode-band and seed sensitivity")
    end
    for fraction in (.1,.25,.5)
        label=lpad(string(round(Int,100fraction)),3,'0')
        add("qg_regime","core_fraction$(label)_bb2to8_s42";core=153.6fraction,
            role="Reduce upper cyclone amplitude at fixed radius and regularization; explicit departure from paper, not QG validity by construction")
    end
    for (label,jet_fraction,core_fraction) in (("half_jet",.5,0.0),("full_jet",1.0,0.0),("full_core",0.0,1.0),("full_both",1.0,1.0))
        add("deep",label*"_bb2to8_s42";jet_fraction,core_fraction,
            role="Same upper winds and seed as formation/reference_bb2to8_s42; modify prescribed deep components")
    end
    for ld in (2.4e6,4e6)
        add("deformation","ld$(round(Int,ld/1000))_bb2to8_s42";ld,
            role="Change deformation radius at fixed dimensional upper winds and resting deep flow")
    end
    add("core","epsilon250km_nx2048";epsilon=250e3,nx=2048,dt=12.5,
        role="Inner-core model sensitivity; compare to epsilon500km_nx2048 at identical resolution")
    add("core","epsilon500km_nx2048";nx=2048,dt=12.5,
        role="Reference core on doubled mesh; fixed cutoff damping decreases physical hyperviscosity by 256")
    add("core","epsilon1000km_nx1024";epsilon=1e6,
        role="Smoother inner-core model at baseline resolution; not convergence toward the singular paper fit")
    add("numerics","half_dt";dt=12.5,
        role="Timestep refinement at unchanged spatial grid and physical hyperviscosity")
    add("numerics","nx2048_fixed_nu";nx=2048,dt=25/32,tau=1e5/256,
        role="Doubled mesh at fixed physical hyperviscosity; expensive explicit-damping timestep")
    add("numerics","weaker_damping";tau=1e6,
        role="Tenfold weaker hyperviscosity at unchanged spatial grid")
    add("numerics","displacement5km";displacement=5e3,
        role="Half displacement parameter; check early linear growth and diagnostic floor")
    k_old=2pi/80e6*floor(1024/3)
    k_new=2pi/100e6*floor(1280/3)
    add("boundary","box100Mm";nx=1280,box=100e6,inner=35e6,outer=45e6,
        tau=1e5*(k_old/k_new)^8,
        role="Larger box and outward taper at unchanged dx and physical hyperviscosity; inspect gauge sensitivity")
    allunique((m.stage,m.id) for m in members) || error("Duplicate experiment")
    stages=unique(m.stage for m in members)
    mkpath(dirname(out));staging=mktempdir(dirname(out);prefix=".rostami-")
    try
        cases=["stage\tid\tcore_speed_parameter_m_s\tLD_m\tdeep_jet_fraction\tdeep_cyclone_fraction\tperturbation\tm_min\tm_max\tseed\tdisplacement_parameter_m\tnx\tdt_s\tduration_s\trole"]
        audit=["stage\tid\tdx_m\tcore_regularization_m\tcore_cells\tLD_cells\tf0_times_duration\tnominal_steps\texplicit_damping_bound\tregime_status"]
        for stage in stages
            mkpath(joinpath(staging,stage)); hashes=String[]
            for m in filter(m->m.stage==stage,members)
                path=joinpath(staging,stage,m.id*".toml")
                open(path,"w") do io
                    println(io,"# Regularized QG analogue of the 2017 circular-flow initial-value experiment; not an exact shallow-water reproduction.")
                    TOML.print(io,m.cfg;sorted=true)
                end
                Validator.load_config(path)
                push!(hashes,bytes2hex(sha256(read(path)))*"\t"*basename(path))
                c=m.cfg; a=c["initial"]; d=c["physics"]; g=c["grid"]; n=c["numerics"]
                push!(cases,join((stage,m.id,a["cyclone_speed_m_s"],d["deformation_radius_m"],d["deep_jet_fraction"],d["deep_cyclone_fraction"],a["perturbation"],a["m_min"],a["m_max"],c["experiment"]["seed"],a["displacement_m"],g["nx"],n["dt_s"],n["duration_s"],m.role),'\t'))
                dx=g["length_m"]/g["nx"]
                push!(audit,join((stage,m.id,dx,a["core_regularization_m"],a["core_regularization_m"]/dx,d["deformation_radius_m"]/dx,PAPER_F0*n["duration_s"],ceil(Int,n["duration_s"]/n["dt_s"]),damping_bound(c),"declarative audit only; inspect projected and evolved hprime/H and vorticity/f"),'\t'))
            end
            write(joinpath(staging,stage,"members.tsv"),join(hashes,"\n")*"\n")
        end
        write(joinpath(staging,"cases.tsv"),join(cases,"\n")*"\n")
        write(joinpath(staging,"physics-audit.tsv"),join(audit,"\n")*"\n")
        write(joinpath(staging,"stages.tsv"),join(["$s\t$(count(m->m.stage==s,members))\tSubmit explicitly after reviewing preceding evidence" for s in stages],"\n")*"\n")
        provenance=["Project.toml","Manifest.toml","analysis-requirements.txt","scripts/make_suite.jl","scripts/make_experiments.jl",
            "scripts/make_stabilization.jl","scripts/make_rostami.jl","scripts/run.jl","scripts/validate_stabilization.jl","scripts/check_rostami_initialization.jl",
            "scripts/radial_stability.py","scripts/stabilization_linear.py","scripts/rostami_audit.py","scripts/summarize_experiments.py","scripts/field_mode_check.py",
            "slurm/common.sh","slurm/ensemble.sbatch","slurm/unity/environment.sh","slurm/unity/submit.sh","slurm/unity/ensemble.sbatch",
            "slurm/unity/experiments_01.sh","slurm/unity/experiments_01.py","slurm/unity/stabilization.sh","slurm/unity/stabilization.py",
            "slurm/unity/stabilization-validation.sbatch","slurm/unity/rostami.sh","slurm/unity/rostami.py"]
        append!(provenance,sort(["src/"*f for f in readdir(joinpath(ROOT,"src")) if endswith(f,".jl")]))
        write(joinpath(staging,"provenance.tsv"),join([bytes2hex(sha256(read(joinpath(ROOT,p))))*"\t"*p for p in provenance],"\n")*"\n")
        meta=Dict("schema_version"=>1,"campaign"=>"rostami2017_regularized_qg_analogue",
            "source_sha256"=>Legacy.source_digest(),"manifest_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml")))),
            "pilot_sha256"=>bytes2hex(sha256(read(pilot_path))),"paper_doi"=>"10.1016/j.icarus.2017.06.006",
            "stages"=>[Dict("name"=>s,"count"=>count(m->m.stage==s,members)) for s in stages],
            "scope"=>"Reconstruct circular upper velocity profiles as balanced QG initial data; no mature hexagon is imported and no exact shallow-water reproduction is claimed.",
            "deep_flow"=>"psi_d = deep_jet_fraction * unperturbed jet streamfunction + deep_cyclone_fraction * unperturbed regularized cyclone streamfunction, with identical taper and gauge; prescribed, not evolving",
            "convection"=>"off in every case; evaluate unforced formation, persistence, drift and numerical sensitivity first",
            "frame"=>Dict("rotation_rate_s"=>PAPER_F0/2,"positive_theta"=>"eastward assumption; register to goflow coordinates and epoch",
                "reference_period_s"=>38363.01,"reference_angular_speed_in_model_s"=>2pi/38363.01-PAPER_F0/2),
            "interpretation"=>Dict("duration"=>"3.125e6 seconds = 1000/f0, not 1000 inertial periods; finite persistence is not a stability theorem",
                "regularization"=>"epsilon varies the inner-core model; compare epsilon250 and epsilon500 at nx2048 before attributing changes to core physics; do not extrapolate toward the singular fit as physical convergence",
                "normalization"=>"10 km is the displacement parameter; broadband coefficients are divided by sqrt(number of modes). At the jet, the nominal delta-psi/jet-speed RMS is 10 km / sqrt(2), not a measured contour displacement.",
                "regime"=>"physics-audit.tsv is a declarative scale/time-step audit, not a QG-validity certification; inspect actual projected/evolved thickness and vorticity/f, especially the cyclone",
                "phase"=>"Measure drift only for a connected coherent jet; retain ambient PV and reference rotation and compare upper winds as well as phase speed",
                "numerics"=>"core/epsilon500km_nx2048 is the fixed-cutoff mesh check; numerics/nx2048_fixed_nu separates resolution from viscosity; outer taper differs from the paper's damping sponge"))
        open(joinpath(staging,"campaign.toml"),"w") do io; TOML.print(io,meta;sorted=true); end
        write(joinpath(staging,"README.txt"),"$(length(members)) regularized QG configurations. No jobs submitted.\nStart with validation, review, then controls and formation. Every later stage requires explicit submission.\nThis is a QG analogue of the 2017 circular-flow experiment, not a saved mature hexagon or full shallow-water reproduction.\nAll convection is off. Paper-like duration is 1000/f0 = 3.125e6 s.\nCore regularization changes the physical model; outer boundary and dissipation differ from the paper.\nSee docs/unity-rostami.md. Never resume another campaign's checkpoints.\n")
        mv(staging,out;force=false)
    finally
        isdir(staging)&&rm(staging;recursive=true)
    end
    (;directory=out,count=length(members),stages=Dict(s=>count(m->m.stage==s,members) for s in stages))
end

function main(args=ARGS)
    if args in (["--help"],["-h"])
        println("julia scripts/make_rostami.jl PILOT.toml NEW_DIRECTORY\nPrepare 32 unforced QG-analogue configurations; no jobs submitted.")
        return
    end
    length(args)==2 || throw(ArgumentError("Use --help"))
    result=make_campaign(args[1],args[2])
    println("Prepared $(result.count) configurations in $(result.directory). No jobs submitted.")
end
end
abspath(PROGRAM_FILE)==(@__FILE__) && RostamiExperiments.main()
