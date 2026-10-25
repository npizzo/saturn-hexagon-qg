#!/usr/bin/env julia

"""Staged lower-layer experiments for the existing annular-jet initializer."""
module LowerLayerExperiments

using TOML, SHA
include(joinpath(@__DIR__, "make_suite.jl"))
const Validator = SuiteGenerator.ConfigurationValidation
const ROOT = dirname(@__DIR__)

function physics_cases()
    cases = [(id="bt", ld=Inf, ratio=0.0)]
    for ld in (0.5e6, 1e6, 2e6, 4e6), ratio in (-0.5, 0.0, 0.5, 1.0)
        rtag = ratio < 0 ? "dn050" : ratio == 0 ? "d000" : ratio == 0.5 ? "dp050" : "dp100"
        push!(cases, (id="ld"*lpad(string(Int(ld/1000)),4,'0')*"_"*rtag, ld=ld, ratio=ratio))
    end
    return cases
end

const FIRST_CASES = ["bt", "ld2000_d000", "ld2000_dp050", "ld4000_d000", "ld4000_dp050", "ld4000_dp100"]
const NUMERIC_ANCHORS = ["ld2000_dp050", "ld4000_dp050"]

function source_digest()
    paths = sort(filter(p -> endswith(p,".jl"), readdir(joinpath(ROOT,"src");join=true)))
    bytes2hex(sha256(join([basename(p)*"\n"*read(p,String) for p in paths],"\n")))
end

function base_config(pilot)
    c = deepcopy(pilot)
    # Explicit protocol settings: physical upper winds and reference frame stay
    # those of the validated pilot. The legacy sphere patch is not substituted.
    c["grid"]["nx"] = 1024
    c["initial"]["perturbation"] = "broadband"
    c["initial"]["displacement_m"] = 10000.0
    c["initial"]["m_max"] = 32
    c["initial"]["mode"] = 6
    c["experiment"]["purpose"] = "selection"
    c["experiment"]["seed"] = 42
    c["boundary"]["sponge_rate_s"] = 0.0
    c["numerics"]["dt_s"] = 50.0
    c["numerics"]["duration_s"] = 500000.0
    c["numerics"]["cutoff_damping_time_s"] = 100000.0
    c["numerics"]["hyper_order"] = 4
    c["convection"]["enabled"] = false
    c["output"]["diagnostic_interval_s"] = 2000.0
    c["output"]["snapshot_interval_s"] = 50000.0
    c["output"]["checkpoint_interval_s"] = 100000.0
    c["output"]["m_max"] = 32
    c["output"]["ntheta"] = 720
    c["output"]["nr"] = 256
    return c
end

"Axisymmetric base-state QG thickness indicator; excludes perturbations and evolution."
function physical_regime(c; samples=8193)
    a, d, b = c["initial"], c["physics"], c["boundary"]
    f0 = 2c["planet"]["rotation_rate_s"]
    rossby = a["jet_speed_m_s"]/(f0*a["jet_width_m"])
    ld = d["deformation_radius_m"]
    isinf(ld) && return (; rossby, min_h_over_H=NaN, max_abs_hprime_over_H=NaN,
                         flag="barotropic_limit_no_finite_layer_inference")
    samples >= 3 || throw(ArgumentError("At least three radial samples required"))
    r = range(0.0,b["outer_radius_m"];length=samples)
    ring(r,U,R,w) = U*w*(tanh((r^2-R^2)/(2R*w))-1)
    window(r) = r <= b["inner_radius_m"] ? 1.0 : r >= b["outer_radius_m"] ? 0.0 :
        (t=(r-b["inner_radius_m"])/(b["outer_radius_m"]-b["inner_radius_m"]); 1-10t^3+15t^4-6t^5)
    delta = [window(x)*(ring(x,a["jet_speed_m_s"],a["jet_radius_m"],a["jet_width_m"]) -
        a["cyclone_speed_m_s"]*a["cyclone_radius_m"]*exp(0.5-x^2/(2a["cyclone_radius_m"]^2)) -
        ring(x,d["deep_jet_speed_m_s"],d["deep_jet_radius_m"],d["deep_jet_width_m"])) for x in r]
    # The taper disk lies within the square. Outside it both unshifted fields are
    # zero; integrating over the disk gives the full-square mean exactly in the
    # continuum limit. The final radial endpoint includes the outside value.
    integrand = r .* delta
    area_mean = 2pi*step(r)*(sum(integrand)-(first(integrand)+last(integrand))/2)/c["grid"]["length_m"]^2
    frac = (delta .- area_mean)./(f0*ld^2)
    min_h_over_H = 1+minimum(frac)
    return (; rossby, min_h_over_H, max_abs_hprime_over_H=maximum(abs,frac),
             flag=min_h_over_H <= 0 ? "nonpositive_inferred_base_depth" : "small_displacement_assumption_requires_review")
end

"""
    make_experiments(pilot_path, new_directory; candidates=[], seeds=42:44,
                     long_duration_s=2e6)

Generate 6 initial growth screens, 11 remaining parameter screens and 8
numerical sensitivity cases. Optional explicit candidate IDs generate seeded
persistence, finite-amplitude competition and matched convection stages.
Portable members.tsv files contain SHA256<TAB>basename; the Unity wrapper
freezes absolute scheduler manifests after transfer. No runs are launched.
"""
function make_experiments(pilot_path, new_directory; candidates=String[], seeds=42:44,
                          long_duration_s=2e6)
    out = abspath(new_directory)
    (ispath(out) || islink(out)) && throw(ArgumentError("Output already exists; use a NEW campaign directory"))
    any(c -> c in ('\t','\n','\r'),out) && throw(ArgumentError("Output path contains a tab or line break"))
    pilot = Validator.load_config(pilot_path)
    pilot["initial"]["jet_speed_m_s"] > 0 || throw(ArgumentError("A positive upper jet is required"))
    pilot["convection"]["bias_max"] > 0 || throw(ArgumentError("A positive bias is required for distinct forcing controls"))
    allunique(candidates) || throw(ArgumentError("Candidate IDs must be unique"))
    seedlist = collect(seeds)
    !isempty(seedlist) && all(s -> s isa Integer && !(s isa Bool) && s >= 0, seedlist) && allunique(seedlist) || throw(ArgumentError("Seeds must be unique nonnegative integers"))
    isfinite(long_duration_s) && long_duration_s >= 2e6 || throw(ArgumentError("Long stages require at least 2e6 s"))
    physics = physics_cases()
    lookup = Dict(p.id=>p for p in physics)
    all(id -> haskey(lookup,id), candidates) || throw(ArgumentError("Unknown candidate ID; use --list-cases"))
    base = base_config(pilot)
    members = NamedTuple[]
    function add(stage, p, tag; family="broadband", seed=42, displacement=10000.0,
                 duration=500000.0, nx=1024, dt=50.0, tau=100000.0,
                 forcing="off", role="")
        c = deepcopy(base)
        id = p.id*"__"*tag
        c["experiment"]["name"] = stage*"__"*id
        c["experiment"]["purpose"] = family == "mode" ? "maintenance" : "selection"
        c["experiment"]["seed"] = seed
        c["physics"]["deformation_radius_m"] = p.ld
        c["physics"]["deep_jet_speed_m_s"] = p.ratio*base["initial"]["jet_speed_m_s"]
        c["initial"]["perturbation"] = family
        c["initial"]["displacement_m"] = displacement
        c["grid"]["nx"] = nx
        c["numerics"]["dt_s"] = dt
        c["numerics"]["duration_s"] = duration
        c["numerics"]["cutoff_damping_time_s"] = tau
        c["convection"]["enabled"] = forcing != "off"
        c["convection"]["bias_max"] = forcing == "unbiased" ? 0.0 : pilot["convection"]["bias_max"]
        if duration >= 2e6
            c["output"]["snapshot_interval_s"] = 200000.0
            c["output"]["checkpoint_interval_s"] = 200000.0
        end
        Validator.validate_config(c)
        push!(members,(;stage,id,physics=p.id,family,seed,forcing,role,cfg=c))
    end
    for id in FIRST_CASES
        add("screen",lookup[id],"bb10km_s42";role="Initial screen: detect rapid growth and measure provisional mode drift")
    end
    for p in physics
        p.id in FIRST_CASES && continue
        add("expanded",p,"bb10km_s42";role="Expand parameter coverage after reviewing the initial six")
    end
    for id in unique(vcat(NUMERIC_ANCHORS,collect(candidates)))
        p = lookup[id]
        add("numerics",p,"dt25";dt=25.0,role="Time-step sensitivity relative to screen")
        add("numerics",p,"nx2048_dt25";nx=2048,dt=25.0,role="Combined mesh and cutoff-dissipation refinement; not fixed-nu convergence")
        add("numerics",p,"amplitude5km";displacement=5000.0,role="Linearity and diagnostic floor relative to screen")
        add("numerics",p,"weaker_damping";tau=1e6,role="Tenfold lower physical hyperviscosity at fixed grid")
    end
    for id in candidates
        p = lookup[id]
        add("maintenance",p,"m6_100km";family="mode",displacement=1e5,duration=long_duration_s,
            role="Seeded sixfold persistence; symmetry-restricted, not a nonlinear stability proof")
        for seed in seedlist
            add("selection",p,"bb100km_s$seed";seed,displacement=1e5,duration=long_duration_s,
                role="Finite-amplitude competition across modes 1:32; phase realization varies")
        end
        for forcing in ("unbiased","biased"), seed in seedlist
            add("convection",p,"m6_$(forcing)_s$seed";family="mode",seed,displacement=1e5,
                duration=long_duration_s,forcing,
                role="Pair with this candidate's maintenance run as the unforced control; identical fresh upper field")
        end
    end
    # Validate the complete protocol before publishing any files.
    allunique([(m.stage,m.id) for m in members]) || error("Duplicate experiment")
    stages = unique([m.stage for m in members])
    mkpath(dirname(out))
    staging = mktempdir(dirname(out);prefix=".saturn-experiments-")
    try
        index = ["stage\tid\tphysics\tfamily\tseed\tforcing\tnx\tdt_s\tduration_s\tdisplacement_m\trole"]
        stage_rows = String[]
        for stage in stages
            mkpath(joinpath(staging,stage))
            rows = String[]
            group = filter(m -> m.stage == stage,members)
            for m in group
                filename = m.id*".toml"
                path = joinpath(staging,stage,filename)
                open(path,"w") do io
                    println(io,"# Staged annular-jet experiment; not an exact-2024 initializer or observed fit.")
                    TOML.print(io,m.cfg;sorted=true)
                end
                Validator.load_config(path)
                push!(rows,bytes2hex(sha256(read(path)))*"\t"*filename)
                push!(index,join((m.stage,m.id,m.physics,m.family,m.seed,m.forcing,
                    m.cfg["grid"]["nx"],m.cfg["numerics"]["dt_s"],m.cfg["numerics"]["duration_s"],
                    m.cfg["initial"]["displacement_m"],m.role),'\t'))
            end
            write(joinpath(staging,stage,"members.tsv"),join(rows,"\n")*"\n")
            push!(stage_rows,stage*"\t"*string(length(group))*"\t"*(stage=="screen" ? "Start here" : "Review preceding evidence before submitting this stage"))
        end
        write(joinpath(staging,"stages.tsv"),join(stage_rows,"\n")*"\n")
        write(joinpath(staging,"cases.tsv"),join(index,"\n")*"\n")
        audit_rows = ["physics\tLD_m\tdeep_jet_m_s\tnominal_width_Ro\tbase_min_h_over_H\tbase_max_abs_hprime_over_H\tinterpretation"]
        for p in physics
            c=deepcopy(base)
            c["physics"]["deformation_radius_m"]=p.ld
            c["physics"]["deep_jet_speed_m_s"]=p.ratio*base["initial"]["jet_speed_m_s"]
            audit=physical_regime(c)
            value(x)=isnan(x) ? "NA" : string(x)
            push!(audit_rows,join((p.id,p.ld,c["physics"]["deep_jet_speed_m_s"],audit.rossby,
                value(audit.min_h_over_H),value(audit.max_abs_hprime_over_H),audit.flag),'\t'))
        end
        write(joinpath(staging,"physics-audit.tsv"),join(audit_rows,"\n")*"\n")
        provenance = ["Project.toml","Manifest.toml","scripts/make_suite.jl","scripts/make_experiments.jl",
            "scripts/run.jl","slurm/common.sh","slurm/ensemble.sbatch",
            "slurm/unity/environment.sh","slurm/unity/submit.sh","slurm/unity/ensemble.sbatch",
            "slurm/unity/experiments_01.sh","slurm/unity/experiments_01.py"]
        append!(provenance,sort(["src/"*f for f in readdir(joinpath(ROOT,"src")) if endswith(f,".jl")]))
        write(joinpath(staging,"provenance.tsv"),join([bytes2hex(sha256(read(joinpath(ROOT,p))))*"\t"*p for p in provenance],"\n")*"\n")
        r = base["initial"]["jet_radius_m"]
        omega = base["planet"]["rotation_rate_s"]
        meta = Dict("schema_version"=>1,"initializer"=>"existing annular jet plus central cyclone",
            "exact_2024_family"=>false,"source_sha256"=>source_digest(),
            "manifest_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml")))),
            "pilot_sha256"=>bytes2hex(sha256(read(pilot_path))),"candidates"=>collect(candidates),
            "stages"=>[Dict("name"=>s,"count"=>count(m->m.stage==s,members)) for s in stages],
            "frame"=>Dict("assumption"=>"positive theta eastward; configured rotation defines coordinate frame",
                "model_rotation_rate_s"=>omega,"system_iii_period_s"=>38362.4,
                "observed_2008_2014_period_s"=>38363.01,
                "observed_angular_speed_model_s"=>2pi/38363.01-omega,
                "observed_speed_at_model_ring_m_s"=>r*(2pi/38363.01-omega),
                "status"=>"conditional comparison; register to the goflow epoch and coordinate convention"),
            "interpretation"=>Dict("nominal_jet_circuit_s"=>2pi*r/base["initial"]["jet_speed_m_s"],
                "screen_claim"=>"rapid-growth screen; a quiet 500000-s trajectory is not stability",
                "long_claim"=>"finite-amplitude persistence and mode competition, followed by matched convection",
                "numerics"=>"2048 keeps cutoff damping time fixed; physical nu decreases approximately 256-fold",
                "diagnostics"=>"runtime jet phase is provisional; validate a connected jet and field coefficients",
                "physical_regime"=>"physics-audit.tsv: axisymmetric base hprime/H=(psi-psid)/(2 Omega LD^2), zero-domain-mean convention; excludes perturbations/evolution; order-one values violate small-displacement QG",
                "forcing"=>"matched seeds do not guarantee identical accepted injections or energy input"))
        open(joinpath(staging,"campaign.toml"),"w") do io; TOML.print(io,meta;sorted=true); end
        write(joinpath(staging,"README.txt"),
            "Prepared $(length(members)) configurations. No jobs launched.\n"*
            "Submit only the six-case screen first; compare numerics before interpreting weak growth.\n"*
            "expanded is the remainder of the 17-case physics matrix, not a required immediate launch.\n"*
            "Candidate stages exist only when explicitly requested by --candidates.\n"*
            "Read physics-audit.tsv: small-LD cases can violate small-interface-displacement QG and are mathematical stress tests.\n"*
            "The maintenance run is also the unforced control for that candidate's six forced convection runs.\n"*
            "Use docs/experiments-01.md for the protocol, limitations and Unity commands.\n"*
            "members.tsv paths are relative to each stage; provenance.tsv paths are relative to the project.\n")
        mv(staging,out;force=false)
    finally
        isdir(staging) && rm(staging;recursive=true)
    end
    return (;directory=out,count=length(members),stages,counts=Dict(s=>count(m->m.stage==s,members) for s in stages))
end

function main(args=ARGS)
    if args == ["--list-cases"]
        foreach(p -> println(p.id,"\tLD=",p.ld," m\tdeep/upper=",p.ratio),physics_cases())
        return
    end
    if args in (["--help"],["-h"])
        println("Usage: julia scripts/make_experiments.jl PILOT.toml NEW_DIRECTORY [--candidates ID,ID] [--seeds 42:44] [--long-duration-s 2000000]")
        println("Default: screen=6, expanded=11, numerics=8. No simulation or submission.")
        println("--candidates adds maintenance, selection and convection stages only for the named physics IDs.")
        println("Use --list-cases for IDs. Portable hashes allow transfer before scheduler preparation.")
        return
    end
    length(args)>=2 || throw(ArgumentError("Expected PILOT.toml NEW_DIRECTORY; use --help"))
    candidates=String[]; seeds=collect(42:44); duration=2e6; seen=Set{String}(); i=3
    while i <= length(args)
        option=args[i]
        option in seen && throw(ArgumentError("Duplicate option $option")); push!(seen,option)
        i < length(args) || throw(ArgumentError("Missing value for $option"))
        if option == "--candidates"
            candidates=String.(split(args[i+1],','))
        elseif option == "--seeds"
            seeds=SuiteGenerator._seeds(args[i+1])
        elseif option == "--long-duration-s"
            duration=parse(Float64,args[i+1])
        else
            throw(ArgumentError("Unknown option $option"))
        end
        i+=2
    end
    result=make_experiments(args[1],args[2];candidates,seeds,long_duration_s=duration)
    println("Prepared ",result.count," configurations in ",result.directory)
    foreach(s -> println(s,": ",result.counts[s]),result.stages)
    println("No simulations or jobs launched. Start with screen; see docs/experiments-01.md.")
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        LowerLayerExperiments.main()
    catch err
        println(stderr,"Experiment preparation failed: ",sprint(showerror,err))
        exit(1)
    end
end
