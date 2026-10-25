#!/usr/bin/env julia
"""Unforced circular-jet stabilization campaign. Never submits jobs."""
module StabilizationExperiments
using TOML, SHA
include(joinpath(@__DIR__, "make_experiments.jl"))
const Legacy = LowerLayerExperiments
const Validator = Legacy.Validator
const ROOT = dirname(@__DIR__)

caseid(ld, alpha, width) = (isinf(ld) ? "bt" : "ld"*string(round(Int,ld/1000))) *
    "_a"*lpad(string(round(Int,100alpha)),3,'0')*"_w"*string(round(Int,width/1000))
case(ld, alpha; width=1e6) = (;id=caseid(ld,alpha,width),ld,alpha,width)
function physics_cases()
    cases = [case(Inf,0.0), case(Inf,1.0)]
    append!(cases,[case(ld,a) for ld in (2e6,3.2e6,4e6) for a in (0.0,0.5,1.0)])
    append!(cases,[case(ld,1.0) for ld in (0.5e6,0.65e6,0.8e6,1e6)])
    append!(cases,[case(3.2e6,a;width=w) for w in (1.5e6,2e6,3e6) for a in (0.0,1.0)])
    cases
end

"Same worst-mode explicit damping bound used by the time integrator (no sponge)."
function damping_bound(c)
    n=c["grid"]["nx"];p=c["numerics"]["hyper_order"]
    c["numerics"]["dt_s"]/c["numerics"]["cutoff_damping_time_s"] *
        (2*(n/2)^2/floor(n/3)^2)^p
end

"Continuum base thickness plus a conservative perturbation bound, not a stability test."
function physical_audit(c)
    a,d,b=c["initial"],c["physics"],c["boundary"]
    f0=2c["planet"]["rotation_rate_s"]
    ld=d["deformation_radius_m"]
    ro=abs(a["jet_speed_m_s"])/(f0*a["jet_width_m"])
    dx=c["grid"]["length_m"]/c["grid"]["nx"]
    isinf(ld) && return (;ro,ld_cells=Inf,base_max=NaN,base_min=NaN,perturb_bound=NaN,
        initial_max_bound=NaN,initial_min_bound=NaN,flag="barotropic_no_layer_depth_inference")
    r=range(0.0,b["outer_radius_m"];length=16385)
    window(x)=x<=b["inner_radius_m"] ? 1.0 : x>=b["outer_radius_m"] ? 0.0 :
        (t=(x-b["inner_radius_m"])/(b["outer_radius_m"]-b["inner_radius_m"]);1-10t^3+15t^4-6t^5)
    upper(x)=a["jet_speed_m_s"]*a["jet_width_m"]*(tanh((x^2-a["jet_radius_m"]^2)/(2a["jet_radius_m"]*a["jet_width_m"]))-1) -
        a["cyclone_speed_m_s"]*a["cyclone_radius_m"]*exp(0.5-x^2/(2a["cyclone_radius_m"]^2))
    delta=[(1-d["deep_match_fraction"])*window(x)*upper(x) for x in r]
    integrand=r.*delta
    gauge=2pi*step(r)*(sum(integrand)-(first(integrand)+last(integrand))/2)/c["grid"]["length_m"]^2
    chi=(delta.-gauge)./(f0*ld^2)
    # Every analytic mode envelope and sech^2 factor is at most one.
    # Factor two includes a conservative bound on the subsequent mean removal.
    nm=a["perturbation"]=="broadband" ? a["m_max"] : a["perturbation"]=="mode" ? 1 : 0
    bound=2abs(a["displacement_m"]*a["jet_speed_m_s"])*sqrt(nm)/(f0*ld^2)
    base_max=maximum(abs,chi);base_min=1+minimum(chi)
    initial_max_bound=base_max+bound;initial_min_bound=base_min-bound
    flag=initial_min_bound<=0 ? "bound_cannot_guarantee_positive_depth" :
        initial_max_bound>=0.2 ? "review_small_displacement_assumption" : "small_initial_thickness_bound_only"
    (;ro,ld_cells=ld/dx,base_max,base_min,perturb_bound=bound,initial_max_bound,initial_min_bound,flag)
end

function make_campaign(pilot_path,new_directory; candidates=String[],seeds=42:44,long_duration_s=5e6)
    out=abspath(new_directory)
    (ispath(out)||islink(out)) && throw(ArgumentError("Choose a NEW campaign directory; existing campaigns are frozen"))
    any(ch->ch in ('\t','\r','\n'),out) && throw(ArgumentError("Invalid output path"))
    allunique(candidates) || throw(ArgumentError("Candidate IDs must be unique"))
    seedlist=collect(seeds)
    !isempty(seedlist) && all(s->s isa Integer && !(s isa Bool) && s>=0,seedlist) && allunique(seedlist) || throw(ArgumentError("Invalid seeds"))
    isfinite(long_duration_s) && long_duration_s>=5e6 || throw(ArgumentError("Candidate duration must be at least 5e6 s; this alone does not establish stability"))
    pilot=Validator.load_config(pilot_path)
    pilot["initial"]["jet_speed_m_s"]>0 || throw(ArgumentError("Positive annular jet required"))
    base=Legacy.base_config(pilot)
    base["physics"]["deep_profile"]="matched_upper"
    base["physics"]["deep_jet_speed_m_s"]=0.0
    base["physics"]["deep_match_fraction"]=1.0
    lookup=Dict(p.id=>p for p in physics_cases())
    all(id->haskey(lookup,id),candidates) || throw(ArgumentError("Unknown candidate; use --list-cases"))
    members=NamedTuple[]
    function add(stage,p,tag;family="broadband",displacement=1e4,nx=1024,dt=50.0,
                 duration=5e5,tau=1e5,seed=42,box=base["grid"]["length_m"],
                 inner=base["boundary"]["inner_radius_m"],outer=base["boundary"]["outer_radius_m"],role="")
        c=deepcopy(base)
        id=p.id*"__"*tag
        c["experiment"]["name"]=stage*"__"*id
        c["experiment"]["purpose"]=stage=="validation" ? "smoke" : family=="mode" ? "maintenance" : "selection"
        c["experiment"]["seed"]=seed
        c["physics"]["deformation_radius_m"]=p.ld
        c["physics"]["deep_match_fraction"]=p.alpha
        c["initial"]["jet_width_m"]=p.width
        c["initial"]["perturbation"]=family
        c["initial"]["displacement_m"]=displacement
        c["grid"]["nx"]=nx
        c["grid"]["length_m"]=box
        c["boundary"]["inner_radius_m"]=inner
        c["boundary"]["outer_radius_m"]=outer
        c["numerics"]["dt_s"]=dt
        c["numerics"]["duration_s"]=duration
        c["numerics"]["cutoff_damping_time_s"]=tau
        c["output"]["jet_search_halfwidth_m"]=min(8e6,max(4e6,2.5p.width))
        c["output"]["diagnostic_interval_s"]=min(2000.0,duration/4)
        c["output"]["snapshot_interval_s"]=min(duration/4,duration>=5e6 ? 2e5 : 5e4)
        c["output"]["checkpoint_interval_s"]=min(duration/2,2e5)
        Validator.validate_config(c)
        damping_bound(c)<=2 || throw(ArgumentError("Explicit damping bound exceeds 2 for $id; reduce timestep"))
        push!(members,(;stage,id,physics=p.id,role,cfg=c))
    end
    # These are software/GPU checks; they do not represent resolved science runs.
    for p in (case(Inf,0.0),case(3.2e6,0.5),case(0.5e6,1.0))
        add("validation",p,"gpu_smoke";nx=256,dt=10.0,duration=400.0,role="Short CPU/GPU agreement and output check")
    end
    for p in (case(Inf,0.0),case(Inf,1.0),case(3.2e6,0.0),case(3.2e6,1.0))
        add("controls",p,"axisymmetric";family="none",displacement=0.0,duration=1e5,
            role="Check axisymmetric numerical drift; BT alpha pair must agree")
    end
    add("controls",case(Inf,0.0),"bb10km_s42";role="Perturbed barotropic instability reference")
    for ld in (2e6,3.2e6,4e6),alpha in (0.0,0.5,1.0)
        add("screen",case(ld,alpha),"bb10km_s42";role="Fixed upper winds; change finite LD and whole deep profile amplitude")
    end
    for ld in (0.5e6,0.65e6,0.8e6,1e6)
        add("transition",case(ld,1.0),"bb10km_s42";dt=25.0,duration=1e6,
            role="Fully matched jet and core; small-LD transition, quiet runs require longer and finer checks")
    end
    for w in (1.5e6,2e6,3e6),alpha in (0.0,1.0)
        add("width",case(3.2e6,alpha;width=w),"bb10km_s42";role="Width sensitivity changes the upper jet; not an observational fit")
    end
    for id in unique(vcat([case(3.2e6,1.0).id,case(0.5e6,1.0).id],collect(candidates)))
        p=lookup[id]
        duration=p.ld<2e6 ? 1e6 : 5e5
        dt=p.ld<2e6 ? 25.0 : 50.0
        add("numerics",p,"half_dt";dt=dt/2,duration,role="Timestep refinement at unchanged grid and physical viscosity")
        add("numerics",p,"nx2048_fixed_cutoff";nx=2048,dt=dt/2,duration,
            role="Double mesh; same cutoff damping time; physical hyperviscosity reduced 256-fold")
        add("numerics",p,"nx2048_fixed_nu";nx=2048,dt=dt/32,duration,tau=1e5/256,
            role="Double mesh at unchanged physical hyperviscosity; compare with fixed-cutoff refinement")
        add("numerics",p,"amplitude5km";dt,duration,displacement=5e3,role="Half-amplitude linearity and diagnostic-floor check")
        add("numerics",p,"weaker_damping";dt,duration,tau=1e6,role="Tenfold weaker hyperviscosity at fixed grid")
    end
    for id in candidates
        p=lookup[id];dt=25.0
        add("maintenance",p,"m6_100km";family="mode",displacement=1e5,duration=long_duration_s,dt,
            role="Seeded sixfold persistence; exact sixfold symmetry is not nonlinear stability")
        for seed in seedlist
            add("competition",p,"bb100km_s$seed";seed,displacement=1e5,duration=long_duration_s,dt,
                role="Finite-amplitude modes 1:32; independent phase seeds; not a published nonlinear state")
        end
    end
    for p in (case(0.5e6,1.0),case(3.2e6,1.0))
        k_old=2pi/base["grid"]["length_m"]*floor(1024/3)
        k_new=2pi/75e6*floor(1280/3)
        add("boundary",p,"box75Mm";nx=1280,box=75e6,inner=30e6,outer=35e6,dt=25.0,
            duration=p.ld<2e6 ? 1e6 : 5e5,tau=1e5*(k_old/k_new)^8,
            role="Larger box and outward continuation at fixed dx and physical hyperviscosity; inspect interface gauge too")
    end
    stages=unique(m.stage for m in members)
    allunique((m.stage,m.id) for m in members) || error("Duplicate case")
    mkpath(dirname(out));staging=mktempdir(dirname(out);prefix=".stabilization-")
    try
        rows=["stage\tid\tphysics\tperturbation\tnx\tdt_s\tduration_s\trole"]
        audits=["stage\tid\tnominal_width_Ro\tLD_cells\tbase_max_abs_hprime_over_H\tbase_min_h_over_H\tperturbation_bound\tinitial_max_abs_bound\tinitial_min_h_bound\tflag"]
        for stage in stages
            mkpath(joinpath(staging,stage));hashes=String[]
            for m in filter(m->m.stage==stage,members)
                path=joinpath(staging,stage,m.id*".toml")
                open(path,"w") do io
                    println(io,"# Unforced circular jet plus cyclone; not a published nonlinear hexagon or observed fit.")
                    TOML.print(io,m.cfg;sorted=true)
                end
                Validator.load_config(path)
                push!(hashes,bytes2hex(sha256(read(path)))*"\t"*basename(path))
                c=m.cfg
                push!(rows,join((m.stage,m.id,m.physics,c["initial"]["perturbation"],c["grid"]["nx"],c["numerics"]["dt_s"],c["numerics"]["duration_s"],m.role),'\t'))
                a=physical_audit(c)
                push!(audits,join((m.stage,m.id,a.ro,a.ld_cells,a.base_max,a.base_min,a.perturb_bound,a.initial_max_bound,a.initial_min_bound,a.flag),'\t'))
            end
            write(joinpath(staging,stage,"members.tsv"),join(hashes,"\n")*"\n")
        end
        write(joinpath(staging,"stages.tsv"),join(["$s\t$(count(m->m.stage==s,members))\tSubmit explicitly after reviewing preceding stage" for s in stages],"\n")*"\n")
        write(joinpath(staging,"cases.tsv"),join(rows,"\n")*"\n")
        write(joinpath(staging,"physics-audit.tsv"),join(audits,"\n")*"\n")
        provenance=["Project.toml","Manifest.toml","analysis-requirements.txt","scripts/make_suite.jl","scripts/make_experiments.jl",
            "scripts/make_stabilization.jl","scripts/run.jl","scripts/validate_gpu.jl","scripts/validate_stabilization.jl",
            "scripts/radial_stability.py","scripts/stabilization_linear.py","scripts/summarize_experiments.py","scripts/field_mode_check.py",
            "slurm/common.sh","slurm/ensemble.sbatch","slurm/unity/environment.sh","slurm/unity/submit.sh",
            "slurm/unity/ensemble.sbatch","slurm/unity/experiments_01.sh","slurm/unity/experiments_01.py",
            "slurm/unity/stabilization.sh","slurm/unity/stabilization.py","slurm/unity/stabilization-validation.sbatch"]
        append!(provenance,sort(["src/"*f for f in readdir(joinpath(ROOT,"src")) if endswith(f,".jl")]))
        write(joinpath(staging,"provenance.tsv"),join([bytes2hex(sha256(read(joinpath(ROOT,p))))*"\t"*p for p in provenance],"\n")*"\n")
        omega=base["planet"]["rotation_rate_s"]
        meta=Dict("schema_version"=>1,"campaign"=>"unforced_whole_profile_stabilization",
            "source_sha256"=>Legacy.source_digest(),"manifest_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml")))),
            "pilot_sha256"=>bytes2hex(sha256(read(pilot_path))),"candidates"=>collect(candidates),
            "stages"=>[Dict("name"=>s,"count"=>count(m->m.stage==s,members)) for s in stages],
            "deep_flow"=>"psi_d=alpha times the unperturbed upper annulus plus cyclone, with the same taper and zero-mean gauge; no deep hexagon",
            "convection"=>"off in every case; no convection-maintenance claim",
            "published_states"=>Dict("pizzo_salmon_2024"=>"reconstruction_required: ambient PV, rotation rate, global matching and physical scaling",
                "rostami_zeitlin_spiga_2017"=>"reconstruction_required: balanced QG counterpart of the RSW state"),
            "frame"=>Dict("rotation_rate_s"=>omega,"positive_theta"=>"eastward assumption; register to goflow convention",
                "reference_period_s"=>38363.01,"reference_angular_speed_in_model_s"=>2pi/38363.01-omega),
            "interpretation"=>Dict("duration"=>"growth screens and candidate persistence; quiet is not stable",
                "regime"=>"continuum initial thickness bounds in physics-audit.tsv; actual projected/evolved thickness and vorticity/f in diagnostics.csv; no certification of QG validity",
                "phase"=>"measure only on a connected coherent jet and cross-check field modes; retain ambient PV and frame",
                "small_LD"=>"weak modes and critical layers require radial and 2D mesh/domain refinement",
                "numerics"=>"both fixed physical hyperviscosity and fixed cutoff damping refinements included"))
        open(joinpath(staging,"campaign.toml"),"w") do io;TOML.print(io,meta;sorted=true);end
        write(joinpath(staging,"README.txt"),"$(length(members)) configurations. No jobs submitted.\nStart with GPU preflight and validation, then controls and screen. Submit other stages only after review.\nAll convection is off. Lower flow includes the central cyclone.\nPublished nonlinear states remain gated; ambient PV and rotation are essential to their reconstruction.\nSee docs/unity-stabilization.md and docs/published-state-readiness.md.\n")
        mv(staging,out;force=false)
    finally
        isdir(staging)&&rm(staging;recursive=true)
    end
    (;directory=out,count=length(members),stages=Dict(s=>count(m->m.stage==s,members) for s in stages))
end

function main(args=ARGS)
    if args==["--list-cases"]
        foreach(p->println(p.id),physics_cases());return
    elseif args==["--help"]
        println("julia scripts/make_stabilization.jl PILOT.toml NEW_DIRECTORY [--candidates ID,ID] [--seeds 42:44] [--long-duration 5000000]")
        println("Default 39 unforced configs; explicit candidates add numerical checks, m6 maintenance and broadband competition.")
        return
    end
    length(args)>=2 || throw(ArgumentError("Use --help"))
    candidates=String[];seeds=42:44;duration=5e6;i=3
    while i<=length(args)
        i<length(args)||throw(ArgumentError("Missing option value"))
        key,value=args[i:i+1]
        if key=="--candidates";candidates=String.(split(value,','))
        elseif key=="--seeds";seeds=Legacy.SuiteGenerator._seeds(value)
        elseif key=="--long-duration";duration=parse(Float64,value)
        else;throw(ArgumentError("Unknown option $key"));end
        i+=2
    end
    result=make_campaign(args[1],args[2];candidates,seeds,long_duration_s=duration)
    println("Prepared $(result.count) configurations in $(result.directory). No jobs submitted.")
end
end
if abspath(PROGRAM_FILE) == @__FILE__
    StabilizationExperiments.main()
end
