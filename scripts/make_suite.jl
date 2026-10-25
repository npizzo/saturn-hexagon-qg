#!/usr/bin/env julia

module SuiteGenerator

using TOML, SHA

# Reuse the package's actual configuration validator without importing FFT/GPU
# dependencies or constructing simulation grids. All these imports are stdlib.
module ConfigurationValidation
using TOML
include(joinpath(@__DIR__, "..", "src", "Convection.jl"))
using .Convection: ConvectionState
include(joinpath(@__DIR__, "..", "src", "Configuration.jl"))
end

export make_suite, main

function _seeds(text)
    values = if occursin(':',text)
        parts=split(text,':')
        length(parts)==2 || throw(ArgumentError("Seed range must be FIRST:LAST"))
        lo,hi=tryparse.(Int,parts)
        lo !== nothing && hi !== nothing && lo <= hi || throw(ArgumentError("Invalid seed range"))
        collect(lo:hi)
    else
        parts=split(text,',')
        parsed=tryparse.(Int,parts)
        all(v -> v !== nothing,parsed) || throw(ArgumentError("Seeds must be comma-separated integers or FIRST:LAST"))
        Int[v for v in parsed]
    end
    !isempty(values) && all(>=(0),values) && allunique(values) || throw(ArgumentError("Seeds must be nonnegative, unique and nonempty"))
    return values
end

"""
    make_suite(pilot_path, new_directory; seeds=42:44, include_emergence=false)

Write validated TOMLs for three physics controls, three forcing controls and
maintenance/selection families. Optional emergence starts the upper layer at
rest and uses only the two forced controls. Nothing is launched. Existing output
paths are refused. `manifest.tsv` uses the scheduler's SHA256-tab-absolute-path
contract; `config_paths.txt` is a plain ordered list for inspection, not the
scheduler manifest. Generate on the target host or regenerate hashes/paths there.
"""
function make_suite(pilot_path::AbstractString, new_directory::AbstractString;
                    seeds=42:44, include_emergence::Bool=false)
    outdir=abspath(new_directory)
    any(c -> c in ('\t','\n','\r'),outdir) && throw(ArgumentError("Output path cannot contain tabs or line breaks"))
    (ispath(outdir) || islink(outdir)) && throw(ArgumentError("Output directory already exists; choose a NEW path"))
    seedlist=collect(seeds)
    !isempty(seedlist) && all(s -> s isa Integer && !(s isa Bool) && s >= 0,seedlist) && allunique(seedlist) || throw(ArgumentError("Seeds must be unique nonnegative integers"))
    pilot=TOML.parsefile(pilot_path)
    ConfigurationValidation.validate_config(pilot)
    isfinite(pilot["physics"]["deformation_radius_m"]) || throw(ArgumentError("Pilot must supply a finite deformation radius for the finite-LD controls"))
    pilot["physics"]["deep_jet_speed_m_s"] != 0 || throw(ArgumentError("Pilot must supply a nonzero deep jet for the deep-flow control"))
    pilot["convection"]["bias_max"] > 0 || throw(ArgumentError("Pilot bias_max must be positive to distinguish biased and unbiased controls"))
    pilot["initial"]["jet_speed_m_s"] > 0 || throw(ArgumentError("Pilot must supply a positive upper jet for maintenance and selection"))
    pilot["initial"]["displacement_m"] > 0 || throw(ArgumentError("Pilot must supply a nonzero maintenance/selection disturbance"))
    pilot["initial"]["m_max"] >= 6 || throw(ArgumentError("Broadband selection must include mode six; initial.m_max must be at least 6"))
    families=include_emergence ? ("maintenance","selection","emergence") : ("maintenance","selection")
    physics_controls=("barotropic_no_deep","finite_no_deep","finite_deep")
    forcing_controls=("off","unbiased","biased")
    members=NamedTuple[]
    for family in families, physics in physics_controls, forcing in forcing_controls, seed in seedlist
        family == "emergence" && forcing == "off" && continue
        cfg=deepcopy(pilot)
        name="$(family)__$(physics)__$(forcing)__seed$(seed)"
        cfg["experiment"]["name"]=name
        cfg["experiment"]["purpose"]=family
        cfg["experiment"]["seed"]=seed
        cfg["physics"]["deformation_radius_m"]=physics == "barotropic_no_deep" ? Inf : pilot["physics"]["deformation_radius_m"]
        cfg["physics"]["deep_jet_speed_m_s"]=physics == "finite_deep" ? pilot["physics"]["deep_jet_speed_m_s"] : 0.0
        cfg["convection"]["enabled"]=forcing != "off"
        cfg["convection"]["bias_max"]=forcing == "biased" ? pilot["convection"]["bias_max"] : 0.0
        if family == "maintenance"
            cfg["initial"]["perturbation"]="mode"
            cfg["initial"]["mode"]=6
        elseif family == "selection"
            cfg["initial"]["perturbation"]="broadband"
        else
            cfg["initial"]["jet_speed_m_s"]=0.0
            cfg["initial"]["cyclone_speed_m_s"]=0.0
            cfg["initial"]["perturbation"]="none"
            cfg["initial"]["displacement_m"]=0.0
        end
        ConfigurationValidation.validate_config(cfg)
        push!(members,(; name,family,physics,forcing,seed,cfg))
    end
    # Complete validation before writing. Publish the prepared directory only
    # after every round-trip check and manifest has succeeded.
    parent=dirname(outdir)
    mkpath(parent)
    staging=mktempdir(parent;prefix=".saturn-suite-")
    paths=String[]
    try
        manifest_rows=String[]
        for member in members
            filename=member.name*".toml"
            staged_path=joinpath(staging,filename)
            open(staged_path,"w") do io
                println(io,"# Generated paired control; illustrative parameters, not an observed fit.")
                TOML.print(io,member.cfg;sorted=true)
            end
            ConfigurationValidation.load_config(staged_path)
            final_path=joinpath(outdir,filename)
            push!(paths,final_path)
            push!(manifest_rows,bytes2hex(sha256(read(staged_path)))*"\t"*final_path)
        end
        write(joinpath(staging,"config_paths.txt"),join(paths,"\n")*"\n")
        write(joinpath(staging,"manifest.tsv"),join(manifest_rows,"\n")*"\n")
        open(joinpath(staging,"README.txt"),"w") do io
            println(io,"Generated $(length(members)) configurations; no simulation or scheduler job has been launched.")
            println(io,"Seeds: ",join(seedlist,", "))
            println(io,"Families: ",join(families,", "))
            println(io,"Physics: barotropic/no deep; finite deformation/no deep; finite deformation/with deep jet.")
            println(io,"Forcing: off; unbiased; biased. Emergence excludes forcing-off states and starts the upper layer at rest.")
            println(io,"Maintenance deliberately seeds mode six. Selection uses the model's broadband initializer.")
            println(io,"Emergence with a deep jet is conditional on that imposed structure; it is not untemplated jet formation.")
            println(io,"All configs passed the package's shared configuration validator and TOML round-trip checks.")
            println(io,"This is validation of inputs, not physical verification, stability, convergence, or runtime estimation.")
            println(io,"manifest.tsv is the scheduler's SHA256<TAB>absolute-config-path list. config_paths.txt is not that manifest.")
            println(io,"Absolute paths assume this directory's final location. Regenerate the manifest after moving hosts/directories.")
            println(io,"Do not treat repeated identical unforced maintenance seeds as independent realizations.")
        end
        chmod(joinpath(staging,"manifest.tsv"),0o440)
        mv(staging,outdir;force=false)
    finally
        isdir(staging) && rm(staging;recursive=true)
    end
    manifest=joinpath(outdir,"manifest.tsv")
    return (; directory=outdir, paths, count=length(paths), manifest,
            manifest_sha256=bytes2hex(sha256(read(manifest))))
end

function main(args=ARGS)
    if args == ["--help"] || args == ["-h"]
        println("Usage: julia scripts/make_suite.jl PILOT.toml NEW_DIRECTORY [--seeds 42:44] [--include-emergence]")
        println("Default: 54 configs = 3 physics × 3 forcing × 2 initial families × 3 seeds.")
        println("--include-emergence adds 18 forced-only, initially quiescent upper-layer cases for the default seeds.")
        println("Reuses the package's stdlib-only input validator. Writes a hash manifest; never launches runs.")
        return 0
    end
    length(args) >= 2 || throw(ArgumentError("Expected PILOT.toml NEW_DIRECTORY; use --help"))
    seeds=collect(42:44)
    include_emergence=false
    seen=Set{String}()
    i=3
    while i <= length(args)
        option=args[i]
        option in seen && throw(ArgumentError("Duplicate option $option"))
        push!(seen,option)
        if option == "--include-emergence"
            include_emergence=true
            i+=1
        elseif option == "--seeds"
            i < length(args) || throw(ArgumentError("--seeds needs FIRST:LAST or a comma-separated list"))
            seeds=_seeds(args[i+1])
            i+=2
        else
            throw(ArgumentError("Unknown option $option"))
        end
    end
    result=make_suite(args[1],args[2];seeds,include_emergence)
    println("Generated ",result.count," validated configs in ",result.directory)
    println("Manifest: ",result.manifest)
    println("Manifest SHA256: ",result.manifest_sha256)
    println("No simulations or jobs launched. Scheduler array range, after review: 1-",result.count)
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(SuiteGenerator.main())
    catch err
        println(stderr,"Suite generation failed: ",sprint(showerror,err))
        exit(1)
    end
end
