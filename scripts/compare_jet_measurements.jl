"""
Offline comparison of jet diagnostics for two runs differing only in timestep.

Usage:
    julia --project=. scripts/compare_jet_measurements.jl BASE_RESULTS REFINED_RESULTS ANALYSIS_OUTPUT

Only matching initial/final snapshot times are compared. Archived diagnostics
are loaded from BASE_RESULTS/src. Cubic interpolation replaces a method only
inside this Julia process; neither production source nor run files are edited.
"""
module JetMeasurementComparison

using JLD2, SHA, TOML, Dates

source_digest(directory) = bytes2hex(sha256(join(
    [name*"\n"*read(joinpath(directory,name),String)
     for name in sort(filter(name -> endswith(name,".jl"),readdir(directory)))], "\n")))

function snapshot_endpoints(path)
    records = NamedTuple[]
    for name in readdir(path)
        occursin(r"^snapshot_[0-9]+\.jld2$",name) || continue
        filename = joinpath(path,name)
        jldopen(filename,"r") do file
            time = file["time_s"]
            isfinite(time) || error("Nonfinite snapshot time: $filename")
            push!(records,(;time=Float64(time),step=file["step"],path=filename))
        end
    end
    length(records) >= 2 || error("At least two snapshots are required in $path")
    sort!(records;by=record -> record.time)
    first(records).time < last(records).time || error("Snapshots need distinct initial/final times")
    return first(records),last(records)
end

function read_run(path)
    config = TOML.parsefile(joinpath(path,"config.toml"))
    metadata = TOML.parsefile(joinpath(path,"metadata.toml"))
    source_digest(joinpath(path,"src")) == metadata["source_sha256"] || error("Archived source hash mismatch in $path")
    bytes2hex(sha256(read(joinpath(path,"Manifest.toml")))) == metadata["manifest_sha256"] || error("Archived Manifest hash mismatch in $path")
    initial,final = snapshot_endpoints(path)
    return (;path,config,metadata,initial,final,domain=JLD2.load(joinpath(path,"domain.jld2")))
end

function write_csv(path,header,rows)
    # All output string values are controlled identifiers; no CSV quoting is
    # needed. Paths and descriptive text are stored separately in TOML.
    open(path,"w") do io
        println(io,join(header,','))
        foreach(row -> println(io,join(row,',')),rows)
    end
end

function measure(diagnostic_module,fields,domain,config,nr,ntheta)
    Base.invokelatest(getfield(diagnostic_module,:polar_diagnostics),
        fields["u"],fields["v"],fields["totalpv"],fields["psi"],fields["s"],domain["x"],domain["y"];
        jet_radius_m=config["initial"]["jet_radius_m"],
        jet_halfwidth_m=config["output"]["jet_search_halfwidth_m"],
        m_max=config["output"]["m_max"],nr,ntheta)
end

function main(args)
    length(args) == 3 || error("Usage: compare_jet_measurements.jl BASE_RESULTS REFINED_RESULTS ANALYSIS_OUTPUT")
    basepath,finepath,output = abspath.(args)
    basepath != finepath || error("Provide distinct result directories")
    base,fine = read_run(basepath),read_run(finepath)
    base.metadata["source_sha256"] == fine.metadata["source_sha256"] || error("Run sources differ")
    base.metadata["manifest_sha256"] == fine.metadata["manifest_sha256"] || error("Run environments differ")
    base.metadata["julia_version"] == fine.metadata["julia_version"] || error("Run Julia versions differ")
    base.metadata["device"] == fine.metadata["device"] || error("Run devices differ")
    original,revised = deepcopy(base.config),deepcopy(fine.config)
    base_dt = pop!(original["numerics"],"dt_s")
    fine_dt = pop!(revised["numerics"],"dt_s")
    original == revised || error("Configurations must differ only in numerics.dt_s")
    0 < fine_dt < base_dt || error("REFINED_RESULTS must use a smaller positive dt_s")
    base.initial.time == fine.initial.time || error("Initial snapshot times differ; no temporal interpolation is performed")
    base.final.time == fine.final.time || error("Final snapshot times differ; no temporal interpolation is performed")
    for key in ("x","y","background_pv","psi_deep","reference_s")
        isequal(base.domain[key],fine.domain[key]) || error("Saved domains differ in $key")
    end
    destinations = [joinpath(output,name) for name in
        ("jet_measurements.csv","jet_measurement_modes.csv","jet_measurement_metadata.toml")]
    any(ispath,destinations) && error("Audit output already exists; choose a new analysis directory or remove only the previous audit files deliberately")

    # Use the exact saved diagnostic implementation, independent of later
    # edits elsewhere in the project. All overrides below are process-local.
    Base.include(@__MODULE__,joinpath(basepath,"src","PolarDiagnostics.jl"))
    diagnostic_module = getfield(@__MODULE__,:PolarDiagnostics)
    isdefined(diagnostic_module,:_sample_cubic) || error("Archived diagnostics lack the cubic sampler")
    initial_base,initial_fine = JLD2.load(base.initial.path),JLD2.load(fine.initial.path)
    final_base,final_fine = JLD2.load(base.final.path),JLD2.load(fine.final.path)
    snapshots = [("initial",base.initial.time,initial_base,initial_fine),
                 ("final",base.final.time,final_base,final_fine)]
    rows,mode_rows = Vector{Vector{Any}}(),Vector{Vector{Any}}()
    settings = [("bilinear",128,360),("cubic_lagrange",128,360),
                ("cubic_lagrange",512,720),("cubic_lagrange",1024,1440)]
    cubic_installed = false
    for (interpolation,nr,ntheta) in settings
        if interpolation == "cubic_lagrange" && !cubic_installed
            Core.eval(diagnostic_module,:(_sample(field,x,y,grid) = _sample_cubic(field,x,y,grid)))
            cubic_installed = true
        end
        for (phase,time,basefields,finefields) in snapshots
            a = measure(diagnostic_module,basefields,base.domain,base.config,nr,ntheta)
            b = measure(diagnostic_module,finefields,fine.domain,fine.config,nr,ntheta)
            delta_phase = angle(b.mode_coefficients[6]*conj(a.mode_coefficients[6]))
            delta_amp = b.mode_amplitudes[6]-a.mode_amplitudes[6]
            other = [m for m in eachindex(a.mode_amplitudes) if m != 6]
            push!(rows,Any[phase,time,interpolation,nr,ntheta,base_dt,fine_dt,
                a.mode_amplitudes[6],b.mode_amplitudes[6],delta_amp,
                100delta_amp/a.mode_amplitudes[6],angle(a.mode_coefficients[6]),
                angle(b.mode_coefficients[6]),delta_phase,-delta_phase/6,
                a.meanjet,b.meanjet,b.meanjet-a.meanjet,a.mean_peak_speed,b.mean_peak_speed,
                a.coverage,b.coverage,maximum(abs.(b.mode_amplitudes[other].-a.mode_amplitudes[other]))])
            for m in eachindex(a.mode_amplitudes)
                push!(mode_rows,Any[phase,time,interpolation,nr,ntheta,m,
                    a.mode_amplitudes[m],b.mode_amplitudes[m],b.mode_amplitudes[m]-a.mode_amplitudes[m],
                    real(a.mode_coefficients[m]),imag(a.mode_coefficients[m]),
                    real(b.mode_coefficients[m]),imag(b.mode_coefficients[m])])
            end
            println((;phase,time,interpolation,nr,ntheta,A6_base_m=a.mode_amplitudes[6],
                A6_refined_m=b.mode_amplitudes[6],delta_A6_m=delta_amp,
                delta_coefficient_phase_rad=delta_phase,delta_mean_radius_m=b.meanjet-a.meanjet))
        end
    end
    for run in (base,fine)
        source_digest(joinpath(run.path,"src")) == run.metadata["source_sha256"] || error("Archived source changed during audit")
    end
    mkpath(output)
    write_csv(destinations[1],["phase","time_s","interpolation","nr","ntheta","base_dt_s","refined_dt_s",
        "A6_base_m","A6_refined_m","A6_refined_minus_base_m","A6_relative_difference_percent",
        "m6_phase_base_rad","m6_phase_refined_rad","m6_wrapped_phase_difference_rad",
        "pattern_orientation_difference_rad","mean_radius_base_m","mean_radius_refined_m",
        "mean_radius_refined_minus_base_m","mean_peak_speed_base_m_s","mean_peak_speed_refined_m_s",
        "coverage_base","coverage_refined","max_other_mode_amplitude_difference_m"],rows)
    write_csv(destinations[2],["phase","time_s","interpolation","nr","ntheta","mode",
        "amplitude_base_m","amplitude_refined_m","amplitude_refined_minus_base_m",
        "coefficient_base_real_m","coefficient_base_imag_m","coefficient_refined_real_m",
        "coefficient_refined_imag_m"],mode_rows)
    provenance = Dict("created_utc"=>string(now(UTC)),"base_results"=>basepath,"refined_results"=>finepath,
        "base_dt_s"=>base_dt,"refined_dt_s"=>fine_dt,"source_sha256"=>base.metadata["source_sha256"],
        "manifest_sha256"=>base.metadata["manifest_sha256"],"run_device"=>base.metadata["device"],
        "run_julia_version"=>base.metadata["julia_version"],"audit_julia_version"=>string(VERSION),
        "base_initial_snapshot"=>base.initial.path,"refined_initial_snapshot"=>fine.initial.path,
        "base_final_snapshot"=>base.final.path,"refined_final_snapshot"=>fine.final.path,
        "initial_time_s"=>base.initial.time,"final_time_s"=>base.final.time,
        "script"=>abspath(@__FILE__),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "comparison"=>"Refined minus base at identical saved times and extraction settings; intermediate snapshots are intentionally excluded.",
        "interpolation"=>"Archived bilinear implementation plus a process-local replacement with its four-point tensor cubic Lagrange sampler. No production source or run file is edited.",
        "limitations"=>"Two timesteps test short-time sensitivity, not convergence order. Polar diagnostic agreement does not establish full-field, spatial, forcing, boundary, long-time, or physical-model convergence.")
    open(destinations[3],"w") do io; TOML.print(io,provenance); end
    return nothing
end

end # module

abspath(PROGRAM_FILE) == abspath(@__FILE__) && JetMeasurementComparison.main(ARGS)
