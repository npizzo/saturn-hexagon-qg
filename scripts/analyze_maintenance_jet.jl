"""
Remeasure saved maintenance snapshots without advancing or modifying the model.

Usage:
  julia --project=. scripts/analyze_maintenance_jet.jl RESULTS [RESULTS...] OUTPUT

Every snapshot is measured on its configured annulus and a wider 7–23 Mm
annulus, using nr=512 and ntheta=720. Final snapshots also use 1024/1440.
The wider search is a sensitivity test, not an assumption of a coherent jet.
"""
module MaintenanceJetAudit

using JLD2, SHA, TOML, Dates, Statistics

source_digest(directory) = bytes2hex(sha256(join(
    [name*"\n"*read(joinpath(directory,name),String)
     for name in sort(filter(name -> endswith(name,".jl"),readdir(directory)))], "\n")))
csvcell(x::AbstractString) = "\""*replace(x,"\""=>"\"\"")*"\""
csvcell(x) = string(x)
csvrow(io,row) = println(io,join(csvcell.(row),','))

function snapshot_records(path)
    records = NamedTuple[]
    for name in readdir(path)
        occursin(r"^snapshot_[0-9]+\.jld2$",name) || continue
        filename = joinpath(path,name)
        jldopen(filename,"r") do file
            time = Float64(file["time_s"])
            isfinite(time) || error("Nonfinite snapshot time: $filename")
            push!(records,(;time,step=file["step"],path=filename))
        end
    end
    isempty(records) && error("No snapshots in $path")
    sort!(records;by=record -> record.time)
    all(diff(getproperty.(records,:time)) .> 0) || error("Snapshot times must be distinct")
    return records
end

function read_run(path)
    config = TOML.parsefile(joinpath(path,"config.toml"))
    metadata = TOML.parsefile(joinpath(path,"metadata.toml"))
    source_digest(joinpath(path,"src")) == metadata["source_sha256"] || error("Archived source hash mismatch in $path")
    bytes2hex(sha256(read(joinpath(path,"Manifest.toml")))) == metadata["manifest_sha256"] || error("Archived Manifest hash mismatch in $path")
    return (;path,id=basename(dirname(path)),config,metadata,
        snapshots=snapshot_records(path),domain=JLD2.load(joinpath(path,"domain.jld2")))
end

function mean_profile_width(r,u)
    all(isfinite,u) || return (;radius=NaN,speed=NaN,fwhm=NaN,closed=false)
    k = argmax(u); peak = u[k]
    peak > 0 || return (;radius=r[k],speed=peak,fwhm=NaN,closed=false)
    half = peak/2
    il = findlast(i -> u[i] <= half,1:k)
    ir_relative = findfirst(i -> u[i] <= half,k:length(u))
    ir = isnothing(ir_relative) ? nothing : k+ir_relative-1
    if isnothing(il) || isnothing(ir) || il == k || ir == k
        return (;radius=r[k],speed=peak,fwhm=NaN,closed=false)
    end
    left = r[il]+(half-u[il])*(r[il+1]-r[il])/(u[il+1]-u[il])
    right = r[ir-1]+(half-u[ir-1])*(r[ir]-r[ir-1])/(u[ir]-u[ir-1])
    return (;radius=r[k],speed=peak,fwhm=right-left,closed=true)
end

"Extra raw-peak information; these thresholds identify sensitivity, not significance."
function peak_details(sampler,grid,u,v,pol,separation,large_jump)
    ntheta,nr = length(pol.theta),length(pol.radii)
    raw_radius,raw_speed,secondary_radius,secondary_speed,ratios,valleys =
        [fill(NaN,ntheta) for _ in 1:6]
    at_inner,at_outer = falses(ntheta),falses(ntheta)
    peak_counts = zeros(Int,ntheta)
    profile = zeros(nr)
    for k in eachindex(pol.theta)
        st,ct = sincos(pol.theta[k])
        for j in eachindex(pol.radii)
            xx,yy = pol.radii[j]*ct,pol.radii[j]*st
            profile[j] = -st*sampler(u,xx,yy,grid)+ct*sampler(v,xx,yy,grid)
        end
        all(isfinite,profile) || continue
        p = argmax(profile)
        raw_radius[k],raw_speed[k] = pol.radii[p],profile[p]
        at_inner[k],at_outer[k] = p == 1,p == nr
        peaks = [j for j in 2:nr-1 if profile[j] > 0 && profile[j] > profile[j-1] && profile[j] >= profile[j+1]]
        peak_counts[k] = length(peaks)
        candidates = [j for j in peaks if abs(pol.radii[j]-pol.radii[p]) >= separation]
        isempty(candidates) && continue
        q = candidates[argmax(profile[candidates])]
        secondary_radius[k],secondary_speed[k] = pol.radii[q],profile[q]
        if profile[p] > 0
            ratios[k] = profile[q]/profile[p]
            valleys[k] = (min(profile[p],profile[q])-minimum(view(profile,min(p,q):max(p,q))))/profile[p]
        end
    end
    jump_to_next = fill(NaN,ntheta)
    for k in 1:ntheta
        kp = mod1(k+1,ntheta)
        if isfinite(pol.jet_radius_theta[k]) && isfinite(pol.jet_radius_theta[kp])
            jump_to_next[k] = abs(pol.jet_radius_theta[kp]-pol.jet_radius_theta[k])
        end
    end
    valid_jumps = filter(isfinite,jump_to_next)
    return (;raw_radius,raw_speed,secondary_radius,secondary_speed,ratios,valleys,
        at_inner,at_outer,peak_counts,jump_to_next,
        max_jump=isempty(valid_jumps) ? NaN : maximum(valid_jumps),
        large_jump_fraction=isempty(valid_jumps) ? NaN : count(>(large_jump),valid_jumps)/length(valid_jumps))
end

function analyze_one(diag,run,snapshot,fields,nr,ntheta,scope,summary_io,traces_io,profiles_io)
    cfg = run.config
    if scope == "configured"
        radius,halfwidth = cfg["initial"]["jet_radius_m"],cfg["output"]["jet_search_halfwidth_m"]
    else
        radius,halfwidth = 15e6,8e6 # Explicitly labeled fixed 7–23 Mm sensitivity band.
    end
    pol = Base.invokelatest(getfield(diag,:polar_diagnostics),
        fields["u"],fields["v"],fields["totalpv"],fields["psi"],fields["s"],run.domain["x"],run.domain["y"];
        jet_radius_m=radius,jet_halfwidth_m=halfwidth,m_max=12,nr,ntheta)
    grid = Base.invokelatest(getfield(diag,:_grid),run.domain["x"],run.domain["y"],size(fields["u"]))
    separation = max(4max(grid.dx,grid.dy),0.5cfg["initial"]["jet_width_m"])
    detail = Base.invokelatest(peak_details,getfield(diag,:_sample_cubic),grid,
        fields["u"],fields["v"],pol,separation,cfg["initial"]["jet_width_m"])
    profile = mean_profile_width(pol.radii,pol.mean_utheta)
    configured_lo = cfg["initial"]["jet_radius_m"]-cfg["output"]["jet_search_halfwidth_m"]
    configured_hi = cfg["initial"]["jet_radius_m"]+cfg["output"]["jet_search_halfwidth_m"]
    good_radius = filter(isfinite,pol.jet_radius_theta)
    outside_fraction = isempty(good_radius) ? NaN : count(r -> r < configured_lo || r > configured_hi,good_radius)/length(good_radius)
    dominant = pol.jet_detected ? argmax(pol.mode_amplitudes) : NaN
    dominant_amplitude = pol.jet_detected ? maximum(pol.mode_amplitudes) : NaN
    common = Any[run.id,cfg["convection"]["enabled"],basename(snapshot.path),snapshot.time,snapshot.step,scope,nr,ntheta]
    row = Any[common...,
        pol.coverage,pol.jet_detected,pol.meanjet,pol.mean_peak_speed,
        profile.radius,profile.speed,profile.fwhm,profile.closed,
        mean(detail.at_inner),mean(detail.at_outer),mean(detail.peak_counts .> 1),
        count(r -> isfinite(r) && r >= .9,detail.ratios)/ntheta,
        count(k -> isfinite(detail.ratios[k]) && detail.ratios[k] >= .9 && detail.valleys[k] >= .05,1:ntheta)/ntheta,
        detail.max_jump,detail.large_jump_fraction,outside_fraction,
        dominant,dominant_amplitude,angle(pol.mode_coefficients[6])]
    for m in 1:12
        append!(row,(pol.mode_amplitudes[m],real(pol.mode_coefficients[m]),imag(pol.mode_coefficients[m])))
    end
    csvrow(summary_io,row)
    for k in eachindex(pol.theta)
        csvrow(traces_io,Any[common...,pol.theta[k],pol.jet_radius_theta[k],
            isfinite(pol.jet_radius_theta[k]),detail.raw_radius[k],detail.raw_speed[k],
            detail.at_inner[k],detail.at_outer[k],detail.peak_counts[k],
            detail.secondary_radius[k],detail.secondary_speed[k],detail.ratios[k],
            detail.valleys[k],detail.jump_to_next[k]])
    end
    for j in eachindex(pol.radii)
        csvrow(profiles_io,Any[common...,pol.radii[j],pol.mean_utheta[j],pol.mean_ur[j],
            pol.mean_pv[j],pol.eddy_uv[j],pol.radial_coverage[j]])
    end
    flush(summary_io);flush(traces_io);flush(profiles_io)
    println((;run=run.id,time=snapshot.time,scope,nr,ntheta,coverage=pol.coverage,
        A6_m=pol.mode_amplitudes[6],dominant,mean_profile_peak_speed=profile.speed,
        mean_profile_fwhm_m=profile.fwhm,boundary_fraction=mean(detail.at_inner .| detail.at_outer),
        ambiguous_fraction=count(k -> isfinite(detail.ratios[k]) && detail.ratios[k] >= .9 && detail.valleys[k] >= .05,1:ntheta)/ntheta,
        max_radius_jump_m=detail.max_jump))
end

function main(args)
    length(args) >= 2 || error("Usage: analyze_maintenance_jet.jl RESULTS [RESULTS...] OUTPUT")
    paths,output = abspath.(args[1:end-1]),abspath(last(args))
    length(unique(paths)) == length(paths) || error("Duplicate input directories")
    runs = read_run.(paths)
    length(unique(getproperty.(runs,:id))) == length(runs) || error("Run directory basenames must be unique")
    all(run -> run.metadata["source_sha256"] == first(runs).metadata["source_sha256"],runs) || error("Saved model sources differ")
    destinations = [joinpath(output,name) for name in
        ("polar_cubic.csv","jettraces.csv","profiles.csv","polar_cubic_metadata.toml")]
    any(ispath,destinations) && error("Audit output files already exist; choose another output directory")
    Base.include(@__MODULE__,joinpath(first(paths),"src","PolarDiagnostics.jl"))
    diag = getfield(@__MODULE__,:PolarDiagnostics)
    Core.eval(diag,:(_sample(field,x,y,grid) = _sample_cubic(field,x,y,grid)))
    common_header = ["run_id","forcing_enabled","snapshot","time_s","step","scope","nr","ntheta"]
    summary_header = [common_header...,
        "coverage","jet_detected","conditional_mean_radius_m","conditional_mean_peak_speed_m_s",
        "mean_profile_peak_radius_m","mean_profile_peak_speed_m_s","mean_profile_fwhm_m","mean_profile_fwhm_closed",
        "raw_max_at_inner_fraction","raw_max_at_outer_fraction","multiple_positive_local_maxima_fraction",
        "strong_separated_secondary_fraction","ambiguous_separated_peaks_fraction",
        "max_adjacent_valid_radius_jump_m","large_adjacent_jump_fraction","trace_outside_configured_annulus_fraction",
        "dominant_mode","dominant_amplitude_m","m6_coefficient_phase_rad"]
    for m in 1:12
        append!(summary_header,["m$(m)_amplitude_m","m$(m)_real_m","m$(m)_imag_m"])
    end
    traces_header = [common_header...,"theta_rad","jet_radius_m","valid_peak","raw_max_radius_m","raw_max_speed_m_s",
        "raw_max_at_inner","raw_max_at_outer","positive_local_maxima_count","secondary_peak_radius_m",
        "secondary_peak_speed_m_s","secondary_to_primary_speed_ratio","interpeak_valley_drop_fraction","jump_to_next_valid_radius_m"]
    profiles_header = [common_header...,"radius_m","mean_utheta_m_s","mean_ur_m_s","mean_totalpv_s",
        "eddy_ur_utheta_m2_s2","velocity_coverage"]
    mkpath(output)
    open(destinations[1],"w") do summary_io
        open(destinations[2],"w") do traces_io
            open(destinations[3],"w") do profiles_io
                csvrow(summary_io,summary_header);csvrow(traces_io,traces_header);csvrow(profiles_io,profiles_header)
                for run in runs
                    for (index,snapshot) in enumerate(run.snapshots)
                        fields = JLD2.load(snapshot.path)
                        for scope in ("configured","wide_sensitivity")
                            analyze_one(diag,run,snapshot,fields,512,720,scope,summary_io,traces_io,profiles_io)
                            if index == length(run.snapshots)
                                analyze_one(diag,run,snapshot,fields,1024,1440,scope,summary_io,traces_io,profiles_io)
                            end
                        end
                    end
                end
            end
        end
    end
    for run in runs
        source_digest(joinpath(run.path,"src")) == run.metadata["source_sha256"] || error("Archived source changed during audit")
    end
    provenance = Dict("created_utc"=>string(now(UTC)),"audit_julia_version"=>string(VERSION),
        "script"=>abspath(@__FILE__),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "source_sha256"=>first(runs).metadata["source_sha256"],
        "interpolation"=>"Archived diagnostics with process-local four-point tensor cubic Lagrange sampler. Source files and simulations are unchanged.",
        "sequence_sampling"=>"nr512/ntheta720 for both configured and wide annuli at every saved actual time; both final annuli repeated at nr1024/ntheta1440.",
        "wide_search_m"=>[7e6,23e6],
        "secondary_peak_definition"=>"Positive interior radial local maximum separated from the raw global maximum by at least max(4 grid cells, half the initialized jet width). Strong secondary means speed ratio >=0.9. Ambiguous separated peaks additionally require the interpeak minimum to be >=0.05 primary peak speed below the smaller peak. These are explicit sensitivity flags, not physical classification or statistical confidence thresholds.",
        "large_jump_definition"=>"Adjacent finite traced radii differing by more than the initialized physical jet width; gaps are never bridged. Angular resolution changes the jump diagnostic.",
        "mean_profile_width"=>"FWHM around the maximum of the full azimuthal-mean tangential velocity profile, using nearest half-height crossings with linear radial interpolation. NaN when either crossing is outside the sampled interval. Multiple peaks can make this width incomplete as a description of the flow.",
        "coverage_warning"=>"Conditional mean traced radius and peak speed use only valid sectors and are biased as global jet measures when coverage is incomplete. Full-circle radial means/profiles are recorded separately. Broader search and full coverage do not establish that maxima belong to a single connected jet.",
        "phase_warning"=>"Mode coefficients and phase remain NaN when the archived diagnostic rejects any angular sector. No missing sectors are filled and no phase is unwrapped across gaps.",
        "runs"=>Dict(run.id=>Dict("results"=>run.path,"forcing_enabled"=>run.config["convection"]["enabled"],
            "manifest_sha256"=>run.metadata["manifest_sha256"],"run_julia_version"=>run.metadata["julia_version"],
            "snapshot_files"=>[s.path for s in run.snapshots],"snapshot_times_s"=>[s.time for s in run.snapshots]) for run in runs))
    open(destinations[4],"w") do io;TOML.print(io,provenance);end
    return nothing
end

end # module

abspath(PROGRAM_FILE) == abspath(@__FILE__) && MaintenanceJetAudit.main(ARGS)
