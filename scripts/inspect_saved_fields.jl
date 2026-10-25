using SaturnPolarQG, JLD2, FFTW, Statistics, LinearAlgebra, TOML, SHA

length(ARGS) == 2 || error("usage: inspect_saved_fields.jl RESULTS ANALYSIS_DIRECTORY")
result, out = abspath.(ARGS)
result == out && error("analysis must be separate from saved results")
mkpath(out)
c = load_config(joinpath(result, "config.toml"))
meta = TOML.parsefile(joinpath(result, "metadata.toml"))
meta["source_sha256"] == SaturnPolarQG.source_digest() || error("analysis source differs from saved run")
bytes2hex(sha256(read(joinpath(result,"Manifest.toml")))) == meta["manifest_sha256"] || error("saved manifest hash differs")
domain = load(joinpath(result, "domain.jld2"))
x, y = domain["x"], domain["y"]
n = length(x); L = c["grid"]["length_m"]; kap2 = c["physics"]["deformation_radius_m"]^-2
kr = reshape((2pi/L).*(0:n÷2), :, 1)
ell = reshape((2pi/L).*vcat(0:n÷2-1, -n÷2:-1), 1, :)
ksq = kr.^2 .+ ell.^2
summary = Dict{String,Any}("source_matches"=>true, "snapshots"=>Any[])
snapshots = sort(filter(s -> startswith(s,"snapshot_") && endswith(s,".jld2"), readdir(result)))
traces = []
for file in snapshots
    f = load(joinpath(result,file))
    for key in ["psi","s","u","v","totalpv"]
        size(f[key]) == (n,n) || error("wrong shape: $file/$key")
        all(isfinite,f[key]) || error("nonfinite field: $file/$key")
    end
    ph = rfft(f["psi"])
    expected = Dict("u"=>irfft(-im.*ell.*ph,n), "v"=>irfft(im.*kr.*ph,n),
                    "s"=>irfft(-(ksq.+kap2).*ph,n))
    errors = Dict(k=>norm(f[k]-a)/max(norm(a),eps()) for (k,a) in expected)
    pv_error = maximum(abs,f["totalpv"]-f["s"]-domain["background_pv"])
    options = (;jet_radius_m=c["initial"]["jet_radius_m"],
        jet_halfwidth_m=c["output"]["jet_search_halfwidth_m"],m_max=c["output"]["m_max"])
    polar(nr,ntheta) = SaturnPolarQG.PolarDiagnostics.polar_diagnostics(
        f["u"],f["v"],f["totalpv"],f["psi"],f["s"],x,y;
        options...,nr,ntheta)
    pol = polar(c["output"]["nr"],c["output"]["ntheta"])
    fine = polar(512,720)
    push!(traces,pol)
    energy = mean(f["u"].^2+f["v"].^2+kap2.*f["psi"].^2)/2
    push!(summary["snapshots"],Dict("file"=>file,"time_s"=>f["time_s"],"step"=>f["step"],
        "all_fields_finite"=>true,"energy_m2_s2"=>energy,
        "enstrophy_s2"=>mean(abs2,f["totalpv"])/2,"mean_s_s"=>mean(f["s"]),
        "relative_inversion_errors"=>errors,"pv_identity_max_error_s"=>pv_error,
        "jet_coverage"=>pol.coverage,"m6_amplitude_m"=>pol.mode_amplitudes[6],
        "m6_finer_sampling_m"=>fine.mode_amplitudes[6],
        "finer_sampling_max_mode_difference_m"=>maximum(abs,pol.mode_amplitudes-fine.mode_amplitudes)))
end
open(joinpath(out,"jet_traces.csv"),"w") do io
    println(io,"theta_rad,"*join(["radius_step_"*string(s["step"])*"_m" for s in summary["snapshots"]],","))
    for j in eachindex(traces[1].theta)
        println(io,join(vcat(traces[1].theta[j],[p.jet_radius_theta[j] for p in traces]),","))
    end
end
open(joinpath(out,"radial_profiles.csv"),"w") do io
    println(io,"radius_m,"*join(["utheta_step_"*string(s["step"])*"_m_s" for s in summary["snapshots"]],","))
    for j in eachindex(traces[1].radii)
        println(io,join(vcat(traces[1].radii[j],[p.mean_utheta[j] for p in traces]),","))
    end
end
last = load(joinpath(result,snapshots[end]))
checkpoint = load(joinpath(result,"checkpoint.jld2"))
summary["checkpoint"] = Dict("time_s"=>checkpoint["time_s"],"step"=>checkpoint["step"],
    "attempted_events"=>checkpoint["convection"].events,
    "next_event_s"=>checkpoint["convection"].next_event,
    "snapshot_s_max_abs_difference_s"=>maximum(abs,irfft(checkpoint["sol"],n)-last["s"]))
open(joinpath(out,"field_checks.toml"),"w") do io
    TOML.print(io,summary)
end
println("Checked ",length(snapshots)," saved snapshots and final checkpoint; wrote ",out)
