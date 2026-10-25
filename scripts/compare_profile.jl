#!/usr/bin/env julia

module ProfileComparison

include(joinpath(@__DIR__, "..", "src", "PolarDiagnostics.jl"))
using .PolarDiagnostics: load_wind_profile
using Statistics: mean

export compare_profiles, write_summary, main

"""
    compare_profiles(observation_csv, simulation_csv)

Interpolate a simulated azimuthal wind profile onto observed radii within the
simulated interval. Both files use the explicit `radius_m,u_theta_m_s` contract.
No extrapolation or frame conversion occurs. RMSE and bias are unweighted over
the compared observed samples; bias is simulation minus observation. Partial
overlap is reported by `coverage`, never silently treated as full coverage.
"""
function compare_profiles(observation_csv::AbstractString, simulation_csv::AbstractString)
    observed = load_wind_profile(observation_csv)
    simulated = load_wind_profile(simulation_csv)
    keep = findall(r -> simulated.radius_min_m <= r <= simulated.radius_max_m,
                   observed.radius_m)
    isempty(keep) && throw(ArgumentError("No observed radius lies within the simulated profile; extrapolation is disabled"))
    radius_m = observed.radius_m[keep]
    observed_u_theta_m_s = observed.u_theta_m_s[keep]
    simulated_u_theta_m_s = zeros(length(keep))
    for (k,r) in enumerate(radius_m)
        j = searchsortedlast(simulated.radius_m,r)
        if j == length(simulated.radius_m)
            simulated_u_theta_m_s[k] = simulated.u_theta_m_s[j]
        else
            fraction = (r-simulated.radius_m[j]) / (simulated.radius_m[j+1]-simulated.radius_m[j])
            simulated_u_theta_m_s[k] = simulated.u_theta_m_s[j] + fraction*(simulated.u_theta_m_s[j+1]-simulated.u_theta_m_s[j])
        end
    end
    residual_m_s = simulated_u_theta_m_s .- observed_u_theta_m_s
    n_observed, n_compared = length(observed.radius_m), length(keep)
    coverage = n_compared/n_observed
    rmse_m_s, bias_m_s = sqrt(mean(abs2,residual_m_s)), mean(residual_m_s)
    return (; n_observed, n_compared, coverage, rmse_m_s, bias_m_s,
            radius_min_m=minimum(radius_m), radius_max_m=maximum(radius_m),
            radius_m, observed_u_theta_m_s, simulated_u_theta_m_s, residual_m_s)
end

function write_summary(io::IO, result)
    println(io,"n_observed,n_compared,coverage,rmse_m_s,bias_m_s,radius_min_m,radius_max_m")
    println(io,join((result.n_observed,result.n_compared,result.coverage,
                     result.rmse_m_s,result.bias_m_s,result.radius_min_m,result.radius_max_m),","))
end

function main(args=ARGS)
    if args == ["--help"] || args == ["-h"]
        println("Usage: julia scripts/compare_profile.jl OBSERVED.csv SIMULATED.csv [SUMMARY.csv]")
        println("Both inputs: radius_m,u_theta_m_s. Inputs must already share units, geometry and reference frame.")
        println("Reports unweighted RMSE and bias (simulation minus observation) over overlapping observed radii.")
        println("No extrapolation. Coverage is n_compared/n_observed. Output defaults to stdout.")
        return 0
    end
    length(args) in (2,3) || throw(ArgumentError("Expected OBSERVED.csv SIMULATED.csv [SUMMARY.csv]; use --help"))
    if length(args) == 3
        output_path = abspath(args[3])
        output_path in abspath.(args[1:2]) && throw(ArgumentError("Output path must differ from both input profiles"))
    end
    result = compare_profiles(args[1],args[2])
    if length(args) == 3
        mkpath(dirname(abspath(args[3])))
        open(io -> write_summary(io,result), args[3], "w")
    else
        write_summary(stdout,result)
    end
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(ProfileComparison.main())
    catch err
        println(stderr,"Profile comparison failed: ",sprint(showerror,err))
        exit(1)
    end
end
