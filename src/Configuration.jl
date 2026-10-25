const CONFIG_KEYS = Dict(
    "experiment" => ["name", "purpose", "seed"],
    "grid" => ["nx", "length_m"],
    "planet" => ["radius_m", "rotation_rate_s"],
    "physics" => ["deformation_radius_m", "deep_jet_speed_m_s", "deep_jet_radius_m", "deep_jet_width_m"],
    "initial" => ["jet_speed_m_s", "jet_radius_m", "jet_width_m", "cyclone_speed_m_s", "cyclone_radius_m", "perturbation", "displacement_m", "mode", "m_max"],
    "boundary" => ["inner_radius_m", "outer_radius_m", "sponge_rate_s"],
    "numerics" => ["dt_s", "duration_s", "hyper_order", "cutoff_damping_time_s", "cfl_limit"],
    "convection" => ["enabled", "r_min_m", "r_max_m", "vortex_radius_m", "q_max_s", "bias_max", "psi_limit_m2_s", "interval_min_s", "interval_max_s"],
    "output" => ["diagnostic_interval_s", "snapshot_interval_s", "checkpoint_interval_s", "m_max", "ntheta", "nr", "jet_search_halfwidth_m"])

function load_config(path)
    c = TOML.parsefile(path)
    validate_config(c)
    return c
end

function validate_config(c)
    Set(keys(c)) == Set(keys(CONFIG_KEYS)) || error("Configuration sections must be exactly $(sort(collect(keys(CONFIG_KEYS))))")
    for (section, allowed) in CONFIG_KEYS
        v = c[section]
        v isa AbstractDict || error("$section must be a table")
        Set(keys(v)) == Set(allowed) || error("Keys in [$section] must be exactly $allowed; unknown/missing fields are rejected")
        for (k, x) in v
            if x isa Real && !(x isa Bool)
                (isfinite(x) || (section == "physics" && k == "deformation_radius_m" && x == Inf) ||
                 (section == "numerics" && k == "cutoff_damping_time_s" && x == Inf)) || error("$section.$k must be finite")
            end
        end
    end
    positive(section, ks) = foreach(k -> c[section][k] isa Real && c[section][k] > 0 || error("$section.$k must be positive"), ks)
    positive("grid", ["nx", "length_m"])
    n = c["grid"]["nx"]; n isa Integer && n >= 32 && iseven(n) || error("grid.nx must be an even integer >= 32")
    positive("planet", ["radius_m", "rotation_rate_s"])
    positive("physics", ["deformation_radius_m", "deep_jet_radius_m", "deep_jet_width_m"])
    positive("initial", ["jet_radius_m", "jet_width_m", "cyclone_radius_m"])
    c["initial"]["displacement_m"] >= 0 || error("initial displacement must be nonnegative")
    c["initial"]["perturbation"] in ("none", "mode", "broadband") || error("initial.perturbation must be none, mode, or broadband")
    for k in ("mode", "m_max")
        c["initial"][k] isa Integer && 1 <= c["initial"][k] <= 32 || error("initial.$k must be an integer in 1:32")
    end
    c["experiment"]["seed"] isa Integer && c["experiment"]["seed"] >= 0 || error("seed must be a nonnegative integer")
    c["experiment"]["purpose"] in ("smoke", "maintenance", "selection", "emergence") || error("unknown experiment purpose")
    positive("numerics", ["dt_s", "duration_s", "cutoff_damping_time_s", "cfl_limit"])
    c["numerics"]["cfl_limit"] <= 0.8 || error("cfl_limit must be <= 0.8")
    c["numerics"]["hyper_order"] isa Integer && 1 <= c["numerics"]["hyper_order"] <= 8 || error("hyper_order must be in 1:8")
    positive("output", ["diagnostic_interval_s", "snapshot_interval_s", "checkpoint_interval_s", "m_max", "ntheta", "nr", "jet_search_halfwidth_m"])
    c["output"]["m_max"] isa Integer && 6 <= c["output"]["m_max"] <= 32 || error("output.m_max must be 6:32")
    for k in ("ntheta", "nr")
        c["output"][k] isa Integer || error("output.$k must be an integer")
    end
    c["output"]["ntheta"] >= 4c["output"]["m_max"] || error("ntheta must resolve diagnostic modes")
    c["output"]["nr"] >= 8 || error("nr must be at least 8")
    b = c["boundary"]; L = c["grid"]["length_m"]; dx = L/n
    0 < b["inner_radius_m"] < b["outer_radius_m"] < L/2 || error("boundary radii must satisfy 0 < inner < outer < L/2")
    b["outer_radius_m"] - b["inner_radius_m"] >= 4dx || error("boundary transition needs at least four grid cells")
    b["sponge_rate_s"] >= 0 || error("sponge rate must be nonnegative")
    c["initial"]["jet_radius_m"] + c["output"]["jet_search_halfwidth_m"] < b["inner_radius_m"] || error("jet search annulus overlaps artificial boundary")
    c["initial"]["jet_width_m"] >= 3dx || error("jet width needs at least three grid cells (convergence still required)")
    c["physics"]["deep_jet_speed_m_s"] == 0 || c["physics"]["deep_jet_width_m"] >= 3dx || error("deep jet width needs at least three grid cells")
    c["initial"]["cyclone_speed_m_s"] == 0 || c["initial"]["cyclone_radius_m"] >= 3dx || error("cyclone radius needs at least three grid cells")
    imax = c["initial"]["perturbation"] == "mode" ? c["initial"]["mode"] : c["initial"]["m_max"]
    c["initial"]["perturbation"] == "none" || 2pi*c["initial"]["jet_radius_m"]/imax >= 6dx || error("initial azimuthal modes require at least six cells per wavelength")
    if c["convection"]["enabled"]
        f = c["convection"]; ConvectionState(f, c["experiment"]["seed"])
        f["vortex_radius_m"] >= 4dx || error("convective vortex radius needs at least four grid cells; increase nx, not the physical radius")
        f["r_max_m"] + 2f["vortex_radius_m"] + dx < b["inner_radius_m"] || error("convective pairs can enter artificial boundary")
    else
        c["convection"]["enabled"] isa Bool || error("convection.enabled must be Boolean")
    end
    return c
end
