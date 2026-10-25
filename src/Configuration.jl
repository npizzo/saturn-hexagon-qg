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

const OPTIONAL_PHYSICS_KEYS = ["deep_profile", "deep_match_fraction", "deep_jet_fraction", "deep_cyclone_fraction"]
const OPTIONAL_INITIAL_KEYS = ["profile", "core_regularization_m", "m_min"]

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
        optional = section == "physics" ? OPTIONAL_PHYSICS_KEYS : section == "initial" ? OPTIONAL_INITIAL_KEYS : String[]
        all(k -> haskey(v, k), allowed) && all(k -> k in allowed || k in optional, keys(v)) ||
            error("Keys in [$section] must include $allowed; only $optional are optional; unknown/missing fields are rejected")
        for (k, x) in v
            x isa Bool && !(section == "convection" && k == "enabled") && error("$section.$k must not be Boolean")
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
    d = c["physics"]
    profile = get(d, "deep_profile", "annular")
    profile in ("annular", "matched_upper", "matched_components") || error("physics.deep_profile must be annular, matched_upper, or matched_components")
    d["deep_jet_speed_m_s"] isa Real || error("physics.deep_jet_speed_m_s must be real")
    if profile == "matched_upper"
        d["deep_jet_speed_m_s"] == 0 || error("matched_upper requires deep_jet_speed_m_s=0; the annular deep jet is not added")
        alpha = get(d, "deep_match_fraction", 1.0)
        alpha isa Real && !(alpha isa Bool) && isfinite(alpha) && alpha >= 0 ||
            error("physics.deep_match_fraction must be a finite nonnegative real")
    elseif profile == "matched_components"
        d["deep_jet_speed_m_s"] == 0 || error("matched_components requires deep_jet_speed_m_s=0; the annular deep jet is not added")
        !haskey(d, "deep_match_fraction") || error("matched_components uses separate jet and cyclone fractions")
        for key in ("deep_jet_fraction", "deep_cyclone_fraction")
            haskey(d,key) || error("matched_components requires physics.$key")
            alpha = d[key]
            alpha isa Real && !(alpha isa Bool) && isfinite(alpha) && alpha >= 0 ||
                error("physics.$key must be a finite nonnegative real")
        end
    else
        !haskey(d, "deep_match_fraction") || error("deep_match_fraction is only valid for deep_profile=matched_upper")
    end
    if profile != "matched_components"
        any(k -> haskey(d,k), ("deep_jet_fraction", "deep_cyclone_fraction")) &&
            error("component fractions are only valid for deep_profile=matched_components")
    end
    positive("initial", ["jet_radius_m", "jet_width_m", "cyclone_radius_m"])
    initial_profile = get(c["initial"], "profile", "annular")
    initial_profile in ("annular", "rostami2017") || error("initial.profile must be annular or rostami2017")
    if initial_profile == "rostami2017"
        haskey(c["initial"],"core_regularization_m") || error("rostami2017 requires initial.core_regularization_m")
        positive("initial", ["core_regularization_m"])
        for k in ("jet_speed_m_s", "cyclone_speed_m_s")
            c["initial"][k] isa Real && c["initial"][k] >= 0 || error("rostami2017 initial.$k must be nonnegative")
        end
    else
        !haskey(c["initial"],"core_regularization_m") || error("core_regularization_m is only valid for initial.profile=rostami2017")
    end
    c["initial"]["displacement_m"] >= 0 || error("initial displacement must be nonnegative")
    c["initial"]["perturbation"] in ("none", "mode", "broadband") || error("initial.perturbation must be none, mode, or broadband")
    for k in ("mode", "m_max")
        c["initial"][k] isa Integer && 1 <= c["initial"][k] <= 32 || error("initial.$k must be an integer in 1:32")
    end
    mmin = get(c["initial"],"m_min",1)
    mmin isa Integer && 1 <= mmin <= c["initial"]["m_max"] || error("initial.m_min must be an integer in 1:m_max")
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
    if initial_profile == "rostami2017"
        core_cells = c["experiment"]["purpose"] == "smoke" ? 1 : 4
        c["initial"]["core_regularization_m"] >= core_cells*dx ||
            error("rostami2017 core regularization needs at least $core_cells grid cells; convergence still required")
    end
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
