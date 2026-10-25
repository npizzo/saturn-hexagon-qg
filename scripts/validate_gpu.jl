using SaturnPolarQG, CUDA, LinearAlgebra, Printf

# Run inside a GPU allocation before launching an ensemble. The short case
# includes stochastic events, host/device transfers, diagnostics and the PDE.
CUDA.functional(true) || error("GPU validation requires a working CUDA allocation; no CPU fallback")
CUDA.allowscalar(false)
path = isempty(ARGS) ? joinpath(@__DIR__, "..", "configs", "smoke.toml") : only(ARGS)
c = load_config(path)
c["experiment"]["purpose"] == "smoke" || error("GPU validation expects a short smoke configuration")
a, b = build_simulation(c; device="cpu"), build_simulation(c; device="gpu")
while a.prob.clock.t < c["numerics"]["duration_s"]
    advance!(a); advance!(b)
end
fa, fb = fields(a), fields(b)
errors = Dict(string(k)=>norm(getfield(fa,k)-getfield(fb,k))/max(norm(getfield(fa,k)), eps()) for k in (:s,:psi,:u,:v))
for k in sort(collect(keys(errors)))
    @printf("%s relative L2 difference: %.6e\n", k, errors[k])
end
maximum(values(errors)) < 1e-8 || error("CPU/GPU short-time agreement failed")
a.prob.clock.t == b.prob.clock.t || error("CPU/GPU clocks differ")
a.convection.events == b.convection.events || error("CPU/GPU event counts differ")
a.convection.next_event == b.convection.next_event || error("CPU/GPU event schedules differ")
abs(invariants(b).mean_s) < 1e-15 || error("GPU mean-PV compensation check failed")
println("CPU/GPU smoke agreement passed; this is not a long-time chaotic trajectory comparison.")
