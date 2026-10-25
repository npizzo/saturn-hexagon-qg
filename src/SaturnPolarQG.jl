module SaturnPolarQG
using GeophysicalFlows, FFTW, CUDA, LinearAlgebra, Statistics, Random
using TOML, JLD2, SHA, Dates, Printf
const QG = GeophysicalFlows.SingleLayerQG
include("Convection.jl")
include("PolarDiagnostics.jl")
using .Convection
using .PolarDiagnostics
include("Configuration.jl")
include("Dynamics.jl")
include("RunIO.jl")
export load_config, validate_config, build_simulation, advance!, fields, invariants,
       save_checkpoint, restore_checkpoint!, run_experiment
end
