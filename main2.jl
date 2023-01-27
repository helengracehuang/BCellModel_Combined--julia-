# Include libraries / packages required
#--------------------------------------------
using DifferentialEquations;

using ArgParse; # Argument parsing

# Benchmarking & profiling packages
using BenchmarkTools;
using StatProfilerHTML;

# Visualization packages
using DataFrames;
using Plots;
using Gadfly;
using Cairo;
using JLD;

# Include source files
#--------------------------------------------
include("ReactionRates2.jl");
include("HelperFunctions.jl");

include("ODE_Receptor3.jl");
include("ODE_NFkB3.jl");
include("ODE_Apoptosis2.jl");
include("ODE_Differentiation.jl");
include("ODE_Proliferation.jl");

include("SimulateFunctions4.jl");

# Argument parsing w/ command line input
#--------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--version", "-v"
            help = "nonthread, thread (static schedule), or spawn (dynamic schedule) version of lineage simulation"
            arg_type = String
            default = "nonthread"
            required = true
        "--initial", "-i"
            help = "reaction rates and steady states for each starting cell after initialization (.jld format)"
            arg_type = String
            default = "initial.jld"
            required = true
        "--cells", "-c"
            help = "simulated outputs from each cell (.jld format)"
            arg_type = String
            default = "cells.jld"
            required = true
        "--output", "-o"
            help = "output file name for simulated cell lineages"
            arg_type = String
            default = "output.txt"
            required = true
        "--reload", "-r"
            help = "whether to reload from previous steady states, provided by the --initial command"
            arg_type = Bool
            action = :store_true
            default = false
            required = false
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
const version = get(parsed_args, "version", "nonthread");
const steady_fn = get(parsed_args, "initial", "initial.jld");
const cells_fn = get(parsed_args, "cells", "cells.jld");
const output_fn = get(parsed_args, "output", "output.txt");
const reload = get(parsed_args, "reload", false);

# Define standard reaction rates
#--------------------------------------------
rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS);
setAllRates!(rates);
const Srates = rates;

# Define ODEs for network
#--------------------------------------------
# Attempt to speed up the solver by splitting into 2 parts
# time-independent ODE (pre-simulation: phase = 1, only receptor & NFkB needs steady state simulation)
function computeNetworkNettFluxes!(nettFlux, concentration, (Srates, reactionFlux), time)
    computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    nothing
end

function computeNetworkNettFluxes_AP1!(nettFlux, concentration, (Srates, reactionFlux), time)
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeNetworkNettFluxes!(nettFlux, concentration, delay, (birthday, Srates, reactionFlux, historicFlux), time)
    computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, time);
    computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, delay, historicFlux, time);
    # computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    computeDiffNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    # computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end

function computeNetworkNettFluxes_AP2!(nettFlux, concentration, (birthday, Srates, reactionFlux, inputCurves), time)
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, time, birthday, inputCurves);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end


# Callback function to detect cell death, mitotic and differentiation events
#----------------------------------------------------------------------
function condition(out, u, t, integrator)
    out[1] = u[CPARP] - 2500;
    out[2] = u[CDH1] - 0.2;
    out[3] = u[IRF4] - 0.65*u[BCL6] - 1.2;
    out[4] = t - CD40L_DELAY;
end

function affect!(integrator, index)
    if (index == 1)
        integrator.u[TOTAL_SPECIES] = 1;
        terminate!(integrator);
    elseif (index == 2)
        if (integrator.u[CYCB] > 2.0)
            integrator.u[MASS] /= 2.0;
            integrator.u[GEN] += 1.0;
            integrator.u[TOTAL_SPECIES] = 2;
            terminate!(integrator);
        end
    elseif (index == 3)
        integrator.u[TOTAL_SPECIES] = 3;
        terminate!(integrator);
    elseif (index == 4)
        integrator.u[CD40L] = CD40L_DOSE;
        integrator.u[ANTIGEN] = 0;
        # integrator.u[ABCR] = 0;
    end
    nothing
end

cellFate = VectorContinuousCallback(condition, affect!, nothing, 4, save_positions=(true, true));

# Set up a structure to hold cells
#-------------------------------------------------------------------------
allCells = Vector{Cell}(undef, FOUNDER_CELL_NUM);
if reload
    allCells = load(steady_fn, "allCells");
    print("Old steady states loaded!")
else
    initializeFounderCells!(Srates, allCells);
    JLD.save(steady_fn, "allCells", allCells);
end

# Define default delay functions (t < tau)
#--------------------------------------------
# const delay(historicFlux, p, t) = (historicFlux .= 0.0);
const delay(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : zeros(TOTAL_SPECIES);

# Define input NEMO curve (according to stimulus type)
#--------------------------------------------
# const NEMOCurve = specifyNEMOCurve(NEMO_shape = NEMO_TYPE);
# const NIKCurve = specifyNEMOCurve(NIK_shape = NIK_TYPE);

print("birthday", '\t', "current_idx", '\t', "parent_idx", '\t', "generation", '\t', "fate", '\t', "fate_t", '\t', "abs_fate_t", '\t', "daughter_1_idx", '\t', "daughter_2_idx", '\n');

# Simulate cell lineages
#--------------------------------------------
if version == "nonthread"
    @time Simulate_nonthreaded!(allCells, Srates, delay);
elseif version == "thread"
    @time Simulate_threaded!(allCells, Srates, delay);
elseif version == "spawn"
    @time Simulate_spawned!(allCells, Srates, delay);
else
    print("Please input the correct version -v: nonthread, thread, or spawn.");
end

# @time Simulate_nonthreaded!(allCells, Srates, delay);
# @time Simulate_threaded!(allCells, Srates, delay);
# @time Simulate_spawned!(allCells, Srates, delay);
JLD.save(cells_fn, "allCells", allCells);

# Output information about all cells (for visualization)
#--------------------------------------------
out = open(output_fn, "w");
write(out, "birthday", '\t', "current_idx", '\t', "parent_idx", '\t', "generation", '\t', "fate", '\t', "fate_t", '\t', "abs_fate_t", '\t', "daughter_1_idx", '\t', "daughter_2_idx", '\n');
for i in 1:length(allCells)
    write(out, string(allCells[i].birthday), '\t', string(allCells[i].current_idx), '\t', string(allCells[i].parent_idx), '\t', string(allCells[i].generation), '\t', string(allCells[i].fate), '\t', string(allCells[i].fate_t), '\t', string(allCells[i].abs_fate_t), '\t', string(allCells[i].daughter_1_idx), '\t', string(allCells[i].daughter_2_idx), '\n');
end
close(out);
