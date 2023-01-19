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

# Include source files
#--------------------------------------------
include("HelperFunctions.jl");

include("ODE_NFkB2.jl");
include("ODE_Apoptosis.jl");
include("ODE_Differentiation.jl");
include("ODE_Proliferation2.jl");

include("SimulateFunctions.jl");

# Argument parsing w/ command line input
#--------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--version", "-v"
            help = "nonthread, thread (static schedule), or spawn (dynamic schedule) version of lineage simulation"
            arg_type = String
            default = "nonthread"
        "--output", "-o"
            help = "output file name for simulated cell lineages"
            arg_type = String
            default = "output.txt"
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
const version = get(parsed_args, "version", "nonthread");
const output_fn = get(parsed_args, "output", "output.txt");

# Define standard reaction rates
#--------------------------------------------
rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS);
setAllRates!(rates);
const Srates = rates;

# Define ODEs for network
#--------------------------------------------
# time-independent ODE (pre-simulation: phase = 1)
function computeNetworkNettFluxes!(nettFlux, concentration, (Srates, reactionFlux), time)
    computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeDiffNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeNetworkNettFluxes!(nettFlux, concentration, delay, (birthday, IKKCurve, Srates, reactionFlux, historicFlux), time)
    computeNFkBNettFluxes!(nettFlux, concentration, delay, reactionFlux, Srates, 2, time, birthday, IKKCurve, historicFlux);
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    computeDiffNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end

# Callback function to detect cell death, mitotic and differentiation events
#----------------------------------------------------------------------
function condition(out, u, t, integrator)
    out[1] = u[TPARP]/u[PARP] - 1;
    out[2] = u[CDH1] - 0.2;
    out[3] = u[IRF4] - 0.65*u[BCL6] - 1.2;
end

function affect!(integrator, index)
    if (index == 1)
        terminate!(integrator);
    elseif (index == 2)
        if (integrator.u[CYCB] > 2.0)
            integrator.u[MASS] /= 2.0;
            integrator.u[GEN] += 1.0;
            terminate!(integrator);
        end
    elseif (index == 3)
        terminate!(integrator);
    end
    nothing
end

cellFate = VectorContinuousCallback(condition, affect!, nothing, 3, save_positions=(false, false));

# Set up a structure to hold cells
#-------------------------------------------------------------------------
allCells = Vector{Cell}(undef, FOUNDER_CELL_NUM);
initializeFounderCells!(Srates, allCells);

# Define default delay functions (t < tau)
#--------------------------------------------
const delay(historicFlux, p, t) = (historicFlux .= 0.0);

# Define input IKK curve (according to stimulus type)
#--------------------------------------------
const IKKCurve = specifyIKKCurve(IKK_shape = IKK_TYPE);

# Simulate cell lineages
#--------------------------------------------
# if version == "nonthread"
#     @time Simulate_nonthreaded!(allCells, Srates, IKKCurve, delay);
# elseif version == "thread"
#     @time Simulate_threaded!(allCells, Srates, IKKCurve, delay);
# elseif version == "spawn"
#     @time Simulate_spawned!(allCells, Srates, IKKCurve, delay);
# else
#     print("Please input the correct version -v: nonthread, thread, or spawn.");
# end

@time Simulate_nonthreaded!(allCells, Srates, IKKCurve, delay);
@time Simulate_threaded!(allCells, Srates, IKKCurve, delay);
@time Simulate_spawned!(allCells, Srates, IKKCurve, delay);

# Output information about all cells (for visualization)
#--------------------------------------------
out = open(output_fn, "w");
write(out, "birthday", '\t', "current_idx", '\t', "parent_idx", '\t', "generation", '\t', "fate", '\t', "fate_t", '\t', "abs_fate_t", '\t', "daughter_1_idx", '\t', "daughter_2_idx", '\n');
for i in 1:length(allCells)
    write(out, string(allCells[i].birthday), '\t', string(allCells[i].current_idx), '\t', string(allCells[i].parent_idx), '\t', string(allCells[i].generation), '\t', string(allCells[i].fate), '\t', string(allCells[i].fate_t), '\t', string(allCells[i].abs_fate_t), '\t', string(allCells[i].daughter_1_idx), '\t', string(allCells[i].daughter_2_idx), '\n');
end
close(out);
