# BCellModel
A combination of Receptor, NFkB, Proliferation, and Apoptosis modules.
We recommend running the model on a server with at least 32 threads. To run the model, one can run the main2.jl file. There are a few options one can set:
* `-v` for the type of multi-thread parallelization, where the option are: `"nonthread"`, `"thread"`, `"spawn"`. `"nonthread"` will run the model linearly and does not parallelize it, while `"thread"` will parallelize the model on static schedule, and `"spawn"` on dynamic schedule.
* `-o` for the destination of output .txt cell lineage file
* `-c` for the destination of output .jld signaling dynamics file (include output of nuclear RelA and cRel)
* `-i` for the destination of output steady state file. In the case it is combined with `-r`, the path will be used to reload from previous steady states, if the parameter distribution & pre-stimulation was already done. This can be used when you would like to rerun a simulation from a .jld file generated previously.
* `-r` for reloading from previous steady states

Example:
```
export JULIA_NUM_THREADS=64 # set number of threads to be used
home_dir="/path/to/dir/BCELL_PROJECT/"
modifier="lineages_125_CD40A_H62" # set the file name for outputs
julia $home_dir'scripts/main2.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' >> $home_dir'job-logs/'$modifier'.out'

```
## Scripts:
- `main2.jl`: function for running the simulations and saving results
- `ConstantParams2.jl`: constant parameters, including stimulus doses, stimulus delay, simulation time, scaling factors, etc.
- `ReactionRates3.jl`: reaction rate parameters for all module
- `ODE_Receptor3.jl`: ODE equations for BCR and CD40 receptor modules
- `ODE_NFkB3.jl`: ODE equations for NFkB module
- `ODE_Apoptosis2.jl`: ODE equations for Apoptosis module
- `ODE_Differentiation.jl`: ODE equations for Differentiation module
- `ODE_Proliferation.jl`: ODE equations for Cell Cycle module
- `SimulateFunctions4+.jl`: pre-simulation and simulation functions
- `HelperFunctions.jl`: helper functions for Michaelis-Menten and Hill functions, as well as parameter distributions
