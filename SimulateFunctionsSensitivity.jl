# reaction rate parameters (added support for fate id as a species) for sensitivity analysis
using Base.Threads;
const SPLOCK = SpinLock(); # spin lock for multi-threading security

# Define structure to hold Cell parameters
#-------------------------------------------------------------------------
mutable struct Cell
    birthday::Float64             # hrs
    current_idx::Int64
    parent_idx::Int64
    generation::Int64
    death_time::Float64                     # (0: surviving, 1: death, 2: division, 3: differentiation)
    div0_time::Float64               # integrator.t at fate (hrs)
    fate_t::Float64           # birthday + fate_t (hrs)
    rates::Matrix{Float64}
    reactionFlux::Matrix{Float64}
    historicFlux::Array{Float64}
    finalConcentration::Array{Float64}
    daughter_1_idx::Int64           # for keeping track of daughter cells when doing depth-first computation
    daughter_2_idx::Int64
end

# Define function to set Cell parameters
#-------------------------------------------------------------------------
function setCell(; birthday = 0.0, current_idx = 0, parent_idx = 0, generation = 1, death_time = 0.0, div0_time = 0.0, fate_t = GLOBAL_END_TIME, rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS), reactionFlux = zeros(Float64, TOTAL_SPECIES, TOTAL_PARAMS), historicFlux = zeros(Float64, TOTAL_SPECIES), finalConcentration = zeros(Float64, TOTAL_SPECIES), daughter_1_idx = 0, daughter_2_idx = 0)
    return Cell(birthday, current_idx, parent_idx, generation, death_time, div0_time, fate_t, rates, reactionFlux, historicFlux, finalConcentration, daughter_1_idx, daughter_2_idx);
end

# Define function to initialize founder cells
#-------------------------------------------------------------------------
function initializeFounderCells!(Srates, allCells)
    @inbounds Threads.@threads for i in 1:FOUNDER_CELL_NUM
        allCells[i] = setCell(current_idx = i, rates = distribute_params(Srates, i));
        # Burn-in period for distributed cells to reach steady-state before simulation
        Pre_simulate!(allCells, i);
        allCells[i].rates[ANTIGEN, INITIALCONC] = ANTIGEN_DOSE;
    end
    nothing
end



# Define function to pre-simulate a founder cell to its steady-state
#-------------------------------------------------------------------------
function Pre_simulate!(allCells, i)
    Bcell = SteadyStateProblem(computeNetworkNettFluxes!, (@view allCells[i].rates[1:end, INITIALCONC]), p=(allCells[i].rates, allCells[i].reactionFlux));
    solution = solve(Bcell, DynamicSS(method, tspan=BURN_IN_PERIOD), abstol=abserr, reltol=relerr, maxiters=1e7);
    allCells[i].rates[1:end, INITIALCONC] = solution.u;
    solution = nothing;
    GC.gc();
    nothing
end

# Define function to simulate a single cell
#-------------------------------------------------------------------------
function Simulate!(allCells, Srates, i, delay)
    @inbounds begin
        has_children = false;
        if CD40L_DELAY == 0
            allCells[i].rates[CD40L, INITIALCONC] = CD40L_DOSE;
        end

        Bcell = DDEProblem(computeNetworkNettFluxes!, allCells[i].rates[1:end, INITIALCONC], delay, (0.0, GLOBAL_END_TIME-allCells[i].birthday), (allCells[i].birthday, allCells[i].rates, allCells[i].reactionFlux, allCells[i].historicFlux); constant_lags=[0.25, 0.75, 1.0, 4.0, 12.0]);
        NFkBsolution = solve(Bcell, MethodOfSteps(method), callback=cellFate, abstol=abserr, reltol=relerr, save_everystep=false, saveat = 0.1, tstops=CD40L_DELAY, maxiters=1e10);
        Bcell2 = ODEProblem(computeNetworkNettFluxes_AP2!, allCells[i].rates[1:end, INITIALCONC], (0.0, GLOBAL_END_TIME-allCells[i].birthday), (allCells[i].birthday, allCells[i].rates, allCells[i].reactionFlux, NFkBsolution));
        solution = solve(Bcell2, method, callback=cellFate, abstol=abserr, reltol=relerr, save_everystep=false, tstops=CD40L_DELAY, maxiters=1e10);

        allCells[i].finalConcentration[1:RECEPTOR_SPECIES+NFKB_SPECIES+DIFF_SPECIES] = NFkBsolution(solution.t[end])[1:RECEPTOR_SPECIES+NFKB_SPECIES+DIFF_SPECIES];
        allCells[i].finalConcentration[RECEPTOR_SPECIES+NFKB_SPECIES+DIFF_SPECIES+1:TOTAL_SPECIES] = solution.u[end][RECEPTOR_SPECIES+NFKB_SPECIES+DIFF_SPECIES+1:TOTAL_SPECIES];

        allCells[i].fate_t = solution.t[end]; # solution.t[end] <= NFkBsolution.t[end] due to definition

        allCells[i].death_time = (solution.u[end][DEATH_TIME] != 0) ? solution.u[end][DEATH_TIME] : 120.0;   # death (for symplified or full apoptosis model)
        allCells[i].div0_time = (solution.u[end][DIV0_TIME] != 0) ? solution.u[end][DIV0_TIME] : 120.0;   # 1st division 
        parameter = RATE_name[(i-1) ÷ ST_num + 1];
        value = (i-1) % ST_num + 1;
        print(string(allCells[i].current_idx), '\t', parameter, '\t', string(value), '\t', string(allCells[i].div0_time), '\t', string(allCells[i].death_time), '\n');

        NFkBsolution = nothing;
        solution = nothing;
        GC.gc();
    end
    # print(has_children, " ", allCells[i].daughter_1_idx, " ", allCells[i].daughter_2_idx)
    return nothing
end

# Define function to simulate a cell lineage (depth-first)
#-------------------------------------------------------------------------
# function Simulate_lineage!(allCells, Srates, i, delay)
#     stack = Int[i]
#     while !isempty(stack)
#         j = pop!(stack)
#         has_children, daughter_1_idx, daughter_2_idx = Simulate!(allCells, Srates, j, delay)
#         if has_children
#             push!(stack, daughter_1_idx)
#             push!(stack, daughter_2_idx)
#         end
#         lock(SPLOCK);

#         cellsSaver = jldopen(cells_fn, "r+");
#         write(cellsSaver, string(j), allCells[j]);
#         close(cellsSaver);

#         print(allCells[j].birthday, '\t', allCells[j].current_idx, '\t', allCells[j].parent_idx, '\t', allCells[j].generation, '\t', allCells[j].fate, '\t', allCells[j].fate_t, '\t', allCells[j].abs_fate_t, '\t', allCells[j].daughter_1_idx, '\t', allCells[j].daughter_2_idx, '\n');
#         unlock(SPLOCK);
#         # GC.gc(true)
#         # ccall(:malloc_trim, Cvoid, (Cint,), 0)
#     end
# end

# Define function to simulate all cell lineages serially
#-------------------------------------------------------------------------
function Simulate_nonthreaded!(allCells, Srates, delay)
    @inbounds for i in 1:FOUNDER_CELL_NUM
        Simulate!(allCells, Srates, i, delay)
    end
end

# Define function to simulate all cell lineages w/ multi-threading (static schedule)
#-------------------------------------------------------------------------
function Simulate_threaded!(allCells, Srates, delay)
    Threads.@threads for i in 1:FOUNDER_CELL_NUM
        Simulate!(allCells, Srates, i, delay)
    end
end

# Define function to simulate all cell lineages w/ multi-threading (dynamic schedule)
#-------------------------------------------------------------------------
function Simulate_spawned!(allCells, Srates, delay)
    tasks = [Threads.@spawn(Simulate!(allCells, Srates, i, delay)) for i in 1:FOUNDER_CELL_NUM]
    out = [fetch(t) for t in tasks]
end
