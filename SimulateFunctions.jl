using Base.Threads;
const SPLOCK = SpinLock(); # spin lock for multi-threading security

# Define structure to hold Cell parameters
#-------------------------------------------------------------------------
mutable struct Cell
    birthday::Float64             # hrs
    current_idx::Int64
    parent_idx::Int64
    generation::Int64
    fate::Int64                     # (0: surviving, 1: death, 2: division, 3: differentiation)
    fate_t::Float64               # integrator.t at fate (hrs)
    abs_fate_t::Float64           # birthday + fate_t (hrs)
    rates::Matrix{Float64}
    reactionFlux::Matrix{Float64}
    historicFlux::Matrix{Float64}
    finalConcentration::Array{Float64}
    daughter_1_idx::Int64           # for keeping track of daughter cells when doing depth-first computation
    daughter_2_idx::Int64
end

# Define function to set Cell parameters
#-------------------------------------------------------------------------
function setCell(; birthday = 0.0, current_idx = 0, parent_idx = 0, generation = 1, fate = 0, fate_t = 0.0, abs_fate_t = GLOBAL_END_TIME, rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS), reactionFlux = zeros(Float64, TOTAL_SPECIES, TOTAL_PARAMS), historicFlux = zeros(Float64, TOTAL_SPECIES, TOTAL_PARAMS), finalConcentration = zeros(Float64, TOTAL_SPECIES), daughter_1_idx = 0, daughter_2_idx = 0)
    return Cell(birthday, current_idx, parent_idx, generation, fate, fate_t, abs_fate_t, rates, reactionFlux, historicFlux, finalConcentration, daughter_1_idx, daughter_2_idx);
end

# Define function to initialize founder cells
#-------------------------------------------------------------------------
function initializeFounderCells!(Srates, allCells)
    @inbounds Threads.@threads for i in 1:FOUNDER_CELL_NUM
        allCells[i] = setCell(current_idx = i, rates = distribute_params(Srates));
        # Burn-in period for distributed cells to reach steady-state before simulation
        Pre_simulate!(allCells, i);
    end
    nothing
end



# Define function to pre-simulate a founder cell to its steady-state
#-------------------------------------------------------------------------
function Pre_simulate!(allCells, i)
    Bcell = SteadyStateProblem(computeNetworkNettFluxes!, (@view allCells[i].rates[1:end, INITIALCONC]), p=(allCells[i].rates, allCells[i].reactionFlux));
    solution = solve(Bcell, DynamicSS(method), tspan=BURN_IN_PERIOD, abstol=abserr, reltol=relerr, maxiters=1e7);
    allCells[i].rates[1:end, INITIALCONC] .= solution.u
    nothing
end

# Define function to simulate a single cell
#-------------------------------------------------------------------------
function Simulate!(allCells, Srates, i, IKKCurve, delay)
    @inbounds begin
        has_children = false;
        Bcell = DDEProblem(computeNetworkNettFluxes!, (@view allCells[i].rates[1:end, INITIALCONC]), delay, (0.0, GLOBAL_END_TIME-allCells[i].birthday), (allCells[i].birthday, allCells[i].rates, allCells[i].reactionFlux, allCells[i].historicFlux); constant_lags=[0.25, 0.75, 1.0, 4.0, 12.0]);
        solution = solve(Bcell, MethodOfSteps(method), callback=cellFate, abstol=abserr, reltol=relerr, save_everystep=false, maxiters=1e10);

        allCells[i].finalConcentration .= solution.u[end];
        allCells[i].fate_t = solution.t[end];
        if solution.u[end][TPARP] ??? 25000.0   # death
            allCells[i].fate = 1;
            allCells[i].abs_fate_t = allCells[i].birthday + allCells[i].fate_t;
        elseif solution.u[end][CDH1] ??? 0.2   # division
            allCells[i].fate = 2;
            allCells[i].abs_fate_t = allCells[i].birthday + allCells[i].fate_t;
            if allCells[i].generation < MAX_GEN   # division has not reached max gen
                # generate daughter cells if the cell divides within (GLOBAL_END_TIME - birthday)
                # calculate division factor for the 2 daughter cells
                division_factor = max(0.01, rand(Normal(1, DAUGHTER_PARTITION_CV)));
                # push daughter 1
                lock(SPLOCK) # lock the array so that it does not cause confusion in multi-threading
                allCells[i].daughter_1_idx = length(allCells)+1
                push!(allCells, setCell(birthday = allCells[i].abs_fate_t, current_idx = allCells[i].daughter_1_idx, parent_idx = allCells[i].current_idx, generation = allCells[i].generation+1, rates = distribute_daughter_params(allCells[i].finalConcentration, allCells[i].rates, division_factor, 1)));
                # push daughter 2
                allCells[i].daughter_2_idx = length(allCells)+1
                push!(allCells, setCell(birthday = allCells[i].abs_fate_t, current_idx = allCells[i].daughter_2_idx, parent_idx = allCells[i].current_idx, generation = allCells[i].generation+1, rates = distribute_daughter_params(allCells[i].finalConcentration, allCells[i].rates, division_factor, 2)));
                unlock(SPLOCK);
                has_children = true;
            end
        elseif solution.u[end][IRF4] - 0.65*solution.u[end][BCL6] ??? 1.2   # differentiation
            allCells[i].fate = 3;
            allCells[i].abs_fate_t = allCells[i].birthday + allCells[i].fate_t;
        else # kept surviving till the end
            allCells[i].fate_t = GLOBAL_END_TIME - allCells[i].birthday;
        end
    end
    print(has_children, " ", allCells[i].daughter_1_idx, " ", allCells[i].daughter_2_idx)
    return has_children, allCells[i].daughter_1_idx, allCells[i].daughter_2_idx
end

# Define function to simulate a cell lineage (depth-first)
#-------------------------------------------------------------------------
function Simulate_lineage!(allCells, Srates, i, IKKCurve, delay)
    stack = Int[i]
    while !isempty(stack)
        j = pop!(stack)
        has_children, daughter_1_idx, daughter_2_idx = Simulate!(allCells, Srates, j, IKKCurve, delay)
        if has_children
            push!(stack, daughter_1_idx)
            push!(stack, daughter_2_idx)
        end
    end
end

# Define function to simulate all cell lineages serially
#-------------------------------------------------------------------------
function Simulate_nonthreaded!(allCells, Srates, IKKCurve, delay)
    @inbounds for i in 1:FOUNDER_CELL_NUM
        Simulate_lineage!(allCells, Srates, i, IKKCurve, delay)
    end
end

# Define function to simulate all cell lineages w/ multi-threading (static schedule)
#-------------------------------------------------------------------------
function Simulate_threaded!(allCells, Srates, IKKCurve, delay)
    Threads.@threads for i in 1:FOUNDER_CELL_NUM
        Simulate_lineage!(allCells, Srates, i, IKKCurve, delay)
    end
end

# Define function to simulate all cell lineages w/ multi-threading (dynamic schedule)
#-------------------------------------------------------------------------
function Simulate_spawned!(allCells, Srates, IKKCurve, delay)
    tasks = [Threads.@spawn(Simulate_lineage!(allCells, Srates, i, IKKCurve, delay)) for i in 1:FOUNDER_CELL_NUM]
    out = [fetch(t) for t in tasks]
end
