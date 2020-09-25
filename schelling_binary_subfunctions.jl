using Random
#module schelling_binary_subfunctions

#export initialize_lattice, assign_periodic_boundary_conditions, initialize_bins, update_counts_joint!,step_schelling,count_neighbors

function initialize_lattice_binary(lattice_length,num_blue_agents)
    lattice = zeros(UInt16,lattice_length,lattice_length)
    loc_empty = findall(iszero,lattice)
    loc_blue = sample(loc_empty, num_blue_agents; replace=false, ordered=false)
    lattice[loc_blue] .= 1;
    loc_red = findall(iszero,lattice)
    lattice[loc_empty] .= 2;
    locs = [loc_blue;loc_red]
    return lattice,locs
end

function update_counts_single!(lattice,bins,bin_IDs,counts_single)
    #num_blue = Array{UInt16,1}(UndefInitializer(), length(bins))
    #num_red = similar(num_blue)
    for ID in bin_IDs
        num_blue_block = sum(lattice[bins.==ID].==1);
        # num_red_block = sum(lattice[bins.==ID].==2);
        counts_single[num_blue_block+1] += 1;
    end
    #num_blue_block = [sum(lattice[bins.==ID].==1) for ID in bin_IDs];
    #num_red_block = [sum(lattice[bins.==ID].==2) for ID in bin_IDs];
    return counts_single
end

function step_schelling_binary(utility_function_blue,
    utility_function_red,
    lattice,
    lattice_length,
    locs,
    num_blue_total,
    num_red_total)
    #function step_schelling(utility_function_blue,utility_function_red,lattice,r_neighbor,c_neighbor,r_agents,c_agents,r_empties,c_empties):
    #Randomly attempt a number of moves equal to the total number of agents
    #Could easily optimize here, but keeping things simple for now
    total_agents = num_blue_total+num_red_total;
    agents_idx = [sample(1:num_blue_total, num_blue_total; replace=true, ordered=true);
        sample((num_blue_total+1):total_agents, num_red_total; replace=true, ordered=true)]
    destination_idx = [sample((num_blue_total+1):total_agents, num_blue_total; replace=true, ordered=true);
        sample(1:num_blue_total, num_red_total; replace=true, ordered=true)]
    random_idx = shuffle(1:length(agents_idx))
    agents_idx = [agents_idx[i] for i in random_idx]
    destination_idx = [destination_idx[i] for i in random_idx]
    for i = 1:total_agents
        if agents_idx[i]<=num_blue_total
            agent_type = 1;
            destination_agent_type = 2;
        else
            agent_type = 2;
            destination_agent_type = 1;
        end
        #agent_type = lattice[r_agents[agents_idx[i]],c_agents[agents_idx[i]]]
        num_neighbors_before_1 = count_neighbors(lattice,lattice_length,agent_type,locs[agents_idx[i]]);
        num_neighbors_switched_1 = count_neighbors(lattice,lattice_length,agent_type,locs[destination_idx[i]]);

        num_neighbors_before_2 = count_neighbors(lattice,lattice_length,destination_agent_type,locs[destination_idx[i]]);
        num_neighbors_switched_2 = count_neighbors(lattice,lattice_length,destination_agent_type,locs[agents_idx[i]]);

        #num_neighbors_0 = count_neighbors(lattice,lattice_length,agent_type,r_agents[agents_idx[i]],c_agents[agents_idx[i]])
        #num_neighbors_1 = count_neighbors(lattice,lattice_length,agent_type,r_empties[empties_idx[i]],c_empties[empties_idx[i]])
        #r_agent = r_agents[agents_idx[i]]
        #c_agent = c_agents[agents_idx[i]]
        #r_empty = r_empties[empties_idx[i]]
        #c_empty = c_empties[empties_idx[i]]

        if agent_type == 1 #Blue
            change_in_utility = utility_function_blue[num_neighbors_switched_1+1] - utility_function_blue[num_neighbors_before_1+1];
            change_in_utility = change_in_utility + utility_function_red[num_neighbors_switched_2+1] - utility_function_red[num_neighbors_before_2+1];
        else
            change_in_utility = utility_function_red[num_neighbors_switched_1+1] - utility_function_red[num_neighbors_before_1+1];
            change_in_utility = change_in_utility + utility_function_blue[num_neighbors_switched_2+1] - utility_function_blue[num_neighbors_before_2+1];
        end
        #print(change_in_utility)
        #print(1/(1+np.exp(-change_in_utility)))
        if (1/(1+exp(-change_in_utility)))>rand(1)[1] #Accept move
            #Update lattice
            lattice[locs[agents_idx[i]]] = destination_agent_type;
            lattice[locs[destination_idx[i]]] = agent_type;
            #Update positions of empties spots and agents locations
            loc_agent_temp = locs[agents_idx[i]];
            locs[agents_idx[i]] = locs[destination_idx[i]];
            locs[destination_idx[i]] = loc_agent_temp;
        end
    end
    return lattice,locs
end

function run_schelling_sim_binary(;lattice_length = 60,
    frac_blue_agents = 0.5,
    num_simulation_steps = 500,
    utility_function_blue = [0, 1, 2, 3, 4, 5, 6, 7, 8] .+ 1,
    utility_function_red = [0, 1, 2, 3, 4, 5, 6, 7, 8],
    bin_length = 5
    )
    num_burnin_steps = 100
    num_snapshots = 10
    num_steps_per_snapshot = Int(num_simulation_steps/num_snapshots)

    num_blue_total = Int(round(frac_blue_agents * (lattice_length^2)))
    num_red_total = lattice_length^2 - num_blue_total
    lattice, locs =
        initialize_lattice_binary(lattice_length, num_blue_total)
    lattice_snapshots = zeros(Int8,size(lattice,1),size(lattice,2),num_snapshots)
    bins, bin_IDs = initialize_bins(lattice, lattice_length, bin_length)
    counts_single = zeros(bin_length^2 + 1,1)
    for step = 1:num_burnin_steps
        lattice, locs = step_schelling_binary(
            utility_function_blue,
            utility_function_red,
            lattice,
            lattice_length,
            locs,
            num_blue_total,
            num_red_total,
        )
    end
    for step = 1:num_simulation_steps
        lattice, locs = step_schelling_binary(
            utility_function_blue,
            utility_function_red,
            lattice,
            lattice_length,
            locs,
            num_blue_total,
            num_red_total,
        )
        update_counts_single!(lattice, bins, bin_IDs, counts_single)
        if mod(step,num_steps_per_snapshot)==0
            lattice_snapshots[:,:,Int(step/num_steps_per_snapshot)] = lattice
        end
    end
    return counts_single, lattice_snapshots,utility_function_blue,utility_function_red
end
