#module schelling_binary_subfunctions

#export initialize_lattice, assign_periodic_boundary_conditions, initialize_bins, update_counts_joint!,step_schelling,count_neighbors

function initialize_lattice_vacancies(lattice_length,num_red_agents,num_blue_agents)
    lattice = zeros(UInt16,lattice_length,lattice_length)
    loc_empty = findall(iszero,lattice)
    loc_blue = sample(loc_empty, num_blue_agents; replace=false, ordered=false)
    lattice[loc_blue] .= 1;
    loc_empty = findall(iszero,lattice)
    loc_red = sample(loc_empty, num_red_agents; replace=false, ordered=false)
    lattice[loc_red] .= 2;
    loc_empty = findall(iszero,lattice)
    locs = [loc_blue;loc_red;loc_empty]
    return lattice,locs
end

function initialize_lattice_three_agents(lattice_length,num_red_agents,num_blue_agents)
    lattice = zeros(UInt16,lattice_length,lattice_length)
    loc_empty = findall(iszero,lattice)
    loc_blue = sample(loc_empty, num_blue_agents; replace=false, ordered=false)
    lattice[loc_blue] .= 1
    loc_empty = findall(iszero,lattice)
    loc_red = sample(loc_empty, num_red_agents; replace=false, ordered=false)
    lattice[loc_red] .= 2
    loc_green = findall(iszero,lattice)
    lattice[loc_green] .= 3
    locs = [loc_blue;loc_red;loc_green]
    return lattice,locs
end

function update_counts_joint!(lattice,bins,bin_IDs,counts_joint)
    #num_blue = Array{UInt16,1}(UndefInitializer(), length(bins))
    #num_red = similar(num_blue)
    for ID in bin_IDs
        num_blue_block = sum(lattice[bins.==ID].==1);
        num_red_block = sum(lattice[bins.==ID].==2);
        counts_joint[num_blue_block+1,num_red_block+1] += 1;
    end
    #num_blue_block = [sum(lattice[bins.==ID].==1) for ID in bin_IDs];
    #num_red_block = [sum(lattice[bins.==ID].==2) for ID in bin_IDs];
    return counts_joint
end

function step_schelling_vacancies(utility_function_blue,
    utility_function_red,lattice,
    lattice_length,
    locs,
    num_blue_total,
    num_red_total,
    num_vacancies)
    #function step_schelling(utility_function_blue,utility_function_red,lattice,r_neighbor,c_neighbor,r_agents,c_agents,r_vacancies,c_vacancies):
    #Randomly attempt a number of moves equal to the total number of agents
    #Could easily optimize here, but keeping things simple for now
    total_agents = num_blue_total+num_red_total;
    agents_idx = sample(1:total_agents, total_agents; replace=true, ordered=true)
    vacancies_idx = sample((total_agents+1):(total_agents+num_vacancies), total_agents; replace=true, ordered=false)

    for i = 1:total_agents
        if agents_idx[i]<=num_blue_total
            agent_type = 1;
        else
            agent_type = 2;
        end
        #agent_type = lattice[r_agents[agents_idx[i]],c_agents[agents_idx[i]]]
        num_neighbors_0 = count_neighbors(lattice,lattice_length,agent_type,locs[agents_idx[i]]);
        num_neighbors_1 = count_neighbors(lattice,lattice_length,agent_type,locs[vacancies_idx[i]]);
        #num_neighbors_0 = count_neighbors(lattice,lattice_length,agent_type,r_agents[agents_idx[i]],c_agents[agents_idx[i]])
        #num_neighbors_1 = count_neighbors(lattice,lattice_length,agent_type,r_vacancies[vacancies_idx[i]],c_vacancies[vacancies_idx[i]])
        #r_agent = r_agents[agents_idx[i]]
        #c_agent = c_agents[agents_idx[i]]
        #r_empty = r_vacancies[vacancies_idx[i]]
        #c_empty = c_vacancies[vacancies_idx[i]]

        if agent_type == 1 #Blue
            change_in_utility = utility_function_blue[num_neighbors_1+1] - utility_function_blue[num_neighbors_0+1];
        else
            change_in_utility = utility_function_red[num_neighbors_1+1] - utility_function_red[num_neighbors_0+1];
        end
        #print(change_in_utility)
        #print(1/(1+np.exp(-change_in_utility)))
        if (1/(1+exp(-change_in_utility)))>rand(1)[1] #Accept move
            #Update lattice
            lattice[locs[agents_idx[i]]] = 0;
            lattice[locs[vacancies_idx[i]]] = agent_type;
            #Update positions of vacancies spots and agents locations
            loc_agent_temp = locs[agents_idx[i]];
            locs[agents_idx[i]] = locs[vacancies_idx[i]];
            locs[vacancies_idx[i]] = loc_agent_temp;
        end
    end
    return lattice,locs
end

function step_schelling_three_agents(utility_function_blue,
    utility_function_red,
    utility_function_green,
    lattice,
    lattice_length,
    locs,
    num_blue_total,
    num_red_total,
    num_green_total)
    #function step_schelling(utility_function_blue,utility_function_red,lattice,r_neighbor,c_neighbor,r_agents,c_agents,r_vacancies,c_vacancies):
    #Randomly attempt a number of moves equal to the total number of agents
    #Could easily optimize here, but keeping things simple for now
    total_agents = num_blue_total+num_red_total+num_green_total;

    agents_blue_idx = sample(1:num_blue_total, num_blue_total; replace=true, ordered=true)
    destinations_blue_idx = sample((num_blue_total+1):total_agents, num_blue_total; replace=true, ordered=false)

    agents_red_idx = sample((num_blue_total+1):(num_blue_total+num_red_total), num_red_total; replace=true, ordered=true)
    destinations_red_idx = sample([(1:num_blue_total);(num_blue_total+num_red_total+1:total_agents)], num_red_total; replace=true, ordered=false)

    agents_green_idx = sample((num_blue_total+num_red_total+1):total_agents, num_green_total; replace=true, ordered=true)
    destinations_green_idx = sample(1:(num_blue_total+num_red_total), num_green_total; replace=true, ordered=false)

    agents_idx = [agents_blue_idx;agents_red_idx;agents_green_idx]
    destinations_idx = [destinations_blue_idx;destinations_red_idx;destinations_green_idx]

    random_idx = shuffle(1:length(agents_idx))
    agents_idx = [agents_idx[i] for i in random_idx]
    destinations_idx = [destinations_idx[i] for i in random_idx]

    for i = 1:total_agents
        if agents_idx[i]<=num_blue_total
            agent_type = 1;
        elseif agents_idx[i]<=num_blue_total+num_red_total
            agent_type = 2;
        else
            agent_type = 3;
        end

        if destinations_idx[i]<=num_blue_total
            destination_agent_type = 1;
        elseif destinations_idx[i]<=num_blue_total+num_red_total
            destination_agent_type = 2;
        else
            destination_agent_type = 3;
        end

        #agent_type = lattice[r_agents[agents_idx[i]],c_agents[agents_idx[i]]]

        num_neighbors_before_1 = count_neighbors(lattice,lattice_length,agent_type,locs[agents_idx[i]]);
        num_neighbors_switched_1 = count_neighbors(lattice,lattice_length,agent_type,locs[destinations_idx[i]]);

        num_neighbors_before_2 = count_neighbors(lattice,lattice_length,destination_agent_type,locs[destinations_idx[i]]);
        num_neighbors_switched_2 = count_neighbors(lattice,lattice_length,destination_agent_type,locs[agents_idx[i]]);

        #num_neighbors_0 = count_neighbors(lattice,lattice_length,agent_type,r_agents[agents_idx[i]],c_agents[agents_idx[i]])
        #num_neighbors_1 = count_neighbors(lattice,lattice_length,agent_type,r_vacancies[vacancies_idx[i]],c_vacancies[vacancies_idx[i]])
        #r_agent = r_agents[agents_idx[i]]
        #c_agent = c_agents[agents_idx[i]]
        #r_empty = r_vacancies[vacancies_idx[i]]
        #c_empty = c_vacancies[vacancies_idx[i]]

        if agent_type == 1 #Blue
            change_in_utility = utility_function_blue[num_neighbors_switched_1+1] - utility_function_blue[num_neighbors_before_1+1];
        elseif agent_type == 2
            change_in_utility = utility_function_red[num_neighbors_switched_1+1] - utility_function_red[num_neighbors_before_1+1];
        else
            change_in_utility = utility_function_green[num_neighbors_switched_1+1] - utility_function_green[num_neighbors_before_1+1];
        end

        if destination_agent_type == 1 #Blue
            change_in_utility = change_in_utility + utility_function_blue[num_neighbors_switched_2+1] - utility_function_blue[num_neighbors_before_2+1];
        elseif destination_agent_type == 2 #Red
            change_in_utility = change_in_utility + utility_function_red[num_neighbors_switched_2+1] - utility_function_red[num_neighbors_before_2+1];
        else
            change_in_utility = change_in_utility + utility_function_green[num_neighbors_switched_2+1] - utility_function_green[num_neighbors_before_2+1];
        end
        #print(change_in_utility)
        #print(1/(1+np.exp(-change_in_utility)))
        if (1/(1+exp(-change_in_utility)))>rand(1)[1] #Accept move
            #Update lattice
            lattice[locs[agents_idx[i]]] = destination_agent_type;
            lattice[locs[destinations_idx[i]]] = agent_type;
            #Update positions of vacancies spots and agents locations
            loc_agent_temp = locs[agents_idx[i]];
            locs[agents_idx[i]] = locs[destinations_idx[i]];
            locs[destinations_idx[i]] = loc_agent_temp;
        end
    end
    return lattice,locs
end

function run_schelling_sim_binary_vacancies(;lattice_length = 60,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 500,
    utility_function_blue = [0, 1, 2, 3, 4, 5, 6, 7, 8] .+ 1,
    utility_function_red = [0, 1, 2, 3, 4, 5, 6, 7, 8],
    bin_length = 5
    )
    num_burnin_steps = 100
    num_snapshots = 10
    num_steps_per_snapshot = Int(num_simulation_steps/num_snapshots)
    num_red_total = Int(round(frac_red_agents * (lattice_length^2)))
    num_blue_total = Int(round(frac_blue_agents * (lattice_length^2)))
    lattice, locs =
        initialize_lattice_vacancies(lattice_length, num_red_total, num_blue_total)
    lattice_snapshots = zeros(Int8,size(lattice,1),size(lattice,2),num_snapshots)
    num_vacancies = length(locs) - num_red_total - num_blue_total
    bins, bin_IDs = initialize_bins(lattice, lattice_length, bin_length)
    counts_joint = zeros(bin_length^2 + 1, bin_length^2 + 1)
    for step = 1:num_burnin_steps
        lattice, locs = step_schelling_vacancies(
            utility_function_blue,
            utility_function_red,
            lattice,
            lattice_length,
            locs,
            num_blue_total,
            num_red_total,
            num_vacancies,
        )
    end
    for step = 1:num_simulation_steps
        lattice, locs = step_schelling_vacancies(
            utility_function_blue,
            utility_function_red,
            lattice,
            lattice_length,
            locs,
            num_blue_total,
            num_red_total,
            num_vacancies,
        )
        update_counts_joint!(lattice, bins, bin_IDs, counts_joint)
        if mod(step,num_steps_per_snapshot)==0
            lattice_snapshots[:,:,Int(step/num_steps_per_snapshot)] = lattice
        end
    end
    return counts_joint, lattice_snapshots,utility_function_blue,utility_function_red
end


function run_schelling_sim_trinary(;lattice_length = 60,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 500,
    utility_function_blue = [0, 1, 2, 3, 4, 5, 6, 7, 8] .+ 1,
    utility_function_red = [0, 1, 2, 3, 4, 5, 6, 7, 8],
    utility_function_green = [0, 1, 2, 3, 4, 5, 6, 7, 8] .- 1,
    bin_length = 5
    )
    num_burnin_steps = 100
    num_snapshots = 10
    num_steps_per_snapshot = Int(num_simulation_steps/num_snapshots)
    num_red_total = Int(round(frac_red_agents * (lattice_length^2)))
    num_blue_total = Int(round(frac_blue_agents * (lattice_length^2)))
    lattice, locs =
        initialize_lattice_three_agents(lattice_length, num_red_total, num_blue_total)
    lattice_snapshots = zeros(Int8,size(lattice,1),size(lattice,2),num_snapshots)
    num_green_total = length(locs) - num_red_total - num_blue_total
    bins, bin_IDs = initialize_bins(lattice, lattice_length, bin_length)
    counts_joint = zeros(bin_length^2 + 1, bin_length^2 + 1)
    for step = 1:num_burnin_steps
        lattice, locs = step_schelling_three_agents(
            utility_function_blue,
            utility_function_red,
            utility_function_green,
            lattice,
            lattice_length,
            locs,
            num_blue_total,
            num_red_total,
            num_green_total,
        )
    end
    for step = 1:num_simulation_steps
        lattice, locs = step_schelling_three_agents(
            utility_function_blue,
            utility_function_red,
            utility_function_green,
            lattice,
            lattice_length,
            locs,
            num_blue_total,
            num_red_total,
            num_green_total,
        )
        update_counts_joint!(lattice, bins, bin_IDs, counts_joint)
        if mod(step,num_steps_per_snapshot)==0
            lattice_snapshots[:,:,Int(step/num_steps_per_snapshot)] = lattice
        end
    end
    return counts_joint,lattice_snapshots,utility_function_blue,utility_function_red,utility_function_green
end
