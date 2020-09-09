#module schelling_binary_subfunctions

#export initialize_lattice, assign_periodic_boundary_conditions, initialize_bins, update_counts_joint!,step_schelling,count_neighbors




function initialize_lattice(lattice_length,num_red_agents,num_blue_agents)
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

function assign_periodic_boundary_conditions(loc,lattice_length)
    #Takes in an array of locations and shifts coordinates that are less than 1 or greater than lattice_length by the appropriate amount to account for periodic boundary conditions
    #Test to see if it is better to use the if statements of just modulus all the elements
    if loc[1] < 1
        r = loc[1]+lattice_length
    elseif loc[1] > lattice_length
        r = loc[1]-lattice_length
    else
        r = loc[1]
    end

    if loc[2] < 1
        c= loc[2]+lattice_length
    elseif loc[2] > lattice_length
        c = loc[2]-lattice_length
    else
        c = loc[2]
    end

    return CartesianIndex((r,c))
    #     for i in eachindex[locs]
#         r,c = loc[1],loc[2];
#         if loc[1]==0||loc[1]==lattice_length
#             r = mod(loc[1],lattice_length);
#         end
#         if loc[2]==0||loc[2]==lattice_length
#             c = mod(loc[1],lattice_length);
#         end
#         locs[]
#     end
end

function initialize_bins(lattice,lattice_length,bin_length)
    bins = similar(lattice);
    idx =  CartesianIndices(bins);
    shift_bin_length = CartesianIndex(bin_length-1,bin_length-1);
    bin_counter = 1;
    for r=1:bin_length:lattice_length
        for c = 1:bin_length:lattice_length
            top_left_idx = CartesianIndex(r,c);
            bins[top_left_idx:(top_left_idx+shift_bin_length)].=bin_counter;
            bin_counter += 1;
        end
    end
    bin_IDs = unique(bins)
    return bins, bin_IDs
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

function count_neighbors(lattice,lattice_length,agent_type,loc)
    I1 = oneunit(loc);
    loc_neighbors = loc-I1:loc+I1;
    #deleteat!(loc_neighbors,5);
    #loc_neighbors = [assign_periodic_boundary_conditions(loc_neighbor,lattice_length) for loc_neighbor in loc_neighbors];
    num_neighbors = 0; #Negative 1 to account for also detecting the central person as a neighbor
    for loc_neighbor in loc_neighbors
        if lattice[assign_periodic_boundary_conditions(loc_neighbor,lattice_length)]==agent_type && loc_neighbor!=loc
            num_neighbors += 1;
        end
    end

    #neighborhood_0 = lattice[r_agents[agents_idx[i]]-1:r_agents[agents_idx[i]]+1,
    #                         c_agents[agents_idx[i]]-1:c_agents[agents_idx[i]]+1]
    #num_neighbors = np.sum(lattice[r_neighbors,c_neighbors] == agent_type)
    return num_neighbors

    #num_neighbors_0 = np.sum(lattice[r_neighbor[r_agents[agents_idx[i]],c_agents[agents_idx[i]],:],
    #                                 c_neighbor[r_agents[agents_idx[i]],c_agents[agents_idx[i]],:]] == agent_type)
    #num_neighbors_1 = np.sum(lattice[r_neighbor[r_empties[empties_idx[i]],c_empties[empties_idx[i]],:],
    #                                 c_neighbor[r_empties[empties_idx[i]],c_empties[empties_idx[i]],:]] == agent_type)
    #neighbors_0 = np.sum(lattice[r_agents[agents_idx[i]] + [-1, 1, 0, 0, -1, -1, 1, 1],
    #                             c_agents[agents_idx[i]] + [0, 0, -1, 1, -1, 1, -1, 1]] == agent_type)
    #neighbors_1 = np.sum(lattice[r_empties[empties_idx[i]] + [-1, 1, 0, 0, -1, -1, 1, 1],
    #                             c_empties[empties_idx[i]] + [0, 0, -1, 1, -1, 1, -1, 1]] == agent_type)
end

function step_schelling(utility_function_blue,utility_function_red,lattice,lattice_length,locs,num_blue_total,num_red_total,num_empties)
#function step_schelling(utility_function_blue,utility_function_red,lattice,r_neighbor,c_neighbor,r_agents,c_agents,r_empties,c_empties):
    #Randomly attempt a number of moves equal to the total number of agents
    #Could easily optimize here, but keeping things simple for now
    total_agents = num_blue_total+num_red_total;
    agents_idx = sample(1:total_agents, total_agents; replace=true, ordered=true)
    empties_idx = sample((total_agents+1):(total_agents+num_empties), total_agents; replace=true, ordered=false)

    for i = 1:total_agents
        if agents_idx[i]<=num_blue_total
            agent_type = 1;
        else
            agent_type = 2;
        end
        #agent_type = lattice[r_agents[agents_idx[i]],c_agents[agents_idx[i]]]
        num_neighbors_0 = count_neighbors(lattice,lattice_length,agent_type,locs[agents_idx[i]]);
        num_neighbors_1 = count_neighbors(lattice,lattice_length,agent_type,locs[empties_idx[i]]);
        #num_neighbors_0 = count_neighbors(lattice,lattice_length,agent_type,r_agents[agents_idx[i]],c_agents[agents_idx[i]])
        #num_neighbors_1 = count_neighbors(lattice,lattice_length,agent_type,r_empties[empties_idx[i]],c_empties[empties_idx[i]])
        #r_agent = r_agents[agents_idx[i]]
        #c_agent = c_agents[agents_idx[i]]
        #r_empty = r_empties[empties_idx[i]]
        #c_empty = c_empties[empties_idx[i]]

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
            lattice[locs[empties_idx[i]]] = agent_type;
            #Update positions of empties spots and agents locations
            loc_agent_temp = locs[agents_idx[i]];
            locs[agents_idx[i]] = locs[empties_idx[i]];
            locs[empties_idx[i]] = loc_agent_temp;
        end
    end
    return lattice,locs
end




#end
