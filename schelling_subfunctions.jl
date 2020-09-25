#subfunctions shared across both binary and trinary cases
function assign_periodic_boundary_conditions(loc,lattice_length)
    #Important function to optimize

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

function count_neighbors(lattice,lattice_length,agent_type,loc)
    I1 = oneunit(loc);
    loc_neighbors = loc-I1:loc+I1;
    num_neighbors = 0; #Negative 1 to account for also detecting the central person as a neighbor
    for loc_neighbor in loc_neighbors
        if lattice[assign_periodic_boundary_conditions(loc_neighbor,lattice_length)]==agent_type && loc_neighbor!=loc
            num_neighbors += 1;
        end
    end
    return num_neighbors
end
