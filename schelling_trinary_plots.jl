using Plots

gr()

function utility_function_plot(utility_function_blue,utility_function_red;utility_function_green = [],legend_choice = :best)
    # trace_blue = scatter(;x=0:8, y=utility_function_blue, mode="lines+markers")
    # trace_red = scatter(;x=0:8, y=utility_function_red, mode="lines+markers")
    if isempty(utility_function_green)
        plot([0:8],[utility_function_blue,utility_function_red],
            label = ["Blue agent" "Red agent"],
            linewidth = 2,
            linecolor = [:blue :red],
            legend = legend_choice,
        )
    else
        plot([0:8],[utility_function_blue,utility_function_red,utility_function_green],
            label = ["Blue agent" "Red agent" "Green agent"],
            linewidth = 2,
            linecolor = [:blue :red :green],
            legend = legend_choice,
        )
    end
end

function lattice_plot_trinary(lattice,counts_joint)
    lattice_length = size(lattice,1)
    p = plot([0.5;lattice_length+0.5;lattice_length+0.5;0.5;0.5],[0.5;0.5;lattice_length+0.5;lattice_length+0.5;0.5],
        framestyle = :none,
        aspect_ratio = 1,
        line = :black,
        hover = false,
        legend = false
        )
    bin_length = Int(sqrt(size(counts_joint,1)-1))
    num_lines = Int((lattice_length/bin_length) - 1)
    for i=1:num_lines
        coord = (i)*bin_length+0.5
        plot!(p,[0, lattice_length],[coord, coord],
            linecolor = :black,
            linewidth = 2,
            hover = false,
            alpha = 0.4,
            legend = false
            )
        plot!(p,[coord, coord],[0, lattice_length],
            linecolor = :black,
            linewidth = 2,
            hover = false,
            alpha = 0.4,
            legend = false
            )
    end

    blue_locs = findall(x->x==1,lattice)
    red_locs = findall(x->x==2,lattice)
    green_locs = findall(x->x==3,lattice)

    sz = 0.4
    blue_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in blue_locs]
    blue_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in blue_locs]
    blue_vertices_x = reduce(vcat,blue_vertices_x)
    blue_vertices_y = reduce(vcat,blue_vertices_y)
    plot!(blue_vertices_x,blue_vertices_y,seriestype = :shape,
        hover = nothing,
        fillcolor = :blue,
        linewidth = 0,
        legend = false
        )

    red_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in red_locs]
    red_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in red_locs]
    red_vertices_x = reduce(vcat,red_vertices_x)
    red_vertices_y = reduce(vcat,red_vertices_y)
    plot!(red_vertices_x,red_vertices_y,seriestype = :shape,
        hover = nothing,
        fillcolor = :red,
        linewidth = 0,
        legend = false
        )

    if ~isempty(green_locs)
        green_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in green_locs]
        green_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in green_locs]
        green_vertices_x = reduce(vcat,green_vertices_x)
        green_vertices_y = reduce(vcat,green_vertices_y)
        plot!(green_vertices_x,green_vertices_y,seriestype = :shape,
            hover = nothing,
            fillcolor = :green,
            linewidth = 0,
            legend = false
            )
    end
    p

end

function headache_plot_trinary(counts_joint;a_label = "Blue",
    b_label = "Red",
    c_label = "Green")
    #a is the top, rows of counts_joint
    #b is the bottom-left, columns of counts_joint
    #c is the bottom-right, size of block - rows - columns of counts_joint

    H = calculate_H_trinary_function(counts_joint)
    data_arr = H

    max_occupancy = size(data_arr,1)-1;
    joint_density = [];
    data_vect = Array{Float64,1}(undef,0);
    for r = 1:size(data_arr,1)
        for c = 1:size(data_arr,2)
            if ~isnan(data_arr[r,c])
                joint_density = cat(joint_density,[r-1 c-1],dims = 1);
                append!(data_vect,[data_arr[r,c]]);
            end
        end
    end
    joint_frac = joint_density./max_occupancy;
    empty_frac = 1 .- joint_frac[:,1] .- joint_frac[:,2];
    delta_frac = 1/(max_occupancy);
    p = plot([0,0.5,1,0],[0,sin(pi/3),0,0],legend = false,
        aspect_ratio =  1,
        axis = nothing,
        grid = false,
        framestyle = :none,
        colorbar_entry = false,
        hover = false,
        colorbar = true
        )
    text_shift = 0.075
    annotate!(p,0-text_shift,0-text_shift,b_label,hover = false,colorbar = true)
    annotate!(p,0.5,sin(pi/3)+text_shift,a_label,hover = false,colorbar = true)
    annotate!(p,1+text_shift,0-text_shift,c_label,hover = false,colorbar = true)
    x_coords_all = Array{Float64,1}(undef,0)
    y_coords_all = Array{Float64,1}(undef,0)
    for i = 1:length(data_vect)
        a_coords =  [joint_frac[i,1]+delta_frac/2; #top-left
            joint_frac[i,1]+delta_frac/2; #top-right
            joint_frac[i,1]; #right
            joint_frac[i,1]-delta_frac/2; #bottom-right
            joint_frac[i,1]-delta_frac/2; #bottom-left
            joint_frac[i,1]; #left
            joint_frac[i,1]+delta_frac/2] #top-left
        b_coords = [joint_frac[i,2];#top-left
            joint_frac[i,2]-delta_frac/2;#top-right
            joint_frac[i,2]-delta_frac/2;#right
            joint_frac[i,2];#bottom-right
            joint_frac[i,2]+delta_frac/2;#bottom-left
            joint_frac[i,2]+delta_frac/2;#left
            joint_frac[i,2]]#top-left
        c_coords =  [empty_frac[i,1]-delta_frac/2;#top-left
            empty_frac[i];#top-right
            empty_frac[i]+delta_frac/2;#right
            empty_frac[i]+delta_frac/2;#bottom-right
            empty_frac[i];#bottom-left
            empty_frac[i]-delta_frac/2;#left
            empty_frac[i]-delta_frac/2]#top-left
        a_coords[a_coords.<0] .=0;b_coords[b_coords.<0] .=0;c_coords[c_coords.<0] .=0
        a_coords[a_coords.>1] .=1;b_coords[b_coords.>1] .=1;c_coords[c_coords.>1] .=1
        y_coords = sin(pi/3).*a_coords
        x_coords = cos(pi/3).*a_coords .+ c_coords
        for i = 1:length(x_coords)
            if a_coords[i]+c_coords[i]>1
                x_coords[i] = x_coords[i] - (a_coords[i]+c_coords[i]-1)
            end
        end
        append!(x_coords_all,[x_coords;NaN])
        append!(y_coords_all,[y_coords;NaN])
    end
    plot!(p,x_coords_all,y_coords_all,
        seriestype = :shape,
        fill_z = data_vect,
        seriescolor = :inferno,
        linealpha = 0,
        alpha = 1,
        # colorbar = true,
        legend = false,
        hover = false,
        colorbar = true
        )
    xlims!(p,(0-1*text_shift, 1+1*text_shift))
    ylims!(p,(0-1*text_shift, sin(pi/3)+2*text_shift))
    p
end

function probability_distribution_plot_trinary(counts_joint;a_label = "Blue",
    b_label = "Red",
    c_label = "Green")
    #a is the top, rows of counts_joint
    #b is the bottom-left, columns of counts_joint
    #c is the bottom-right, size of block - rows - columns of counts_joint

    P = counts_joint./sum(counts_joint)
    data_arr = P

    max_occupancy = size(data_arr,1)-1;
    joint_density = [];
    data_vect = Array{Float64,1}(undef,0);
    for r = 1:size(data_arr,1)
        for c = 1:size(data_arr,2)
            if ~isnan(data_arr[r,c])
                joint_density = cat(joint_density,[r-1 c-1],dims = 1);
                append!(data_vect,[data_arr[r,c]]);
            end
        end
    end
    joint_frac = joint_density./max_occupancy;
    empty_frac = 1 .- joint_frac[:,1] .- joint_frac[:,2];
    delta_frac = 1/(max_occupancy);
    p = plot([0,0.5,1,0],[0,sin(pi/3),0,0],legend = false,
        aspect_ratio =  1,
        axis = nothing,
        grid = false,
        framestyle = :none,
        colorbar_entry = false,
        hover = false,
        colorbar = true
        )
    text_shift = 0.075
    annotate!(p,0-text_shift,0-text_shift,b_label,10,hover = false,colorbar = true)
    annotate!(p,0.5,sin(pi/3)+text_shift,a_label,10,hover = false,colorbar = true)
    annotate!(p,1+text_shift,0-text_shift,c_label,10,hover = false,colorbar = true)
    x_coords_all = Array{Float64,1}(undef,0)
    y_coords_all = Array{Float64,1}(undef,0)
    for i = 1:length(data_vect)
        a_coords =  [joint_frac[i,1]+delta_frac/2; #top-left
            joint_frac[i,1]+delta_frac/2; #top-right
            joint_frac[i,1]; #right
            joint_frac[i,1]-delta_frac/2; #bottom-right
            joint_frac[i,1]-delta_frac/2; #bottom-left
            joint_frac[i,1]; #left
            joint_frac[i,1]+delta_frac/2] #top-left
        b_coords = [joint_frac[i,2];#top-left
            joint_frac[i,2]-delta_frac/2;#top-right
            joint_frac[i,2]-delta_frac/2;#right
            joint_frac[i,2];#bottom-right
            joint_frac[i,2]+delta_frac/2;#bottom-left
            joint_frac[i,2]+delta_frac/2;#left
            joint_frac[i,2]]#top-left
        c_coords =  [empty_frac[i,1]-delta_frac/2;#top-left
            empty_frac[i];#top-right
            empty_frac[i]+delta_frac/2;#right
            empty_frac[i]+delta_frac/2;#bottom-right
            empty_frac[i];#bottom-left
            empty_frac[i]-delta_frac/2;#left
            empty_frac[i]-delta_frac/2]#top-left
        a_coords[a_coords.<0] .=0;b_coords[b_coords.<0] .=0;c_coords[c_coords.<0] .=0
        a_coords[a_coords.>1] .=1;b_coords[b_coords.>1] .=1;c_coords[c_coords.>1] .=1
        y_coords = sin(pi/3).*a_coords
        x_coords = cos(pi/3).*a_coords .+ c_coords
        for i = 1:length(x_coords)
            if a_coords[i]+c_coords[i]>1
                x_coords[i] = x_coords[i] - (a_coords[i]+c_coords[i]-1)
            end
        end
        append!(x_coords_all,[x_coords;NaN])
        append!(y_coords_all,[y_coords;NaN])
    end
    plot!(p,x_coords_all,y_coords_all,
        seriestype = :shape,
        fill_z = data_vect,
        seriescolor = :blues,
        markerstrokealpha = 0,
        linealpha = 0,
        alpha = 1,
        # colorbar = true,
        legend = false,
        hover = false,
        colorbar = true
        )
    xlims!(p,(0-1*text_shift, 1+1*text_shift))
    ylims!(p,(0-1*text_shift, sin(pi/3)+2*text_shift))
    p
end

function transition_forecast_plot_trinary(counts_joint;
    initial_number_blue = 5,
    initial_number_red = 5,
    t,
    a_label = "Blue",
    b_label = "Red",
    c_label = "Green")
    #a is the top, rows of counts_joint
    #b is the bottom-left, columns of counts_joint
    #c is the bottom-right, size of block - rows - columns of counts_joint
    H = calculate_H_trinary_function(counts_joint)
    future_dist = forecast_trinary(H,initial_number_blue,initial_number_red,t)
    data_arr = future_dist

    max_occupancy = size(data_arr,1)-1;
    joint_density = [];
    data_vect = Array{Float64,1}(undef,0);
    for r = 1:size(data_arr,1)
        for c = 1:size(data_arr,2)
            if ~isnan(data_arr[r,c])
                joint_density = cat(joint_density,[r-1 c-1],dims = 1);
                append!(data_vect,[data_arr[r,c]]);
            end
        end
    end
    joint_frac = joint_density./max_occupancy;
    empty_frac = 1 .- joint_frac[:,1] .- joint_frac[:,2];
    delta_frac = 1/(max_occupancy);
    p = plot([0,0.5,1,0],[0,sin(pi/3),0,0],legend = false,
        aspect_ratio =  1,
        axis = nothing,
        grid = false,
        framestyle = :none,
        colorbar_entry = false,
        hover = false,
        colorbar = true
        )
    text_shift = 0.075
    annotate!(p,0-text_shift,0-text_shift,b_label,10,hover = false,colorbar = true)
    annotate!(p,0.5,sin(pi/3)+text_shift,a_label,10,hover = false,colorbar = true)
    annotate!(p,1+text_shift,0-text_shift,c_label,10,hover = false,colorbar = true)
    x_coords_all = Array{Float64,1}(undef,0)
    y_coords_all = Array{Float64,1}(undef,0)
    for i = 1:length(data_vect)
        a_coords =  [joint_frac[i,1]+delta_frac/2; #top-left
            joint_frac[i,1]+delta_frac/2; #top-right
            joint_frac[i,1]; #right
            joint_frac[i,1]-delta_frac/2; #bottom-right
            joint_frac[i,1]-delta_frac/2; #bottom-left
            joint_frac[i,1]; #left
            joint_frac[i,1]+delta_frac/2] #top-left
        b_coords = [joint_frac[i,2];#top-left
            joint_frac[i,2]-delta_frac/2;#top-right
            joint_frac[i,2]-delta_frac/2;#right
            joint_frac[i,2];#bottom-right
            joint_frac[i,2]+delta_frac/2;#bottom-left
            joint_frac[i,2]+delta_frac/2;#left
            joint_frac[i,2]]#top-left
        c_coords =  [empty_frac[i,1]-delta_frac/2;#top-left
            empty_frac[i];#top-right
            empty_frac[i]+delta_frac/2;#right
            empty_frac[i]+delta_frac/2;#bottom-right
            empty_frac[i];#bottom-left
            empty_frac[i]-delta_frac/2;#left
            empty_frac[i]-delta_frac/2]#top-left
        a_coords[a_coords.<0] .=0;b_coords[b_coords.<0] .=0;c_coords[c_coords.<0] .=0
        a_coords[a_coords.>1] .=1;b_coords[b_coords.>1] .=1;c_coords[c_coords.>1] .=1
        y_coords = sin(pi/3).*a_coords
        x_coords = cos(pi/3).*a_coords .+ c_coords
        for i = 1:length(x_coords)
            if a_coords[i]+c_coords[i]>1
                x_coords[i] = x_coords[i] - (a_coords[i]+c_coords[i]-1)
            end
        end
        append!(x_coords_all,[x_coords;NaN])
        append!(y_coords_all,[y_coords;NaN])
    end
    plot!(p,x_coords_all,y_coords_all,
        seriestype = :shape,
        fill_z = data_vect,
        seriescolor = :blues,
        markerstrokealpha = 0,
        linealpha = 0,
        alpha = 1,
        # colorbar = true,
        legend = false,
        hover = false,
        colorbar = true
        )
    xlims!(p,(0-1*text_shift, 1+1*text_shift))
    ylims!(p,(0-1*text_shift, sin(pi/3)+2*text_shift))
    p
end


function segregation_indices_plot_trinary(counts_joint)
    H_theil, D = calculate_seg_indices(counts_joint)
    bar(["Entropy","Dissimilarity"], [H_theil,D],
        legend = false
        )
    # plot(bar_plot)
end

function frustration_plot_trinary(counts_joint;a_label = "Blue",
    b_label = "Red",
    c_label = "Green")
    #a is the top, rows of counts_joint
    #b is the bottom-left, columns of counts_joint
    #c is the bottom-right, size of block - rows - columns of counts_joint
    N_max = size(counts_joint,1)-1
    H = calculate_H_trinary_function(counts_joint)
    N1 = collect(0:N_max) * ones(Int,1,N_max+1)
    n1 = N1./N_max
    N2 = ones(Int,N_max+1,1) * collect(0:N_max)'
    n2 = N2./N_max
    y_int,V1 = fit_linear(n1[:],H[:])
    y_int,V2 = fit_linear(n2[:],H[:])
    f = H - V1.*n1 - V2 .* n2
    data_arr = f

    max_occupancy = size(data_arr,1)-1;
    joint_density = [];
    data_vect = Array{Float64,1}(undef,0);
    for r = 1:size(data_arr,1)
        for c = 1:size(data_arr,2)
            if ~isnan(data_arr[r,c])
                joint_density = cat(joint_density,[r-1 c-1],dims = 1);
                append!(data_vect,[data_arr[r,c]]);
            end
        end
    end
    joint_frac = joint_density./max_occupancy;
    empty_frac = 1 .- joint_frac[:,1] .- joint_frac[:,2];
    delta_frac = 1/(max_occupancy);
    p = plot([0,0.5,1,0],[0,sin(pi/3),0,0],legend = false,
        aspect_ratio =  1,
        axis = nothing,
        grid = false,
        framestyle = :none,
        colorbar_entry = false,
        hover = false,
        colorbar = true
        )
    text_shift = 0.075
    annotate!(p,0-text_shift,0-text_shift,b_label,hover = false,colorbar = true)
    annotate!(p,0.5,sin(pi/3)+text_shift,a_label,hover = false,colorbar = true)
    annotate!(p,1+text_shift,0-text_shift,c_label,hover = false,colorbar = true)
    x_coords_all = Array{Float64,1}(undef,0)
    y_coords_all = Array{Float64,1}(undef,0)
    for i = 1:length(data_vect)
        a_coords =  [joint_frac[i,1]+delta_frac/2; #top-left
            joint_frac[i,1]+delta_frac/2; #top-right
            joint_frac[i,1]; #right
            joint_frac[i,1]-delta_frac/2; #bottom-right
            joint_frac[i,1]-delta_frac/2; #bottom-left
            joint_frac[i,1]; #left
            joint_frac[i,1]+delta_frac/2] #top-left
        b_coords = [joint_frac[i,2];#top-left
            joint_frac[i,2]-delta_frac/2;#top-right
            joint_frac[i,2]-delta_frac/2;#right
            joint_frac[i,2];#bottom-right
            joint_frac[i,2]+delta_frac/2;#bottom-left
            joint_frac[i,2]+delta_frac/2;#left
            joint_frac[i,2]]#top-left
        c_coords =  [empty_frac[i,1]-delta_frac/2;#top-left
            empty_frac[i];#top-right
            empty_frac[i]+delta_frac/2;#right
            empty_frac[i]+delta_frac/2;#bottom-right
            empty_frac[i];#bottom-left
            empty_frac[i]-delta_frac/2;#left
            empty_frac[i]-delta_frac/2]#top-left
        a_coords[a_coords.<0] .=0;b_coords[b_coords.<0] .=0;c_coords[c_coords.<0] .=0
        a_coords[a_coords.>1] .=1;b_coords[b_coords.>1] .=1;c_coords[c_coords.>1] .=1
        y_coords = sin(pi/3).*a_coords
        x_coords = cos(pi/3).*a_coords .+ c_coords
        for i = 1:length(x_coords)
            if a_coords[i]+c_coords[i]>1
                x_coords[i] = x_coords[i] - (a_coords[i]+c_coords[i]-1)
            end
        end
        append!(x_coords_all,[x_coords;NaN])
        append!(y_coords_all,[y_coords;NaN])
    end
    plot!(p,x_coords_all,y_coords_all,
        seriestype = :shape,
        fill_z = data_vect,
        seriescolor = :inferno,
        linealpha = 0,
        alpha = 1,
        # colorbar = true,
        legend = false,
        hover = false,
        colorbar = true
        )
    xlims!(p,(0-1*text_shift, 1+1*text_shift))
    ylims!(p,(0-1*text_shift, sin(pi/3)+2*text_shift))
    p
end

function vexation_plot_trinary(counts_joint;a_label = "Blue",
    b_label = "Red",
    c_label = "Green")
    #a is the top, rows of counts_joint
    #b is the bottom-left, columns of counts_joint
    #c is the bottom-right, size of block - rows - columns of counts_joint

    N_max = size(counts_joint,1)-1
    H = calculate_H_trinary_function(counts_joint)
    N1 = collect(0:N_max) * ones(Int,1,N_max+1)
    n1 = N1./N_max
    N2 = ones(Int,N_max+1,1) * collect(0:N_max)'
    n2 = N2./N_max
    y_int,V1 = fit_linear(n1[:],H[:])
    y_int,V2 = fit_linear(n2[:],H[:])
    data_arr = V1.*n1 + V2.*n2

    max_occupancy = size(data_arr,1)-1;
    joint_density = [];
    data_vect = Array{Float64,1}(undef,0);
    for r = 1:size(data_arr,1)
        for c = 1:size(data_arr,2)
            if ~isnan(data_arr[r,c])
                joint_density = cat(joint_density,[r-1 c-1],dims = 1);
                append!(data_vect,[data_arr[r,c]]);
            end
        end
    end
    joint_frac = joint_density./max_occupancy;
    empty_frac = 1 .- joint_frac[:,1] .- joint_frac[:,2];
    delta_frac = 1/(max_occupancy);
    p = plot([0,0.5,1,0],[0,sin(pi/3),0,0],legend = false,
        aspect_ratio =  1,
        axis = nothing,
        grid = false,
        framestyle = :none,
        colorbar_entry = false,
        hover = false,
        colorbar = true
        )
    text_shift = 0.075
    annotate!(p,0-text_shift,0-text_shift,b_label,hover = false,colorbar = true)
    annotate!(p,0.5,sin(pi/3)+text_shift,a_label,hover = false,colorbar = true)
    annotate!(p,1+text_shift,0-text_shift,c_label,hover = false,colorbar = true)
    x_coords_all = Array{Float64,1}(undef,0)
    y_coords_all = Array{Float64,1}(undef,0)
    for i = 1:length(data_vect)
        a_coords =  [joint_frac[i,1]+delta_frac/2; #top-left
            joint_frac[i,1]+delta_frac/2; #top-right
            joint_frac[i,1]; #right
            joint_frac[i,1]-delta_frac/2; #bottom-right
            joint_frac[i,1]-delta_frac/2; #bottom-left
            joint_frac[i,1]; #left
            joint_frac[i,1]+delta_frac/2] #top-left
        b_coords = [joint_frac[i,2];#top-left
            joint_frac[i,2]-delta_frac/2;#top-right
            joint_frac[i,2]-delta_frac/2;#right
            joint_frac[i,2];#bottom-right
            joint_frac[i,2]+delta_frac/2;#bottom-left
            joint_frac[i,2]+delta_frac/2;#left
            joint_frac[i,2]]#top-left
        c_coords =  [empty_frac[i,1]-delta_frac/2;#top-left
            empty_frac[i];#top-right
            empty_frac[i]+delta_frac/2;#right
            empty_frac[i]+delta_frac/2;#bottom-right
            empty_frac[i];#bottom-left
            empty_frac[i]-delta_frac/2;#left
            empty_frac[i]-delta_frac/2]#top-left
        a_coords[a_coords.<0] .=0;b_coords[b_coords.<0] .=0;c_coords[c_coords.<0] .=0
        a_coords[a_coords.>1] .=1;b_coords[b_coords.>1] .=1;c_coords[c_coords.>1] .=1
        y_coords = sin(pi/3).*a_coords
        x_coords = cos(pi/3).*a_coords .+ c_coords
        for i = 1:length(x_coords)
            if a_coords[i]+c_coords[i]>1
                x_coords[i] = x_coords[i] - (a_coords[i]+c_coords[i]-1)
            end
        end
        append!(x_coords_all,[x_coords;NaN])
        append!(y_coords_all,[y_coords;NaN])
    end
    plot!(p,x_coords_all,y_coords_all,
        seriestype = :shape,
        fill_z = data_vect,
        seriescolor = :inferno,
        linealpha = 0,
        alpha = 1,
        # colorbar = true,
        legend = false,
        hover = false,
        colorbar = true
        )
    xlims!(p,(0-1*text_shift, 1+1*text_shift))
    ylims!(p,(0-1*text_shift, sin(pi/3)+2*text_shift))
    p
end

function forecast_trinary(H,initial_number_blue,initial_number_red,t)
    #Let's build a big transition matrix
    N_max = size(H,1)-1
    Ns = Array{Int16,2}(undef,0,3)

    for N1 = 0:N_max
        for N2 = 0:N_max-N1
            N3 = N_max - N1 - N2
            Ns = [Ns;N1 N2 N3]
        end
    end
    ns = Ns./N_max
    T_matrix = zeros(size(Ns,1),size(Ns,1))
    for i = 1:size(Ns,1) #For each initial state
        #Calculate the probability of going to a future state in one step
        #For each state, there are 7 options,
        #Switch blue to red
        if Ns[i,1]>0 #Blue is chosen
            #BLUE TO RED
            dH_dN_b_r = N_max.*(H[Ns[i,1]+1-1,Ns[i,2]+1+1]-H[Ns[i,1]+1,Ns[i,2]+1])
            #BLUE TO OTHER
            dH_dN_b_o = N_max.*(H[Ns[i,1]+1-1,Ns[i,2]+1]-H[Ns[i,1]+1,Ns[i,2]+1])
            R_b_r = ns[i,1].*exp.(-dH_dN_b_r)./(1.0.+exp.(-dH_dN_b_r)+exp(-dH_dN_b_o));
            R_b_o = ns[i,1].*exp.(-dH_dN_b_o)./(1.0.+exp.(-dH_dN_b_r)+exp(-dH_dN_b_o));
            idx_b_r = [[Ns[i,1]-1;Ns[i,2]+1;Ns[i,3]] == Ns[j,:] for j = 1:size(Ns,1)]
            idx_b_o = [[Ns[i,1]-1;Ns[i,2];Ns[i,3]+1] == Ns[j,:] for j = 1:size(Ns,1)]
            T_matrix[i,idx_b_r] .= R_b_r
            T_matrix[i,idx_b_o] .= R_b_o
        else
            R_b_r = 0
            R_b_o = 0
        end
        if Ns[i,2]>0
            #RED TO BLUE
            dH_dN_r_b = N_max.*(H[Ns[i,1]+1+1,Ns[i,2]+1-1]-H[Ns[i,1]+1,Ns[i,2]+1])
            #RED TO OTHER
            dH_dN_r_o = N_max.*(H[Ns[i,1]+1,Ns[i,2]+1-1]-H[Ns[i,1]+1,Ns[i,2]+1])
            R_r_b = ns[i,2].*exp.(-dH_dN_r_b)./(1.0.+exp.(-dH_dN_r_b)+exp(-dH_dN_r_o));
            R_r_o = ns[i,2].*exp.(-dH_dN_r_o)./(1.0.+exp.(-dH_dN_r_b)+exp(-dH_dN_r_o));
            idx_r_b = [[Ns[i,1]+1;Ns[i,2]-1;Ns[i,3]] == Ns[j,:] for j = 1:size(Ns,1)]
            idx_r_o = [[Ns[i,1];Ns[i,2]-1;Ns[i,3]+1] == Ns[j,:] for j = 1:size(Ns,1)]
            T_matrix[i,idx_r_b] .= R_r_b
            T_matrix[i,idx_r_o] .= R_r_o
        else
            R_r_b = 0
            R_r_o = 0
        end
        if Ns[i,3]>0
            #OTHER TO BLUE
            dH_dN_o_b = N_max.*(H[Ns[i,1]+1+1,Ns[i,2]+1]-H[Ns[i,1]+1,Ns[i,2]+1])
            #OTHER TO RED
            dH_dN_o_r = N_max.*(H[Ns[i,1]+1,Ns[i,2]+1+1]-H[Ns[i,1]+1,Ns[i,2]+1])
            R_o_b = ns[i,3].*exp.(-dH_dN_o_b)./(1.0.+exp.(-dH_dN_o_b)+exp(-dH_dN_o_r));
            R_o_r = ns[i,3].*exp.(-dH_dN_o_r)./(1.0.+exp.(-dH_dN_o_b)+exp(-dH_dN_o_r));
            idx_o_b = [[Ns[i,1]+1;Ns[i,2];Ns[i,3]-1] == Ns[j,:] for j = 1:size(Ns,1)]
            idx_o_r = [[Ns[i,1];Ns[i,2]+1;Ns[i,3]-1] == Ns[j,:] for j = 1:size(Ns,1)]
            T_matrix[i,idx_o_b] .= R_o_b
            T_matrix[i,idx_o_r] .= R_o_r
        else
            R_o_b = 0
            R_o_r = 0
        end
        R_no_switch = 1 - R_b_r - R_b_o - R_r_b - R_r_o - R_o_b - R_o_r
        #Now assign them to the correct entries in the transition matrix.
        T_matrix[i,i] = R_no_switch
    end
    replace!(T_matrix,NaN=>0)
    initial_state = zeros(size(Ns,1),1)
    initial_number_other = N_max - initial_number_blue - initial_number_red
    idx_initial_state = [[initial_number_blue;initial_number_red;initial_number_other] == Ns[j,:] for j = 1:size(Ns,1)]
    initial_state[idx_initial_state] .= 1
    t = Int(t)
    final_state = (T_matrix^t)*initial_state
    future_dist = NaN.*similar(H)
    for i=1:size(Ns,1)
        future_dist[Ns[i,1]+1,Ns[i,2]+1] = final_state[i]
    end
    return future_dist
end

function rate_plot_trinary(counts_joint,
    switching_process
    ;
    a_label = "Blue",
    b_label = "Red",
    c_label = "Green")
    #a is the top, rows of counts_joint
    #b is the bottom-left, columns of counts_joint
    #c is the bottom-right, size of block - rows - columns of counts_joint
    N_max = size(counts_joint,1)-1
    H = calculate_H_trinary_function(counts_joint)

    if switching_process == 1 #Switching blue to red
        R = NaN.*similar(H)
        for N1=1:N_max
            for N2 = 0:N_max-1
                delta_H = N_max.*(H[N1+1-1,N2+1+1] .- H[N1+1,N2+1])
                R[N1+1,N2+1] = 1.0./(1.0.+exp.(delta_H))
            end
        end
    elseif switching_process == 2 #Swtiching blue to green/empty
        R = NaN.*similar(H)
        for N1=1:N_max
            for N2 = 0:N_max
                delta_H = N_max.*(H[N1+1-1,N2+1] .- H[N1+1,N2+1])
                R[N1+1,N2+1] = 1.0./(1.0.+exp.(delta_H))
            end
        end
    else #switching red to empty
        R = NaN.*similar(H)
        for N1=0:N_max
            for N2 = 1:N_max
                delta_H = N_max.*(H[N1+1,N2+1-1] .- H[N1+1,N2+1])
                R[N1+1,N2+1] = 1.0./(1.0.+exp.(delta_H))
            end
        end
    end
    data_arr = R

    max_occupancy = size(data_arr,1)-1;
    joint_density = [];
    data_vect = Array{Float64,1}(undef,0);
    for r = 1:size(data_arr,1)
        for c = 1:size(data_arr,2)
            if ~isnan(data_arr[r,c])
                joint_density = cat(joint_density,[r-1 c-1],dims = 1);
                append!(data_vect,[data_arr[r,c]]);
            end
        end
    end
    joint_frac = joint_density./max_occupancy;
    empty_frac = 1 .- joint_frac[:,1] .- joint_frac[:,2];
    delta_frac = 1/(max_occupancy);
    p = plot([0,0.5,1,0],[0,sin(pi/3),0,0],legend = false,
        aspect_ratio =  1,
        axis = nothing,
        grid = false,
        framestyle = :none,
        colorbar_entry = false,
        hover = false,
        colorbar = true
        )
    text_shift = 0.075
    annotate!(p,0-text_shift,0-text_shift,b_label,hover = false,colorbar = true)
    annotate!(p,0.5,sin(pi/3)+text_shift,a_label,hover = false,colorbar = true)
    annotate!(p,1+text_shift,0-text_shift,c_label,hover = false,colorbar = true)
    x_coords_all = Array{Float64,1}(undef,0)
    y_coords_all = Array{Float64,1}(undef,0)
    for i = 1:length(data_vect)
        a_coords =  [joint_frac[i,1]+delta_frac/2; #top-left
            joint_frac[i,1]+delta_frac/2; #top-right
            joint_frac[i,1]; #right
            joint_frac[i,1]-delta_frac/2; #bottom-right
            joint_frac[i,1]-delta_frac/2; #bottom-left
            joint_frac[i,1]; #left
            joint_frac[i,1]+delta_frac/2] #top-left
        b_coords = [joint_frac[i,2];#top-left
            joint_frac[i,2]-delta_frac/2;#top-right
            joint_frac[i,2]-delta_frac/2;#right
            joint_frac[i,2];#bottom-right
            joint_frac[i,2]+delta_frac/2;#bottom-left
            joint_frac[i,2]+delta_frac/2;#left
            joint_frac[i,2]]#top-left
        c_coords =  [empty_frac[i,1]-delta_frac/2;#top-left
            empty_frac[i];#top-right
            empty_frac[i]+delta_frac/2;#right
            empty_frac[i]+delta_frac/2;#bottom-right
            empty_frac[i];#bottom-left
            empty_frac[i]-delta_frac/2;#left
            empty_frac[i]-delta_frac/2]#top-left
        a_coords[a_coords.<0] .=0;b_coords[b_coords.<0] .=0;c_coords[c_coords.<0] .=0
        a_coords[a_coords.>1] .=1;b_coords[b_coords.>1] .=1;c_coords[c_coords.>1] .=1
        y_coords = sin(pi/3).*a_coords
        x_coords = cos(pi/3).*a_coords .+ c_coords
        for i = 1:length(x_coords)
            if a_coords[i]+c_coords[i]>1
                x_coords[i] = x_coords[i] - (a_coords[i]+c_coords[i]-1)
            end
        end
        append!(x_coords_all,[x_coords;NaN])
        append!(y_coords_all,[y_coords;NaN])
    end
    plot!(p,x_coords_all,y_coords_all,
        seriestype = :shape,
        fill_z = data_vect,
        seriescolor = :inferno,
        linealpha = 0,
        alpha = 1,
        # colorbar = true,
        legend = false,
        hover = false,
        colorbar = true
        )
    xlims!(p,(0-1*text_shift, 1+1*text_shift))
    ylims!(p,(0-1*text_shift, sin(pi/3)+2*text_shift))
    p
end


function intro_headache_trinary_figures(lattice_snapshots,counts_joint,utility_function_blue,utility_function_red,utility_function_green)
    p1 = utility_function_plot(utility_function_blue,utility_function_red,utility_function_green = utility_function_green)
    p2 = lattice_plot_trinary(lattice_snapshots[:,:,end],counts_joint)
    p3 = probability_distribution_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Green"
        )
    p4 = headache_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Green"
        )
    p5 = frustration_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Green"
        )
    p6 = vexation_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Green"
        )
    plot(p1,p2,p3,p4,p5,p6,
        layout = (2,3),
        title = ["Utility" "Simulation snapshot" "Histogram" "H(n1,n2,n3)" "f(n2,n2,n3)" "V1*n1+V2*n2"],
    )
end


function intro_headache_binary_w_vacancies_figures(lattice_snapshots,counts_joint,utility_function_blue,utility_function_red)
    p1 = utility_function_plot(utility_function_blue,utility_function_red)
    p2 = lattice_plot_trinary(lattice_snapshots[:,:,end],counts_joint)
    p3 = probability_distribution_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Empty"
        )
    p4 = headache_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Empty"
        )
    p5 = frustration_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Empty"
        )
    p6 = vexation_plot_trinary(counts_joint,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Empty"
        )
    plot(p1,p2,p3,p4,p5,p6,
        layout = (2,3),
        title = ["Utility" "Simulation snapshot" "Histogram" "H(n1,n2,n3)" "f(n2,n2,n3)" "V1*n1+V2*n2"],
    )
end

function interpretation_headache_binary_w_vacancies_figures(lattice_snapshots,
    counts_joint,
    utility_function_blue,
    utility_function_red,
    num_steps,
    initial_number_blue,
    initial_number_red,
    )
    p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_binary(lattice_snapshots[:,:,end],counts_single)
    p3 = probability_distribution_plot_trinary(counts_joint,
        a_label = "B",
        b_label = "R",
        c_label = "V"
        )
    #Utility,lattice,probability,Headache,Rate1,Rate2,Rate3
    #Evolution
    p4 = headache_plot_trinary(counts_joint,
        a_label = "B",
        b_label = "R",
        c_label = "V"
        )
    p5 = rate_plot_trinary(counts_joint,
        1;
        a_label = "B",
        b_label = "R",
        c_label = "V")
    p6 = rate_plot_trinary(counts_joint,
        2;
        a_label = "B",
        b_label = "R",
        c_label = "V")
    p7 = rate_plot_trinary(counts_joint,
        3;
        a_label = "B",
        b_label = "R",
        c_label = "V")
    p8 = transition_forecast_plot_trinary(counts_joint;
        initial_number_blue = initial_number_blue,
        initial_number_red = initial_number_red,
        t = num_steps,
        a_label = "B",
        b_label = "R",
        c_label = "V")

    l = @layout[grid(2,1) grid(2,1) a; b c d]

    plot(p1,p2,p3,p4,p5,p6,p7,p8,
        layout = l,
        title = ["Utility" "Snapshot" "Histogram" "H" "Rate b->r" "Rate b->v" "Rate r->v" string(num_steps)],
    )
end

function interpretation_headache_trinary_figures(lattice_snapshots,
    counts_joint,
    utility_function_blue,
    utility_function_red,
    utility_function_green,
    num_steps,
    initial_number_blue,
    initial_number_red,
    )
    p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_trinary(lattice_snapshots[:,:,end],counts_single)
    p3 = probability_distribution_plot_trinary(counts_joint,
        a_label = "B",
        b_label = "R",
        c_label = "G"
        )
    #Utility,lattice,probability,Headache,Rate1,Rate2,Rate3
    #Evolution
    p4 = headache_plot_trinary(counts_joint,
        a_label = "B",
        b_label = "R",
        c_label = "G"
        )
    p5 = rate_plot_trinary(counts_joint,
        1;
        a_label = "B",
        b_label = "R",
        c_label = "G")
    p6 = rate_plot_trinary(counts_joint,
        2;
        a_label = "B",
        b_label = "R",
        c_label = "G")
    p7 = rate_plot_trinary(counts_joint,
        3;
        a_label = "B",
        b_label = "R",
        c_label = "G")
    p8 = transition_forecast_plot_trinary(counts_joint;
        initial_number_blue = initial_number_blue,
        initial_number_red = initial_number_red,
        t = num_steps,
        a_label = "B",
        b_label = "R",
        c_label = "G")

    l = @layout[grid(2,1) grid(2,1) a; b c d]

    plot(p1,p2,p3,p4,p5,p6,p7,p8,
        layout = l,
        title = ["Utility" "Snapshot" "Histogram" "H" "Rate b->r" "Rate b->v" "Rate r->v" string(num_steps)],
    )
end


# function make_schelling_trinary_figures_vacancies(lattice_snapshots,counts_joint,utility_function_blue,utility_function_red)
#     p1 = utility_function_plot(utility_function_blue,utility_function_red)
#     p2 = lattice_plot_trinary(lattice_snapshots[:,:,end],counts_joint)
#     p3 = headache_plot_trinary(counts_joint,
#         a_label = "Blue",
#         b_label = "Red",
#         c_label = "Empty"
#         )
#     p4 = segregation_indices_plot_trinary(counts_joint)
#     plot(p1,p2,p3,p4,
#         layout = (2,2),
#         title = ["Utility" "Simulation snapshot" "Headache" "Segregation indices"],
#     )
# end
#
# function make_schelling_trinary_figures_three_agents(lattice_snapshots,
#     counts_joint,
#     utility_function_blue,
#     utility_function_red,
#     utility_function_green)
#     p1 = utility_function_plot(utility_function_blue,
#         utility_function_red,
#         utility_function_green = utility_function_green)
#     p2 = lattice_plot_trinary(lattice_snapshots[:,:,end],counts_joint)
#     p3 = headache_plot_trinary(counts_joint,
#         a_label = "Blue",
#         b_label = "Red",
#         c_label = "Green"
#         )
#     p4 = segregation_indices_plot(counts_joint)
#     plot(p1,p2,p3,p4,
#         layout = (2,2),
#         title = ["Utility" "Simulation snapshot" "Headache" "Segregation indices"],
#     )
# end

function plot_all_three_schelling_sims(
    lattice_snapshots_binary,counts_single_binary,utility_function_blue,utility_function_red,
    lattice_snapshots_binary_w_vacancies,counts_joint_binary_w_vacancies,utility_function_blue_w_vacancies,utility_function_red_w_vacancies,
    lattice_snapshots_trinary,counts_joint_trinary,utility_function_blue_trinary,utility_function_red_trinary,utility_function_green_trinary)

    p1 = utility_function_plot(utility_function_blue_binary,
        utility_function_red_binary,
        )
    p2 = utility_function_plot(utility_function_blue_binary_w_vacancies,
        utility_function_red_w_vacancies,
        )
    p3 = utility_function_plot(utility_function_blue_trinary,
        utility_function_red_trinary,
        utility_function_green = utility_function_green_trinary)

    p4 = lattice_plot_binary(lattice_snapshots_binary[:,:,end],counts_single_binary)
    p5 = lattice_plot_trinary(lattice_snapshots_binary_w_vacancies[:,:,end],counts_joint_binary_w_vacancies)
    p6 = lattice_plot_trinary(lattice_snapshots_trinary[:,:,end],counts_joint_trinary)

    plot(p1,p2,p3,p4,p5,p6,
        layout = (2,3),
        title = ["Binary" "Binary w/ Vacanices" "Trinary" "" "" ""],
    )
end

function segregation_indices_distributions_comparison_binary_w_vacancies(lattice_snapshots_A,counts_joint_A,utility_function_blue_A,utility_function_red_A,
    lattice_snapshots_B,counts_joint_B,utility_function_blue_B,utility_function_red_B,
    lattice_snapshots_C,counts_joint_C,utility_function_blue_C,utility_function_red_C
    )
    #Each
    p_utility_A = utility_function_plot(utility_function_blue_A,utility_function_red_A,legend_choice = false)
    p_utility_B = utility_function_plot(utility_function_blue_B,utility_function_red_B,legend_choice = false)
    p_utility_C = utility_function_plot(utility_function_blue_C,utility_function_red_C,legend_choice = false)
    xlabel!(p_utility_C,"# Neighbors")
    ylabel!(p_utility_A,"A")
    ylabel!(p_utility_B,"B")
    ylabel!(p_utility_C,"C")

    p_lattice_A = lattice_plot_trinary(lattice_snapshots_A[:,:,end],counts_joint_A)
    p_lattice_B = lattice_plot_trinary(lattice_snapshots_B[:,:,end],counts_joint_B)
    p_lattice_C = lattice_plot_trinary(lattice_snapshots_C[:,:,end],counts_joint_C)

    p_distribution_A = probability_distribution_plot_trinary(counts_joint_A,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Vacant")
    p_distribution_B = probability_distribution_plot_trinary(counts_joint_B,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Vacant")
    p_distribution_C = probability_distribution_plot_trinary(counts_joint_C,
        a_label = "Blue",
        b_label = "Red",
        c_label = "Vacant")

    H_theil_A, D_A = calculate_seg_indices(counts_joint_A)
    H_theil_B, D_B = calculate_seg_indices(counts_joint_B)
    H_theil_C, D_C = calculate_seg_indices(counts_joint_C)

    p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
        bar_position = :dodge,
        label = ["Entropy" "Dissimilarity"],
        legend = :topright
        )
    ylims!(p_seg_indices_bar,(0,1))
    xlabel!(p_seg_indices_bar,"City")

    l = @layout[grid(3,3) a{0.2w}]
    plot(p_utility_A,p_lattice_A,p_distribution_A,
        p_utility_B,p_lattice_B,p_distribution_B,
        p_utility_C,p_lattice_C,p_distribution_C,
        p_seg_indices_bar, layout = l,
        title = ["Utility" "Snapshot" "Histogram" "" "" "" "" "" "" "Seg. index"]
        )

    # l = @layout[grid(2,3) a{0.2w}]
    # plot(p_utility_A,p_lattice_A,p_distribution_A,
    #     p_utility_B,p_lattice_B,p_distribution_B,
    #     # p_utility_C,p_lattice_C,p_distribution_C,
    #     p_seg_indices_bar, layout = l,
    #     title = ["Utility" "Snapshot" "Histogram" "" "" "" "" "" "" "Seg. index"]
    #     )
end

function segregation_indices_distributions_comparison_trinary(lattice_snapshots_A,counts_joint_A,utility_function_blue_A,utility_function_red_A,utility_function_green_A,
    lattice_snapshots_B,counts_joint_B,utility_function_blue_B,utility_function_red_B,utility_function_green_B,
    lattice_snapshots_C,counts_joint_C,utility_function_blue_C,utility_function_red_C,utility_function_green_C
    )
    #Each
    p_utility_A = utility_function_plot(utility_function_blue_A,utility_function_red_A,utility_function_green = utility_function_green_A,legend_choice = false)
    p_utility_B = utility_function_plot(utility_function_blue_B,utility_function_red_B,utility_function_green = utility_function_green_B,legend_choice = false)
    p_utility_C = utility_function_plot(utility_function_blue_C,utility_function_red_C,utility_function_green = utility_function_green_C,legend_choice = false)
    xlabel!(p_utility_C,"# Neighbors")
    ylabel!(p_utility_A,"A")
    ylabel!(p_utility_B,"B")
    ylabel!(p_utility_C,"C")

    p_lattice_A = lattice_plot_trinary(lattice_snapshots_A[:,:,end],counts_joint_A)
    p_lattice_B = lattice_plot_trinary(lattice_snapshots_B[:,:,end],counts_joint_B)
    p_lattice_C = lattice_plot_trinary(lattice_snapshots_C[:,:,end],counts_joint_C)

    p_distribution_A = probability_distribution_plot_trinary(counts_joint_A,
        a_label = "B",
        b_label = "R",
        c_label = "G")
    p_distribution_B = probability_distribution_plot_trinary(counts_joint_B,
        a_label = "B",
        b_label = "R",
        c_label = "G")
    p_distribution_C = probability_distribution_plot_trinary(counts_joint_C,
        a_label = "B",
        b_label = "R",
        c_label = "G")

    H_theil_A, D_A = calculate_seg_indices(counts_joint_A)
    H_theil_B, D_B = calculate_seg_indices(counts_joint_B)
    H_theil_C, D_C = calculate_seg_indices(counts_joint_C)

    p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
        bar_position = :dodge,
        label = ["Entropy" "Dissimilarity"],
        legend = :topright
        )
    ylims!(p_seg_indices_bar,(0,1))
    xlabel!(p_seg_indices_bar,"City")

    l = @layout[grid(3,3) a{0.2w}]
    plot(p_utility_A,p_lattice_A,p_distribution_A,
        p_utility_B,p_lattice_B,p_distribution_B,
        p_utility_C,p_lattice_C,p_distribution_C,
        p_seg_indices_bar, layout = l,
        title = ["Utility" "Snapshot" "Histogram" "" "" "" "" "" "" "Seg. index"]
        )

    # l = @layout[grid(2,3) a{0.2w}]
    # plot(p_utility_A,p_lattice_A,p_distribution_A,
    #     p_utility_B,p_lattice_B,p_distribution_B,
    #     # p_utility_C,p_lattice_C,p_distribution_C,
    #     p_seg_indices_bar, layout = l,
    #     title = ["Utility" "Snapshot" "Histogram" "" "" "" "" "" "" "Seg. index"]
    #     )
end
