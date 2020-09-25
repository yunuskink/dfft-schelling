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

function lattice_plot_binary(lattice,counts_single)
    lattice_length = size(lattice,1)
    p = plot([0.5;lattice_length+0.5;lattice_length+0.5;0.5;0.5],[0.5;0.5;lattice_length+0.5;lattice_length+0.5;0.5],
        framestyle = :none,
        aspect_ratio = 1,
        line = :black,
        hover = false,
        legend = false
        )
    bin_length = Int(sqrt(length(counts_single)-1))
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
    p

end

function headache_plot_binary(counts_single;
    x_label = "# blue agents",
    linecolor = :blue,
    legend_labels = [""]
    )
    #counts_single is an array with each column representing a different experiment
    H = Array{Float64,2}(undef,length(counts_single),0)
    V_n = Array{Float64,2}(undef,length(counts_single),0)
    f_n = Array{Float64,2}(undef,length(counts_single),0)
    N_max = length(counts_single)-1
    n_s = 0:N_max
    n_s = collect(n_s)./N_max
    H_tmp = calculate_H_binary(counts_single)
    y_int, V_i = fit_linear(n_s,H_tmp)
    H_tmp = H_tmp .- y_int
    y_int, V_i = fit_linear(n_s,H_tmp)
    H = hcat(H,H_tmp)
    V_tmp = V_i.*n_s
    V_n = hcat(V_n,V_tmp)
    f_tmp = H_tmp.-V_tmp
    f_n = hcat(f_n,f_tmp)

    p = plot(H,label = "H",linecolor = linecolor,
        linewidth = 2)
    plot!(p,f_n,linestyle = :dash,label = "f",linecolor = linecolor,
        linewidth = 2)
    plot!(p,V_n,linestyle = :dot,label = "Vn",linecolor = linecolor,
        linewidth = 2)

    xlabel!(p,x_label)
    p
end

function rate_plot_binary(counts_single;
    x_label = "# blue agents",
    linecolor = :blue,
    legend_labels = [""]
    )
    #counts_single is an array with each column representing a different experiment
    H = calculate_H_binary(counts_single)
    N_max = length(counts_single)-1
    delta_H = N_max.*(H[1:end-1] .- H[2:end])
    R = 1.0./(1.0.+exp.(delta_H))
    p = plot(R,label = "H",
        linecolor = linecolor,
        linewidth = 2)
    xlabel!(p,x_label)
    p
end

function transition_matrix_plot_binary(counts_single,
    t
    )
    #counts_single is an array with each column representing a different experiment
    H = calculate_H_binary(counts_single)
    N_max = length(counts_single)-1
    dH_dN_increasing = N_max.*(H[2:end] .- H[1:end-1])
    dH_dN_decreasing = N_max.*(H[1:end-1] .- H[2:end])
    n_s = collect(0:N_max)./N_max
    decreasing_N_rate = n_s[2:end].*exp.(-dH_dN_decreasing)./(1.0.+exp.(-dH_dN_decreasing));
    increasing_N_rate = (1.0.-n_s[1:end-1]).*exp.(-dH_dN_increasing)./(1.0.+exp.(-dH_dN_increasing));
    constant_N_rate = 1.0 .- [0;decreasing_N_rate] .- [increasing_N_rate;0];
    T_matrix = zeros(size(H,1),size(H,1))
    for i=1:length(decreasing_N_rate)
        T_matrix[i+1,i] = decreasing_N_rate[i]
    end
    for i=1:length(increasing_N_rate)
        T_matrix[i,i+1] = increasing_N_rate[i]
    end
    for i=1:length(constant_N_rate)
        T_matrix[i,i] = constant_N_rate[i]
    end
    # T_matrix = diag(constant_N_rate) + diag(increasing_N_rate,1) + diag(decreasing_N_rate,-1);
    # println(increasing_N_rate)
    # println(constant_N_rate)
    ns_label = [string(i) for i = 0:N_max]
    t = Int(t)
    T_matrix = T_matrix^t
    T_matrix = rotl90(T_matrix)
    T_matrix = reverse(T_matrix,dims = 1)
    # T_matrix = reverse(T_matrix,dims = 1)
    # T_matrix = reverse(T_matrix,dims = 2)
    p = heatmap(ns_label,ns_label,T_matrix,aspect_ratio = 1)

    # p = plot(R,label = "H",
    #     linecolor = linecolor,
    #     linewidth = 2)
    # xlabel!(p,x_label)
    p
end


function segregation_indices_plot_binary(counts_single)
    H_theil, D = calculate_seg_indices_binary(counts_single)
    bar(["Entropy","Dissimilarity"], [H_theil,D],
        legend = false
        )
    # plot(bar_plot)
end

function intro_headache_binary_figures(lattice_snapshots,
    counts_single,
    utility_function_blue,
    utility_function_red)
    p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_binary(lattice_snapshots[:,:,end],counts_single)
    p3 = bar([0:(length(counts_single)-1)],counts_single,legend = false)
    p4 = headache_plot_binary(counts_single)

    plot(p1,p2,p3,p4,
        layout = (2,2),
        title = ["Utility" "Snapshot" "Histogram" "DFFT"],
    )
end

function interpretation_headache_binary_figures(lattice_snapshots,
    counts_single,
    utility_function_blue,
    utility_function_red,
    t)
    p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_binary(lattice_snapshots[:,:,end],counts_single)
    p3 = bar([0:(length(counts_single)-1)],counts_single,legend = false)
    p4 = headache_plot_binary(counts_single)
    p5 = rate_plot_binary(counts_single)
    p6 = transition_matrix_plot_binary(counts_single,t)

    plot(p1,p2,p3,p4,p5,p6,
        layout = (2,3),
        title = ["Utility" "Snapshot" "Histogram" "DFFT" "Rate b->r" string(t)],
    )
end


function compositional_invariance_binary_figure(lattice_snapshots_1,counts_single_1,fraction_blue_1,
    lattice_snapshots_2,counts_single_2,fraction_blue_2,
    lattice_snapshots_3,counts_single_3,fraction_blue_3,
    utility_function_blue,utility_function_red)

    p_utility = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    xlabel!(p_utility,"# Neighbors")
    # ylabel!(p_utility_A,"A")
    # ylabel!(p_utility_B,"B")
    # ylabel!(p_utility_C,"C")

    p_lattice_1 = lattice_plot_binary(lattice_snapshots_1[:,:,end],counts_single_1)
    p_lattice_2 = lattice_plot_binary(lattice_snapshots_2[:,:,end],counts_single_2)
    p_lattice_3 = lattice_plot_binary(lattice_snapshots_3[:,:,end],counts_single_3)

    p_distribution_1 = bar([0:(length(counts_single_1)-1)],counts_single_1,legend = false)
    p_distribution_2 = bar([0:(length(counts_single_2)-1)],counts_single_2,legend = false)
    p_distribution_3 = bar([0:(length(counts_single_3)-1)],counts_single_3,legend = false)
    xlabel!(p_distribution_3,"# Blue")

    H_1 = calculate_H_binary(counts_single_1)
    H_2 = calculate_H_binary(counts_single_2)
    H_3 = calculate_H_binary(counts_single_3)
    p_headaches_separate = plot([0:(length(H_1)-1)],
        [H_1 H_2 H_3],
        label = [fraction_blue_1 fraction_blue_2 fraction_blue_3],
        )

    H_arr = [H_1 H_2 H_3]
    H_mean = similar(H_1)
    for i=1:size(H_arr,1)
        H_tmp = 0
        num_not_nan = 0
        for j=1:size(H_arr,2)
            if ~isnan(H_arr[i,j])
                H_tmp += H_arr[i,j]
                num_not_nan += 1
            end
        end
        H_mean[i] = H_tmp./num_not_nan
    end
    N_max = length(H_1)-1
    C_1,V_1 = fit_linear((0:N_max),H_1.-H_mean)
    C_2,V_2 = fit_linear((0:N_max),H_2.-H_mean)
    C_3,V_3 = fit_linear((0:N_max),H_3.-H_mean)

    f_1 = H_1 .- V_1.*(0:N_max)
    f_2 = H_2 .- V_2.*(0:N_max)
    f_3 = H_3 .- V_3.*(0:N_max)

    f_2_diff = f_2.-f_1
    f_3_diff = f_3.-f_1
    f_2 = f_2 .- mean(f_2_diff[.!isnan.(f_2_diff)])
    f_3 = f_3 .- mean(f_3_diff[.!isnan.(f_3_diff)])

    p_DFFT_decomposition = plot(0:N_max,
        [f_1 f_2 f_3],
        label = [fraction_blue_1 fraction_blue_2 fraction_blue_3]
        )
    xlabel!(p_DFFT_decomposition,"# Blue")
    # plot!(p_DFFT_decomposition,
    #     0:N_max,
    #     [collect(C_1.+V_1.*(0:N_max)) collect(C_2.+V_2.*(0:N_max)) collect(C_3.+V_3.*(0:N_max))],
    #     )
    # p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
    #     bar_position = :dodge,
    #     label = ["Entropy" "Dissimilarity"],
    #     legend = :topright
    #     )
    # ylims!(p_seg_indices_bar,(0,1))
    # xlabel!(p_seg_indices_bar,"City")

    l = @layout[a{0.2w} grid(3,2) grid(2,1)]
    plot(p_utility,
        p_lattice_1,p_distribution_1,
        p_lattice_2,p_distribution_2,
        p_lattice_3,p_distribution_3,
        p_headaches_separate,
        p_DFFT_decomposition,
        layout = l,
        title = ["Utility" "Snapshot" "Histogram" "" "" "" "" "H(N_b)" "f(N_b)"]
        )

    #Show the following figures...
    #1) Utility functions
    #2) Snapshots of the three cities
    #3) Three Headache functions
    #4) Three aligned frustration functions and three linear vexation parts
end

function sample_size_invariance_binary_figure(lattice_snapshots_1,counts_single_1,fraction_blue_1,
    lattice_snapshots_2,counts_single_2,fraction_blue_2,
    lattice_snapshots_3,counts_single_3,fraction_blue_3,
    utility_function_blue,utility_function_red)

    p_utility = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    xlabel!(p_utility,"# Neighbors")
    # ylabel!(p_utility_A,"A")
    # ylabel!(p_utility_B,"B")
    # ylabel!(p_utility_C,"C")

    p_lattice_1 = lattice_plot_binary(lattice_snapshots_1[:,:,end],counts_single_1)
    p_lattice_2 = lattice_plot_binary(lattice_snapshots_2[:,:,end],counts_single_2)
    p_lattice_3 = lattice_plot_binary(lattice_snapshots_3[:,:,end],counts_single_3)

    p_distribution_1 = bar([0:(length(counts_single_1)-1)],counts_single_1,legend = false)
    p_distribution_2 = bar([0:(length(counts_single_2)-1)],counts_single_2,legend = false)
    p_distribution_3 = bar([0:(length(counts_single_3)-1)],counts_single_3,legend = false)
    xlabel!(p_distribution_3,"# Blue")

    H_1 = calculate_H_binary(counts_single_1)
    H_1 = H_1 .- mean(H_1[.!isnan.(H_1)])
    H_2 = calculate_H_binary(counts_single_2)
    H_2 = H_2 .- mean(H_2[.!isnan.(H_2)])
    H_3 = calculate_H_binary(counts_single_3)
    H_3 = H_3 .- mean(H_3[.!isnan.(H_3)])
    N_max_1 = length(counts_single_1)-1
    N_max_2 = length(counts_single_2)-1
    N_max_3 = length(counts_single_3)-1
    p_headaches_separate = plot([0:(length(H_1)-1)]./N_max_1,
        [H_1],
        label = [bin_length_1],
        )
    plot!(p_headaches_separate,
        [0:(length(H_2)-1)]./N_max_2,
        [H_2],
        label = [bin_length_2],
        )
    plot!(p_headaches_separate,
        [0:(length(H_3)-1)]./N_max_3,
        [H_3],
        label = [bin_length_3],
        )

    xlabel!(p_headaches_separate,"Fraction Blue")
    # plot!(p_DFFT_decomposition,
    #     0:N_max,
    #     [collect(C_1.+V_1.*(0:N_max)) collect(C_2.+V_2.*(0:N_max)) collect(C_3.+V_3.*(0:N_max))],
    #     )
    # p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
    #     bar_position = :dodge,
    #     label = ["Entropy" "Dissimilarity"],
    #     legend = :topright
    #     )
    # ylims!(p_seg_indices_bar,(0,1))
    # xlabel!(p_seg_indices_bar,"City")

    l = @layout[a{0.2w} grid(3,2) a]
    plot(p_utility,
        p_lattice_1,p_distribution_1,
        p_lattice_2,p_distribution_2,
        p_lattice_3,p_distribution_3,
        p_headaches_separate,
        layout = l,
        title = ["Utility" string(bin_length_1,"x",bin_length_1) "Histogram" string(bin_length_2,"x",bin_length_2) "" string(bin_length_3,"x",bin_length_3) "" "H(n_b)"]
        )

    #Show the following figures...
    #1) Utility functions
    #2) Snapshots of the three cities
    #3) Three Headache functions
    #4) Three aligned frustration functions and three linear vexation parts
end


function segregation_indices_distributions_comparison_binary_figure(lattice_snapshots_A,counts_single_A,utility_function_blue_A,utility_function_red_A,
    lattice_snapshots_B,counts_single_B,utility_function_blue_B,utility_function_red_B,
    lattice_snapshots_C,counts_single_C,utility_function_blue_C,utility_function_red_C
    )
    #Each
    p_utility_A = utility_function_plot(utility_function_blue_A,utility_function_red_A,legend_choice = false)
    p_utility_B = utility_function_plot(utility_function_blue_B,utility_function_red_B,legend_choice = false)
    p_utility_C = utility_function_plot(utility_function_blue_C,utility_function_red_C,legend_choice = false)
    xlabel!(p_utility_C,"# Neighbors")
    ylabel!(p_utility_A,"A")
    ylabel!(p_utility_B,"B")
    ylabel!(p_utility_C,"C")

    p_lattice_A = lattice_plot_binary(lattice_snapshots_A[:,:,end],counts_single_A)
    p_lattice_B = lattice_plot_binary(lattice_snapshots_B[:,:,end],counts_single_B)
    p_lattice_C = lattice_plot_binary(lattice_snapshots_C[:,:,end],counts_single_C)

    p_distribution_A = bar([0:(length(counts_single_C)-1)],counts_single_A,legend = false)
    p_distribution_B = bar([0:(length(counts_single_C)-1)],counts_single_B,legend = false)
    p_distribution_C = bar([0:(length(counts_single_C)-1)],counts_single_C,legend = false)
    xlabel!(p_distribution_C,"# Blue")

    H_theil_A, D_A = calculate_seg_indices_binary(counts_single_A)
    H_theil_B, D_B = calculate_seg_indices_binary(counts_single_B)
    H_theil_C, D_C = calculate_seg_indices_binary(counts_single_C)

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
end
