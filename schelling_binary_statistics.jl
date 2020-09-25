#using SpecialFunctions
#using Polynomials
#using StatsBase
# using NaNMath

# function calculate_probability_mass(densities)
#
#     return probability_mass
# end
function calculate_H_binary(counts_single)
    N_max = length(counts_single)-1
    P = counts_single./(sum(counts_single))
    # H = [-log(P[N1+1,N2+1]*factorial(big(N1))*factorial(big(N2))/factorial(big(N_max))) for N1 = 0:N_max, N2 = 0:N_max]
    H = [-log(P[N+1]) - loggamma(N+1) - loggamma(N_max-N+1) + loggamma(N_max+1) for N = 0:N_max]
    H = H./N_max
    replace!(H, Inf => NaN)
    return H
end

function calculate_f_V_binary(counts_single_arr)
    N_max = size(counts_single,1)-1
    P = similar(counts_single_arr)
    for i=1:size(P,2)
        P[:,i] = counts_single_arr[:,i]./(sum(counts_single_arr[:,i]))
    end
    # H = [-log(P[N1+1,N2+1]*factorial(big(N1))*factorial(big(N2))/factorial(big(N_max))) for N1 = 0:N_max, N2 = 0:N_max]
    H = [-log(P[N+1]) - loggamma(N+1) + loggamma(N_max+1) for N = 0:N_max]
    H = H./N_max
    replace!(H, Inf => NaN)
    return H
end

function fit_linear(x,y)
    # M = [ones(length(x),1) x]
    # sol = y\M
    idx = .!isnan.(y)
    x = x[idx]
    y = y[idx]
    x_shift = x .- mean(x)
    y_shift = y .- mean(y)
    slope = (x_shift'*x_shift)\x_shift'*y_shift
    y_int = mean(y) - mean(x)*slope
    # slope = sol[2]
    # println()
    return y_int,slope
end

function calculate_seg_indices_binary(counts_single)
    # CHECK THIS CHECK THIS CHECK THIS!!!!
    # D is the dissimilarity index calculated assuming empty sites are truly empty houses with no people
    N_max = length(counts_single)-1
    num_blue_arr = [N*counts_single[N+1] for N = 0:N_max]
    num_red_arr = [(N_max-N)*counts_single[N+1] for N = 0:N_max]
    # D = 0.5.*sum(abs.(num_blue_arr./sum(num_blue_arr) - num_red_arr./sum(num_red_arr)))
    D = 0.5.*sum(abs.(num_blue_arr./sum(num_blue_arr) - num_red_arr./sum(num_red_arr)))
    # H_their is the Theil entropy index calculated as if empty sites are a type of other person
    #Uses equation from Table 1 from Reardon and Firebaugh 2002
    # num_empty_arr = [(N_max - N1 - N2)*counts_joint[N1+1,N2+1] for N1 = 0:N_max, N2 = 0:N_max]
    total_blue = sum(num_blue_arr)
    total_red = sum(num_red_arr)
    # total_empty = sum(num_empty_arr)
    # total = total_blue+total_red+total_empty
    total = total_blue+total_red
    prop_blue = total_blue/total
    prop_red = total_red/total
    # prop_empty = total_empty/total
    E = prop_blue * log(1/prop_blue) + prop_red * log(1/prop_red)

    right_sum_blue = 0
    for N = 0:N_max
        if N>0
            right_sum_blue+=counts_single[N+1]*N_max/total*N/N_max/prop_blue*log(N/N_max/prop_blue)
        end
    end
    right_sum_red = 0
    for N = 0:N_max
        if N>0
            right_sum_red+=counts_single[N_max-N+1]*N_max/total*N/N_max/prop_red*log(N/N_max/prop_red)
        end
    end

    # right_sum_empty = 0
    # for N1 = 0:N_max, N2 = 0:N_max
    #     if N1+N2<N_max;
    #         right_sum_empty+=counts_joint[N1+1,N2+1]*N_max/total*(N_max-N1-N2)/N_max/prop_empty*log((N_max-N1-N2)/N_max/prop_empty)
    #     end
    # end
    H_theil = 1/E * (prop_blue*right_sum_blue + prop_red*right_sum_red)
    return (H_theil,D)
end

# H_theil, D = calculate_seg_indices(counts_joint)
#
# #Testing indices with very segregated counts_joint
# counts_joint_seg = 1000*[1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 0]
# H_theil_seg, D_seg = calculate_seg_indices(counts_joint_seg)
# #Testing indices with very integregated counts_joint
# counts_joint_int = 1000*[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]
# H_theil_int, D_int = calculate_seg_indices(counts_joint_int)
