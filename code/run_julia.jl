using JSON3
using LinearAlgebra
using NPZ
using Plots
using Random

include("main.jl")
include("von_mises.jl")

Random.seed!(42)

function one_run(distribution, arange, e, grid, N, n_runs)
    n = size(arange, 1 )
    E, Z = zeros(n_runs, N, n), zeros(n_runs, N, n)
    accepted_theta, rejected_theta = zeros(n), zeros(n)

    for i in 1:n
        ## change size of distributon
        a_i = arange[i];
        b_i = a_i/e;

        for i in 1:size(distribution,1)
            distribution[i].a = a_i;
            distribution[i].b = b_i;
            fix_A!(distribution[i]);
        end

        ## distance matrix between all the centres
        S = [x.center for x in distribution]
        distance = pairwise(PeriodicEuclidean(grid), S);

        ## matrix of neighbouring ellipses, which are less than 2*a apart
        in_proximity = findall(x -> (x<=2*a_i) & (x!=0), distance);

        ## get vicinities of all ellipses
        results = [vicinity(i, in_proximity, distribution) for i in 1:N];
        (vicinities, indices) = (getindex.(results, 1), getindex.(results, 2));

        ## find new state
        (E[:,:,i], Z[:,:,i], acc, rej) = metropolis(distribution, in_proximity, grid, n_runs, 0);
        accepted_theta[i] = size(acc, 1)
        rejected_theta[i] = size(rej, 1)
    end

    # return E[end,:,:], Z[end,:,:], accepted_theta, rejected_theta
    return E[end,:,:], Z[end,:,:], accepted_theta, rejected_theta
end

## a0 is starting ellipse size, e is ration between semi axis, N number of ellipses,
## grid is grid size, n_runs is the number of runs of ellipse's orientational relaxation
## n_metropol is the number of rotations while doing the monte carlo, da is the growth
## step and a_max is the maximum size of ellipses
function many_runs(a0 = 4., e = 5/2, da = 0.1, a_max = 8., N = 128, grid = [100,100], n_runs = 10, n_metropol = 1000)
    initial = hcat( [rand()*grid[1], rand()*grid[2]], [rand()*grid[1], rand()*grid[2]]);
    samples = generate_samples(initial, grid, N);

    ## change circles into ellipses
    b0 = a0/e;                          # small semi-axis

    ## set a_range for growing ellipses
    # a_range = [15.: 0.1: 30. ...]
    a_range = [a0: da:a_max ...];
    n_a = size(a_range, 1);

    # E, Z = zeros(n_metropol, n_a, n_runs, n_a), zeros(n_metropol, n_runs, n_a);
    EN, Coord = Matrix{Float64}[], Matrix{Float64}[];
    accepted_theta, rejected_theta = zeros(n_a, n_runs), zeros(n_a, n_runs);
    correlation = zeros(100, n_runs)

    for i in 1:n_runs
        Random.seed!(i);

        ## set distribution
        distribution = [Ellipse(convert(Vector{Float64}, vector), a0, b0) for vector in eachcol(samples)];

        (E, Z, accepted_theta[:,i], rejected_theta[:,i]) = one_run(distribution, a_range, e, grid, N, n_metropol);
        push!(EN, E);
        push!(Coord, Z);
        correlation[:,i] = angle_correlation(distribution, 3*a_max, grid)
        
        # dists[i] = distribution
    end

    ## write data into files
#    data = Dict("distribution" => dists);
#    JSON3.write("one_run.json", data);

    ## do a bit of statistic
    all_energy = [0.5*sum(EN[i], dims=1) for i in 1:size(EN,1)];
    return_energy = mean(all_energy);
    err_energy = std(all_energy);

    all_coord = [sum(Coord[i], dims=1)/size(Coord[1][:,1], 1) for i in 1:size(Coord,1)];
    return_coord = mean(all_coord);
    err_coord = std(all_energy);

    all_accepted = mean(accepted_theta, dims=2)';
    err_accept = std(accepted_theta, dims=2)';
    all_rejected = mean(rejected_theta, dims=2)';
    err_reject = std(rejected_theta, dims=2)';

    all_correlation = mean(correlation, dims=2)';
    err_correlation = std(correlation, dims=2)';
    print(size(a_range));
    print(size(all_accepted));
    print(size(return_energy));
    print(size(return_coord));
    print(size(all_correlation));
 
    npzwrite("run_100.npz", 
        Dict("a_range" => a_range, "rejected_theta" => all_rejected, 
             "error_rejected" => err_reject, "accepted_theta" => all_accepted, 
             "error_accepted" => err_accept, "E" => return_energy, 
             "error_energy" => err_energy, "Z" => return_coord, "error_coord" => err_coord, 
             "correlation" => all_correlation, "error_correlation" => err_correlation, 
             "grid" => grid)
        );
end