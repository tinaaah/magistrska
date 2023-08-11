using ArgParse, JSON3, LinearAlgebra, NPZ, Plots, Random
using Base.Iterators: flatten

include("main.jl")
include("von_mises.jl")

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
		"--a_0"
		help = "starting ellipse size"
		arg_type = Float64
		required = true
	
		"--a_max"
		help = "maximum size of ellipses"
		arg_type = Float64
		required = true
	
		"--d_a"
		help = "growth step"
		arg_type = Float64
		default = 0.1
	
		"--e"
		help = "ratio between semi axi"
		arg_type = Float64
		required = true
	
		"--N"
		help = "number of ellipses"
		arg_type = Int
		default = 1024
	
		"--x_size"
		help = "x length grid size"
		arg_type = Int
		default = 1000

		"--y_size"
		help = "y width grid size"
		arg_type = Int
		default = 1000

		"--n_runs"
		help = "number of growth iterations from the same starting system"
		arg_type = Int
		default = 100
	
		"--n_metropol"
		help = "number of rotations while performing the monte carlo"
		arg_type = Int
		default = 10000
	
		"--seed"
		help = "seed for randomizer"
		arg_type = Int
		default = 42
	end
	return parse_args(s)
end

function one_run(distribution, arange, e, grid, N, n_metropol)
    n = size(arange, 1 )
    # E, Z = zeros(n_runs, N, n), zeros(n_runs, N, n)
    E, Z = zeros(N, n), zeros(N, n)
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
        S = [x.center for x in distribution];
        distance = pairwise(PeriodicEuclidean(grid), S);

        ## matrix of neighbouring ellipses, which are less than 2*a apart
        in_proximity = findall(x -> (x<=2*a_i) & (x!=0), distance);

        ## get vicinities of all ellipses
        results = [vicinity(i, in_proximity, distribution) for i in 1:N];
        (vicinities, indices) = (getindex.(results, 1), getindex.(results, 2));

        ## find new state
        (E[:,i], Z[:,i], acc, rej) = metropolis2(distribution, in_proximity, grid, n_metropol, 0);
        accepted_theta[i] = size(acc, 1);
        rejected_theta[i] = size(rej, 1);
    end
    # return E, Z, Iaccepted_theta, rejected_theta
    return vec(sum(E, dims=1))*0.5, vec(sum(Z, dims=1))*0.5, accepted_theta, rejected_theta
    gc()
end

function many_runs(a0, e, da, a_max, N, grid, n_runs, n_metropol, seed0)
    initial = hcat( [rand()*grid[1], rand()*grid[2]], [rand()*grid[1], rand()*grid[2]]);
    samples = generate_samples(initial, grid, N);

    ## change circles into ellipses
	b0 = a0/e;                          # small semi-axis

    ## set a_range for growing ellipses
    a_range = [a0: da:a_max ...];
    n_a = size(a_range, 1);

    ## for multithreading
    EN = [zeros(n_a) for _ in 1:Threads.nthreads()];
    EN_squared = [zeros(n_a) for _ in 1:Threads.nthreads()];
    Coord = [zeros(n_a) for _ in 1:Threads.nthreads()];
    Coord_squared = [zeros(n_a) for _ in 1:Threads.nthreads()];

    accepted_theta = [zeros(n_a) for _ in 1:Threads.nthreads()];
    accepted_squared = [zeros(n_a) for _ in 1:Threads.nthreads()];

    rejected_theta = [zeros(n_a) for _ in 1:Threads.nthreads()];
    rejected_squared = [zeros(n_a) for _ in 1:Threads.nthreads()];
    
    correlation = [zeros(100) for _ in 1:Threads.nthreads()];
    correlation_squared = [zeros(100) for _ in 1:Threads.nthreads()];

    Threads.@threads for i in 1:n_runs
    # for i in 1:n_runs
        ## every run start from a distribution with same orientations
        Random.seed!(seed0);
		distribution = [Ellipse(convert(Vector{Float64}, vector), a0, b0) for vector in eachcol(samples)];

        ## deterministic random so that it can be recreated
        Random.seed!(i);

        (e_i, z_i, accept_i, reject_i) = one_run(distribution, a_range, e, grid, N, n_metropol);

        EN[Threads.threadid()] += e_i;
        EN_squared[Threads.threadid()] += e_i.^2;

        Coord[Threads.threadid()] += z_i;
        Coord_squared[Threads.threadid()] += z_i.^2;

        accepted_theta[Threads.threadid()] += accept_i;
        accepted_squared[Threads.threadid()] += accept_i.^2;

        rejected_theta[Threads.threadid()] += reject_i;
        rejected_squared[Threads.threadid()] += reject_i.^2;

        corr = angle_correlation(distribution, 3*a_max, grid);
        correlation[Threads.threadid()] += corr;
        correlation_squared[Threads.threadid()] += corr.^2;
    end

    ## sum all threads
    E = sum(EN, dims=1)[1];
    E2 = sum(EN_squared, dims=1)[1];
 
    ## calculate average and standard deviation
    # mean_energy = E./n_runs;
    # err_energy = sqrt.(E2./n_runs - mean_energy.^2);
 
    ## do the same for number of neighbours    
    Z = sum(Coord, dims=1)[1];
    Z2 = sum(Coord_squared, dims=1)[1];
 
    # mean_coord = Z./n_runs;
    # err_coord = sqrt.(Z2./n_runs - mean_coord.^2);
 
    ## do the same for number of accepted rotations
    accept = sum(accepted_theta, dims=1)[1];
    accept2 = sum(accepted_squared, dims=1)[1];
 
    # mean_accept = accept./n_runs;
    # err_accept = sqrt.(accept2./n_runs - mean_accept.^2);
 
    ## do the same for number of rejected rotations
    reject = sum(rejected_theta, dims=1)[1];
    reject2 = sum(rejected_squared, dims=1)[1];
 
    # mean_reject = reject./n_runs;
    # err_reject = sqrt.(reject2./n_runs - mean_reject.^2);
 
    ## do the same for correlations
    corr = sum(correlation, dims=1)[1];
    corr2 = sum(correlation_squared, dims=1)[1];
 
    # mean_corr = corr./n_runs;
    # err_corr = sqrt.(corr2./n_runs - mean_corr.^2);

    ## save as numpy arrays
    eps = trunc(Int, e)
	out_file = "run_e_$eps.npz"
# 	npzwrite(out_file, 
#    	Dict("a_range" => a_range, "rejected_theta" => mean_reject, 
# 	 	 "error_rejected" => err_reject, "accepted_theta" => mean_accept, 
#       	 "error_accepted" => err_accept, "E" => mean_energy, 
# 	  	 "error_energy" => err_energy, "Z" => mean_coord, 
# 	  	 "error_coord" => err_coord, "correlation" => mean_corr, 
# 	  	 "error_correlation" => err_corr, "grid" => grid)
#     );
 	npzwrite(out_file, 
    	Dict("a_range" => a_range, 
        "E" => E, "E2" => E2, "Z" => Z, "Z2" => Z2, "accept" => accept, "accept2" => accept2,
        "reject" => reject, "reject2" => reject2, "corr" => corr, "corr2" => corr2, 
        "grid" => grid)
    );
end

function backup_runs(a0, e, da, a_max, N, grid, n_runs, n_metropol, seed0)
    initial = hcat( [rand()*grid[1], rand()*grid[2]], [rand()*grid[1], rand()*grid[2]]);
    samples = generate_samples(initial, grid, N);

    ## change circles into ellipses
	b0 = a0/e;                          # small semi-axis

    ## set a_range for growing ellipses
    # a_range = [15.: 0.1: 30. ...]
    a_range = [a0: da:a_max ...];
    n_a = size(a_range, 1);

    # EN_results = Matrix{Float64}[]
    # Coord_results = Matrix{Float64}[]

    ## for multithreading
    EN_results = [Matrix{Float64}[] for _ in 1:Threads.nthreads()];
    Coord_results = [Matrix{Float64}[] for _ in 1:Threads.nthreads()];
    as = [Float64[] for _ in 1:Threads.nthreads()];

    accepted_theta, rejected_theta = zeros(n_a, n_runs), zeros(n_a, n_runs);
    correlation = zeros(100, n_runs);

    Threads.@threads for i in 1:n_runs
    # for i in 1:n_runs
        ## every run start from a distribution with same orientations
        Random.seed!(seed0);
		distribution = [Ellipse(convert(Vector{Float64}, vector), a0, b0) for vector in eachcol(samples)];

        ## deterministic random so that it can be recreated
        Random.seed!(i);

        (E, Z, accepted_theta[:,i], rejected_theta[:,i]) = one_run(distribution, a_range, e, grid, N, n_metropol);

        # push!(EN_results, E);
        # push!(Coord_results, Z);

        push!(EN_results[Threads.threadid()], E)
        push!(Coord_results[Threads.threadid()], Z)
        push!(as[Threads.threadid()], i)

        correlation[:,i] = angle_correlation(distribution, 3*a_max, grid)
    end

    EN = collect(flatten(EN_results));
    Coord = collect(flatten(Coord_results));

    # Coord = zip(flatten(as), flatten(Coord_results))

    ## do a bit of statistic
    all_energy = [0.5*sum(EN[i], dims=1) for i in 1:size(EN,1)];
    return_energy = mean(all_energy);
    err_energy = std(all_energy);

    all_coord = [sum(Coord[i], dims=1)/size(Coord[1][:,1], 1) for i in 1:size(Coord,1)];
    return_coord = mean(all_coord);
    err_coord = std(all_coord);

    norm_accepted = accepted_theta./(accepted_theta+rejected_theta)
    all_accepted = mean(norm_accepted, dims=2)';
    err_accept = std(norm_accepted, dims=2)';

    norm_rejected = rejected_theta./(accepted_theta+rejected_theta)
    all_rejected = mean(norm_rejected, dims=2)';
    err_reject = std(norm_rejected, dims=2)';

    all_correlation = mean(correlation, dims=2)';
    err_correlation = std(correlation, dims=2)';

    eps = trunc(Int, e)
	out_file = "run_e_$eps.npz"
	# npzwrite(out_file, 
 #    	Dict("a_range" => a_range, "rejected_theta" => all_rejected, 
	# 		"error_rejected" => err_reject, "accepted_theta" => all_accepted, 
 #     	    "error_accepted" => err_accept, "E" => return_energy, 
	#  		"error_energy" => err_energy, "Z" => return_coord, 
	#  		"error_coord" => err_coord, "correlation" => all_correlation, 
	#  		"error_correlation" => err_correlation, "grid" => grid)
 #    );
 #    
	npzwrite(out_file, 
     	Dict("a_range" => a_range, "accepted_theta" => norm_accepted, 
     	    "rejected_theta" => norm_rejected, 
     	    "E" => all_energy, "Z" => all_coord, 
	  		"correlation" => correlation, 
	  		"grid" => grid)
    );
     
 
    # println("#grid = $grid")
    # println("#a_range, energy, err_energy, coord, err_coord, rej, err_rej, acc, err_acc")
    # 
    # output1 = reduce(hcat, (a_range, return_energy', err_energy', return_coord', err_coord', all_rejected', err_reject', all_accepted', err_accept'))
    # Base.print_matrix(stdout, output1)
    # 
    # output2 = reduce(hcat, (all_correlation', err_correlation') )
    # Base.print_matrix(stdout, output2)
end


function create_gif(a0, amax, da, e, grid, N, n_metropol)
	initial = hcat( [rand()*grid[1], rand()*grid[2]], [rand()*grid[1], rand()*grid[2]]);
    samples = generate_samples(initial, grid, N);

    ## change circles into ellipses
    b0 = a0/e;                  # small semi-axis

    ## set a_range for growing ellipses
    a_range = [a0: da:amax ...];
    n_a = size(a_range, 1);
	
	data = Dict([]);
    Threads.@threads for i in 1:n_a
        ## change size of distributon
        a_i = a_range[i];
        b_i = a_i/e;

        distribution = [Ellipse(convert(Vector{Float64}, vector), a_i, b_i) for vector in eachcol(samples)];

        for j in 1:size(distribution,1)
            distribution[j].a = a_i;
            distribution[j].b = b_i;
            fix_A!(distribution[j]);
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
        E_curr, Z_curr, acc, rej = metropolis(distribution, in_proximity, grid, n_metropol, 0);
	
		## display state every iteration
		data["$i state"] = distribution;
    end
	data = sort(data, by = x -> parse(Int, (split(x, " ")[1])))
	JSON3.write("n_$n_metropol.json", data);
	return data
end

function main()
	pa = parse_commandline();
	println("Parsed args:");
	for (arg, val) in pa
		print("$arg => $val ");
		println(typeof(arg));
	end
	Random.seed!(pa["seed"]);
	grid = [pa["x_size"], pa["y_size"]];
	many_runs(pa["a_0"], pa["e"], pa["d_a"], pa["a_max"], pa["N"], grid, pa["n_runs"], pa["n_metropol"], pa["seed"]);
end

main()