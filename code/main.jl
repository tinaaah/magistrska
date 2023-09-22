using DelimitedFiles 
using Distances
using LinearAlgebra 
using Optim
using Random
using Statistics

## A = R(theta) * A' * R^T(theta)
function Rotate(theta, matrix)
	R =	[cos(theta) -sin(theta);
		 sin(theta) cos(theta)];
	R_t = transpose(R);
	R*matrix*R_t
end

## generate distribution
function generate_samples(input_samples, grid, N, k0=0)
    samples = input_samples;
	
    while size(samples, 2) < N
		if k0==0
        	k = div(size(samples, 2), 2) + 1;
		else
			k = k0;
		end;
 
		## generate random k candidates
 		candidates = zeros(2, k);
 		for i in 1:k
			candidates[1, i] = rand() * grid[1]
			candidates[2, i] = rand() * grid[2]
 		end
 
		## calculate distance matrix between candidates and existing points
 		distances = pairwise(PeriodicEuclidean(grid), samples, candidates);
 
 		## select the one with the largest minimum distance
 		r_min = minimum(distances, dims=1);
 		i = argmax(vec(r_min));
 
 		samples = hcat(samples, candidates[:,i]);
	end
    samples
end

## define ellipses with a, b, center vector and generate random angle and A matrix
mutable struct Ellipse
    center::Vector{Float64};
    a::Float64;
    b::Float64;
	angle::Float64;
	A::Matrix{Float64};
	energy::Float64;
	coord::Int64;

    function Ellipse(center::Vector{Float64}, a::Float64, b::Float64)
        angle = rand() * 2*pi;
		A = Rotate(angle, diagm([a^2, b^2]));
		energy, coord = 0, 0;
        new(center, a, b, angle, A, energy, coord);
    end

end

## Fix ellipses after rotating/changing size
function fix_A!(ellipse::Ellipse)
    ellipse.A = Rotate(ellipse.angle, diagm([ellipse.a^2, ellipse.b^2]))
end

## calculate distance vector
function dist_vect(point1, point2, grid)
	grid2 = grid*0.5; 
	map(mod, point1 - point2 + grid2, grid) - grid2
end

## input is an ellipse: centre, width(a), height(b) and angle
## f = lambda*(1-lambda)*r^T * A^(-1) * r
function f(x, E1, E2, grid)
    A1 = E1.A;
    A2 = E2.A;

	dr = dist_vect(E1.center, E2.center, grid); 
    C = inv( (1-x)*A1 + x*A2 ); 

    x*(1-x)*dr' * C * dr
end

## mu is maximum value of f
function mu(E1, E2, grid)
	# x = vcat(range(0, 1, length=100));
	# maximum( map((vect) -> f(vect, E1, E2, grid), x) )
	objective(x) = -f(x, E1, E2, grid); 
	result = optimize(objective, 0.0, 1.0);
	max_value = -result.minimum
end

## energy is contact function between all the neighbours
function energy(ellipse_a, vicinity_a, grid)
	new_E, new_Z = 0, 0;
	if size(vicinity_a, 1) != 0
		all_mu = [mu(ellipse_a, ellipse_b, grid) for ellipse_b in vicinity_a];
		E = -log.(@. ifelse(all_mu < 1, all_mu, 1));
		new_Z = size(findall(x -> x!=0, E), 1);
		new_E = sum(E);
	end;
	return new_E, new_Z
end

## orientational correlation function
function psi_2(ellipse_a, a_vicinity)
	if size(a_vicinity, 1) < 2
		psi = 1;
	else
		psi = sum(exp.(2*im .* (ellipse_a.angle .- [ellipse_b.angle for ellipse_b in a_vicinity]))) / length(a_vicinity);
	end
	return psi
end

function angle_correlation(distribution, r_max, grid)
	radius_range = vcat(range(0, r_max, length=100))
    g = zeros(size(radius_range, 1));
	N = size(distribution, 1);

    for i in 1:size(radius_range,1)
        d_a = radius_range[i];

		## cget the vicinities of all ellipses
		S = [x.center for x in distribution]
		distance = pairwise(PeriodicEuclidean(grid), S)

		## matrix of neighbouring ellipses, which are less than 2*d_a apart
		in_proximity = findall(x -> (x<=2*d_a) & (x!=0), distance)

		## get vicinities of all ellipses
		results = [vicinity(i, in_proximity, distribution) for i in 1:N];
		(vicinities, indices) = (getindex.(results, 1), getindex.(results, 2))

        # Calculate all psi_2 at current r
        psi2 = [psi_2(distribution[j], vicinities[j]) for j in 1:N];

        # g(r) for each j
        g[i] = mean([abs(psi2[j] * conj(psi2[k])) for k in 1:N, j in 1:N]);
    end
    return g
end

## helper function
function vicinity(i, in_proximity, distribution)
	temp = [in_proximity[findall(x-> x[2] == i, in_proximity)][:][t][1] for t in 1:size(findall(x-> x[2] == i, in_proximity), 1) ];
	return distribution[temp], temp
	end

function metropolis(distribution, in_proximity, grid, n=100, T=0)
	N = size(distribution, 1);
	last_i = 0;

	t, E, Z = 0, zeros(n, N), zeros(n, N);
	accepted_theta, rejected_theta = Vector{Union{Float64,Missing}}(missing, n), Vector{Union{Float64,Missing}}(missing, n);

	## get vicinities of all ellipses
	results = [vicinity(i, in_proximity, distribution) for i in 1:N];
	(vicinities, indices) = (getindex.(results, 1), getindex.(results, 2))

	##energy of the starting system    
	result = [energy(distribution[i], vicinities[i], grid) for i in 1:N];
	(E[1,:], Z[1,:]) = (getindex.(result, 1), getindex.(result, 2))

	for i in 2:n 
		E[i,:] = E[i-1,:];
		Z[i,:] = Z[i-1,:];

		## select a random ellipse and copy it
		j = rand(1:N);
		j_ellipse = deepcopy(distribution[j]);
		

		## generate random angle thera
		j_ellipse.angle = rand_von_Mises(1, j_ellipse.angle, 3)[1];
		fix_A!(j_ellipse);

		## because delta theta can get to 2*pi
		function transform(x)
			x > pi && (x -= 2*pi);
			x < -pi && (x += 2*pi);
			x
		end
		delta_theta = transform(distribution[j].angle - j_ellipse.angle);

		## decide if you want to accept new step or not
		new_E, new_Z = energy(j_ellipse, vicinities[j], grid);
		delta_E = new_E - E[i,j];

		## boltzmann probability distribution
		u = rand();

		if (delta_E<=0 && u<1) | (delta_E>0 && delta_E<=T*log(1/u))
			## accept this step
			distribution[j].angle = j_ellipse.angle;
			fix_A!(distribution[j]);

			## save parameters for later
			## change energy
			E[i,j] = new_E;
			Z[i,j] = new_Z;

			## change energy of all neighbors as well
			for k in indices[j];
				E[i, k], Z[i, k] = energy(distribution[k], vicinities[k], grid);
			end

			## delta theta
			t += 1;
			accepted_theta[t] = delta_theta;
		else
			## save parameters for later
			rejected_theta[i-t] = delta_theta;
			continue
		end

		## If energy reaches zero stop iterating

		last_i = i
		if sum(E[i,:]) == 0
			break
		end
	end

	## calculate energy for all pairs after the final rotation
	result = [energy(distribution[i], vicinities[i], grid) for i in 1:N];
	(E[last_i,:], Z[last_i,:]) = (getindex.(result, 1), getindex.(result, 2))

	for k in 1:size(distribution, 1)
		distribution[k].energy = E[last_i,k];
		distribution[k].coord = Z[last_i,k];
	end

	return E, Z, collect(skipmissing(accepted_theta)), collect(skipmissing(rejected_theta));
end

## because delta theta can get to 2*pi
function transform(x)
	x > pi && (x -= 2*pi);
	x < -pi && (x += 2*pi);
	x
end

## this algorithm returns only the last state
function metropolis2(distribution, in_proximity, grid, n=100, T=0)
	N = size(distribution, 1);
	last_i = 0;

	t=0;
	accepted, rejected = 0, 0;
	# accepted_theta, rejected_theta = Vector{Union{Float64,Missing}}(missing, n), Vector{Union{Float64,Missing}}(missing, n);

	## get vicinities of all ellipses
	results = [vicinity(i, in_proximity, distribution) for i in 1:N];
	(vicinities, indices) = (getindex.(results, 1), getindex.(results, 2))

	##energy of the starting system    
	result = [energy(distribution[i], vicinities[i], grid) for i in 1:N];
	(E, Z) = (getindex.(result, 1), getindex.(result, 2));

	for i in 2:n 
		## select a random ellipse and copy it
		j = rand(1:N);
		j_ellipse = deepcopy(distribution[j]);
		
		## generate random angle thera
		j_ellipse.angle = rand_von_Mises(1, j_ellipse.angle, 3)[1];
		fix_A!(j_ellipse);

		delta_theta = transform(distribution[j].angle - j_ellipse.angle);

		## decide if you want to accept new step or not
		new_E, new_Z = energy(j_ellipse, vicinities[j], grid);
		delta_E = new_E - E[j];

		## boltzmann probability distribution
		u = rand();

		if (delta_E<=0 && u<1) | (delta_E>0 && delta_E<=T*log(1/u))
			## accept this step
			distribution[j].angle = j_ellipse.angle;
			fix_A!(distribution[j]);

			## save parameters for later
			## change energy
			E[j] = new_E;
			Z[j] = new_Z;

			## change energy of all neighbors as well
			for k in indices[j];
				E[k], Z[k] = energy(distribution[k], vicinities[k], grid);
			end

			## delta theta
			t += 1;
			accepted += 1
			# accepted_theta[t] = delta_theta;
		else
			## save parameters for later
			# rejected_theta[i-t] = delta_theta;
			rejected += 1
			continue
		end

		## If energy reaches zero stop iterating
		last_i = i
		if sum(E) == 0
			break
		end
	end

	## calculate energy for all pairs after the final rotation
	result .= [energy(distribution[i], vicinities[i], grid) for i in 1:N];
	(E, Z) = (getindex.(result, 1), getindex.(result, 2))

	for k in 1:size(distribution, 1)
		distribution[k].energy = E[k];
		distribution[k].coord = Z[k];
	end

	# return E, Z, collect(skipmissing(accepted_theta)), collect(skipmissing(rejected_theta));
	return sum(E, dims=1)[1]*0.5, sum(Z, dims=1)[1]/N, accepted, rejected;
end
