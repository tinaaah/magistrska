using JSON3
using LinearAlgebra
using NPZ
using Plots
using Random

include("main.jl")
include("von_mises.jl")

 Random.seed!(42)

## generate distribution
grid = [100, 100];
N = 128; 

## select two random points on a grid
initial = hcat( [rand()*grid[1], rand()*grid[2]], [rand()*grid[1], rand()*grid[2]]);

## distribution of N points
samples = generate_samples(initial, grid, N);

## uncomment if you want plot
#=
    x = samples[1,:]
    y = samples[2,:]

    plotek = plot(x, y, seriestype = :scatter, title = "datapoints", legend=false)
    savefig("test.png")
=#
## turn points into ellipses with random orientations
a, e = 15., 5/2       ## big semi-axis and ratio
b = a/e             ## small semi-axis

distribution = [Ellipse(convert(Vector{Float64}, vector), a, b) for vector in eachcol(samples)]

## let's test the contact function
S = [x.center for x in distribution]
distance = pairwise(PeriodicEuclidean(grid), S)

## matrix of neighbouring ellipses, which are less than 2*a apart
in_proximity = findall(x -> (x<=2*a) & (x!=0), distance)

## get vicinities of all ellipses
results = [vicinity(i, in_proximity, distribution) for i in 1:N];
(vicinities, indices) = (getindex.(results, 1), getindex.(results, 2))

    ## check the viccinity of a random ellipse
    i = rand(1:N)

    i_ellipse = distribution[i]
    i_ellipse.colour = "#f94144"

    temp = [in_proximity[findall(x-> x[2] == i, in_proximity)][:][t][1] for t in 1:size(findall(x-> x[2] == i, in_proximity), 1) ]
    i_neighbourhood = distribution[temp]
    for j_ellipse in i_neighbourhood
        j_ellipse.colour = "#f8961e"
    end

    overlap = [mu(i_ellipse, neighbour, grid) for neighbour in i_neighbourhood]
    print(overlap)

    JSON3.write("output.json", distribution)

#= 
## try to do the metropolis

    data = Dict("distribution" => distribution, "grid" => grid);
    JSON3.write("output.json", data)

    ## correlation from 0 to 3 times the large axis
    correlation = angle_correlation(distribution, 3*a)

    npzwrite("data.npz", Dict("rejected_theta" => rejected_theta, "accepted_theta" => accepted_theta, "E" => E, "Z" => Z, "correlation" => correlation, "grid" => grid))
=# 

## enlarge and rotate ellipses
