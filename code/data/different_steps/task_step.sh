#!/bin/bash

steps=( 0.1 0.5 0.01 0.005 )
for step in "${steps[@]}" ; do

	time julia --threads auto ../../run_julia.jl --n_0 0.15 --n_max 1.15 --d_n "$step" --e 2.5 --N 1024 --x_size 1000 --y_size 1000 --n_runs 8 --n_metropol 10000 --seed 42

done
