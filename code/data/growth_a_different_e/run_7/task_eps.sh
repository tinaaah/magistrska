#!/bin/bash

eps=( 1.5 2.5 4 5 )
for e in "${eps[@]}" ; do

    time julia --threads auto ../../run_julia.jl --a_0 15. --a_max 30. --d_a 0.1 --e "$e" --N 1024 --x_size 1000 --y_size 1000 --n_runs 25 --n_metropol 10000 --seed 42

done
