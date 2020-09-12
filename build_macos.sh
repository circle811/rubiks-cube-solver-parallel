#!/bin/bash

time clang++ -std=c++17 -fconstexpr-steps=100000000 -I include/ -O2 src/test.cpp -o test

time clang++ -std=c++17 -fconstexpr-steps=100000000 -I include/ -O2 -c src/solver.cpp -o solver.o
time clang++ -std=c++17 -fconstexpr-steps=100000000 -I include/ -O2 -c src/cuda_cube.cpp -o cuda_cube_nocuda.o
time clang++ -std=c++17 -fconstexpr-steps=100000000 -I include/ -O2 solver.o cuda_cube_nocuda.o -o solver_nocuda

ln -s solver_nocuda solver
