# Rubik's Cube Solver (Parallel) with Multithreading and GPU

## Introduction
The purpose of this program is to solve Rubik's Cube and to explore the power of parallel computing.

The main programming language is C++ 17, CUDA is also used for Nvidia GPU support.

## Algorithm
There are three algorithms: two phase, optimum X, optimum Y.

Two phase algorithm requires 1GB RAM. It runs very fast but the solution might not optimum.

Optimum X algorithm requires 2GB RAM. It always outputs optimum solution.

Optimum Y algorithm requires 4GB RAM. It always outputs optimum solution, and runs faster than optimum X.

Optimum algorithm has three implements: single thread, multi thread, CUDA.

The CUDA implement requires the same size of GPU RAM as CPU RAM.

## Build

### Linux
To build the program, Clang 10 or higher with C++ 17 support is required.

To use the power of parallel, a Nvidia GPU and CUDA toolkit 10 is required.

GCC is not recommend because it is much slower and exhausts much more RAM.
If you need work with GCC, modify the script by yourself.

```shell
bash ./build_linux.sh
```

### MacOS
To build the program, Clang 10 or higher with C++ 17 support is required.
The build-in Clang of the newest macOS is OK.

```shell
bash ./build_macos.sh
```

### Windows
I have not tested on Windows. It must work with a C++ 17 compiler.

## Cache files
The program will create some cache files at `cache/` dir when needed.
Next time it will load them to save time.

The SHA 256 checksum of all cache files is in `sha256.txt`.

## Run test

Show help.
```shell
./test --help
```

It takes about 20 minutes to run basic test.
```shell
./test  --n_thread 4  --seed 0
```

It takes about 1 hour to run full test.
```shell
./test  --n_thread 4  --seed 0  --full
```

## Run solver

Show help.
```shell
./solver --help
```

### Example

Two Phase
```shell
./solver  --algorithm 2p  --n_thread 4  --2p_n_moves 24  --input example.txt  --output result.txt
```

Optimum Y
```shell
./solver  --algorithm opty  --n_thread 4  --sym_n_moves 6  --n_solution=1  --input example.txt  --output result.txt
```

Thread Optimum Y
```shell
./solver  --algorithm thread_opty  --schedule simple  --n_thread 4  --bfs_count=100000  --input example.txt  --output result.txt
```

CUDA Optimum Y
```shell
./solver  --algorithm cuda_opty  --schedule simple  --n_thread 4  --n_cuda_thread 4096  --bfs_count=100000  --input example.txt  --output result.txt
```

Check and modify `example.txt` to solve your own. 

### Example Superflip

Superflip is the first state found witch need at least 20 steps to solve.

We can use symmetry to reduce the search space. 

Thread Optimum Y
```shell
./solver  --algorithm thread_opty  --schedule simple  --n_thread 4  --bfs_count=100000  --input example_superflip.txt  --output result.txt
```

CUDA Optimum Y
```shell
./solver  --algorithm cuda_opty  --schedule simple  --n_thread 4  --n_cuda_thread 4096  --bfs_count=100000  --input example_superflip.txt  --output result.txt
```

Modify the options and run the program again.

### Example 20

Those are some most difficult sates witch need at least 20 steps to solve.

Thread Optimum Y
```shell
./solver  --algorithm thread_opty  --schedule simple  --n_thread 4  --bfs_count=100000  --input example_20.txt  --output result.txt
```

CUDA Optimum Y
```shell
./solver  --algorithm cuda_opty  --schedule simple  --n_thread 4  --n_cuda_thread 4096  --bfs_count=100000  --input example_20.txt  --output result.txt
```

Modify the options and run the program again.

## Benchmark

Hardware: Intel Xeon E5 2690 and Nvidia RTX 2080 Ti

Time to solve superflip with optimum Y algorithm.

implement               | time (second) | speedup
----------------------- | ------------- | -------
CPU 1 thread (baseline) |       140.723 |   1.000
CPU 2 threads           |        73.000 |   1.928
CPU 4 threads           |        37.168 |   3.786
CPU 8 threads           |        19.653 |   7.160
GPU 4096 threads        |         6.918 |  20.342

## Reference

http://kociemba.org/cube.htm

https://www.cube20.org/
