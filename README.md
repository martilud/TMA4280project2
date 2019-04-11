# TMA4280 - Project 2
Repository for the second assignments in the course TMA4280 - Introduction to Super Computing

## User guide
The project is implemented in C++, and the main source file is stored in the `src/parallel-poisson.cpp`. 

The program takes several command line parameters:
* `-n`: to set the problem size (must be a power of 2)
* `-p`: number of OpenMP threads to use
* `-t`: spesify if any test should be used:
  * `u` for unit test
  * `v` for verification test
* `-f`: to spesify which function to be used in the poisson problem

In order to run the program one can e.g. run the following commands:
```
cd (path)/project2/build  
cmake ..  
make  
mpirun -np 4 ./parallelPoisson -n 128 -p 4 -t v
```

To get more information on how to use the program on a super computer
https://www.hpc.ntnu.no/display/hpc/Running+Jobs
