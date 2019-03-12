# TMA4280 - Project 2
Repository for the second assignments in the course TMA4280 - Introduction to Super Computing

## User guide
All ansers to all questions (except plotting which is done in a seperate Python file for simplisity) are implementet in C++, and can be runned from src/main.cpp (see the `CMakeList.txt`).

main.cpp takes several command line arguments to be able to run the different tasks:
* -q - ("question") 
* -n - number of elements to be summed. not all tasks requieres this.

A typpical example will could be question 5 with method zeta2:

```
cd (path)/project1/build  
cmake ..  
make  
mpirun -np 4 project1/project1 -q 5 -m zeta2 -n 1000
```

