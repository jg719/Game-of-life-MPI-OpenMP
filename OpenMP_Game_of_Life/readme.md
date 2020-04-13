CID: 01300431
E-mail: jg719@ic.ac.uk
Compiler used: Apple clang version 11.0.0 (clang-1100.0.33.17)
OS: macOS Catalina
Number of logical processors: hw.physicalcpu: 2; hw.logicalcpu: 4
Serial execution time (100x100x100): 0.299904s
Parallel execution time (100x100x100) (NOT VECTORISED): 0.0595788s (compiled with "clang++ -Xpreprocessor -fopenmp -lomp ConwaysGame_Parallel.cpp")
Parallel execution time (100x100x100) (VECTORISED):0.0104483s (compiled with "clang++ -Xpreprocessor -fopenmp -lomp -O3 ConwaysGame_Parallel.cpp")


This repository contains a parallel version of Conway's Game of Life using OpenMP. Further changes have been carried out with respect to the the serial version to improve efficiency. To find more details, please refer to the "ConwaysGame_Parallel.cpp" documentation. 

A performance analysis is also included "Performance1.png", "Performance2.png" and "Performance_anaysis.xlsx". This files show a comparison between the performance of the Serial, OpenMP and OpenMP + vectorisation versions. 
