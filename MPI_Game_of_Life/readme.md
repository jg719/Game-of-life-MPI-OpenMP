OVERVIEW

The purpose of this assignment is to implement a scalable version of the famous cellular automaton Game of Life model (https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life). The serial solution of this problem results in exponential increases in computing time as the problem size is changed. This problem is overcome through this parallel implementation. 

Domain decomposition has been used to split the board into blocks which are as even as possible. Each process is in charge of solving the game on their block, and must, therefore, communicate its boundaries to its neighbour processes. This raises the need for message passing between processes, which is satisfied by using Open MPI. 

The code has been executed on Imperial’s CX1 HPC for high-dimensional grids, using up to 32 cores in parallel. This results in an increase in speed-up and yields a reasonable parallel efficiency for up to 16 cores. 

Lastly, the correctness of the code has been ensured through testing. Simple gliders have been generated in such a way that the middle, boundary and corner cases are tested – both for periodic and non-periodic boundaries. 

![alt text](https://github.com/jg719/Game-of-life-MPI-OpenMP/blob/master/MPI_Game_of_Life/post-processing%20and%20animations/animations/periodic_big.gif "Game of Life local simulation on 4 processes")

REPOSITORY STRUCTURE

Source 

	Parallel_Game_of_life.cpp – Source code containing MPI implementation of the game of life. 
	
	a.out – Compiled executable version of the code
	
HPC

	Makefile: Used for building  
	
	my_script.pbs: PBS file used to specify settings when submitting a job to the HPC
	
	Output files: Contains the files returned by the HPC displaying timing information. 
	
Performance 

	Report.pdf: Discussion of the code’s performance 
	
	Data.xlsx: Timing data of the code’s performance
	
Post-processing and animations 

	Results: Contains the files saved for a test run locally for a 100x100 domain in 4 processors. 
	
	Animations: Contains different animations of the game. 
	
	Post_processing.ipynb: Python script used for creating a visualisation of the output files. 


KEY CODE FEATURES

	MPI datatypes

Using such a data structure bypasses the need of making a contiguous copy of non-contiguously stored elements by pointing at them. 
Even though the columns are the only non-contiguous data to be communicated, datatypes have been created for the rest of the elements for consistency. This results in 16 datatypes – corresponding to the send and receive of rows, columns and corners. 

	Adaptive domain decomposition

To assign the dimensions of a subdomain to each process, a brute force algorithm is employed. For a given number of processes, this algorithm returns a rectangular distribution which aims to distribute the processes as evenly as possible (i.e.: 9 processes are arranged in a 3x3 matrix). This approach is perfect for square grids but might be suboptimal for highly uneven domains. Further work would involve adapting this algorithm such that it arranges the processes by looking at a ratio between the domain’s rows and columns. 

	Periodicity

The boundaries of the domain can be switched between periodic and non-periodic by the user. Non-periodic boundaries are implemented by having a boundary process send data to itself – which in Open MPI results in no communication at all and therefore saves the need for unnecessary communications. 

	Post-processing
	
Functions for outputting each processor’s data and information have been implemented. These are fed into a Jupyter notebook that resembles the original domain and creates an animation for the user to visualise how the game evolves. 

	Load balancing analysis

Reduction operations have been implemented with the objective of finding the difference between the processors that take the longest and the shortest periods to execute. Ideally (for a perfectly balanced domain), this difference would be close to zero. However, this is not true for highly uneven domains that result in suboptimal domain decompositions – this reinforces the idea of using an improved domain decomposition algorithm. 


GAME INITIALISATION 

The initial state of the domain can be set in three different ways:
-	Random: A pseudo-random number generator is used to produce a random distribution of dead and alive cells in the domain. 

-	Gliders: A simple function that that places a glider on the top left corner of the domain has been created for testing purposes. 

-	Hardcoded: Patterns can be easily inserted into the code by hand.

Further work would involve allowing the user to input an initial state through a file. 

The source code has been locally compiled using the following command:

	"mpicxx -std=c++17 -O3 Parallel_Game_of_life.cpp"

Note that the vectorisation flag is not used when performing tests on the HPC. 


POST-PROCESSING

Each processor’s data and information is output into the results directory in a .txt file. These files are ready for being processed by the Python script included in the form of a Jupyter notebook as part of this repository. 
The Jupyter notebook creates a GIF animation to enable the user to visualise the iterations of the game of life. This is done by creating an n-D array that contains a set of arrays corresponding to each iteration. An example of these animations and post-processing can be found in the animations directory.  
