// 01300431
// jg719@ic.ac.uk
// Compiler used: Apple clang version 11.0.0 (clang-1100.0.33.17)
// OS: macOS Catalina
// Number of logical processors: hw.physicalcpu: 2; hw.logicalcpu: 4
// Serial execution time (100x100x100): 0.299904s
// Parallel execution time (100x100x100) (NOT VECTORISED): 0.0595788s (compiled with "clang++ -Xpreprocessor -fopenmp -lomp ConwaysGame_Parallel.cpp")
// Parallel execution time (100x100x100) (VECTORISED):0.0104483s (compiled with "clang++ -Xpreprocessor -fopenmp -lomp -O3 ConwaysGame_Parallel.cpp")

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>

#include <omp.h>
#include <chrono>  

//Note that this is a parallel implementation of the Game of Life with a periodic grid
// Parallelism has been achieved through OpenMP, resulting in a speed-up of up to a factor of 5 wrt serial
std::vector<std::vector<unsigned short int> > grid, new_grid; // These will store 0's and 1's. Unsigned short int is the datatype that occupies the least space in memory 
// This results in a huge speed-up (it takes half of the time to execute)
int imax, jmax; // If the grid is squared, then jmax can be replaced by imax, saving 1 variable

// Include function prototypes
int num_neighbours(int ii, int jj); 
void grid_to_file(int it);
void do_iteration(void);


int main(){
	srand(time(NULL)); 
	int i, j, n, max_steps = 100; 
	double start, end, ave = 0;
	imax = 100;
	jmax = imax;
	for(int u = 0; u<10; u++){
		start = omp_get_wtime(); // Start time measurement
		srand(time(NULL)); // Seed pesudo-random number generator
		grid.resize(imax, std::vector<unsigned short int>(jmax)); // Shape domain as matrix of imax x jmax (rows x cols) 
		new_grid.resize(imax, std::vector<unsigned short int>(jmax)); 

		// The possibility of storing the matrices in 1D has been considered. However, this has proved to not lead to any improvements in efficieny
		omp_set_num_threads(4); // Number of threads to run on 
		// For loops below are used to initialise the grid 
		#pragma omp parallel for collapse(2) private(i, j)  // Make a copy of i and j into the thread's stack for faster access
		// Collapsing these for loops results in a very slight increase in efficiency (around 2%)
		for (i = 0; i < imax; i++){
			for (j = 0; j < jmax; j++){
				grid[i][j] = (rand() % 2);
			}
		}

		// For loop below performs iterations
		for (n = 0; n < max_steps; n++){ // Should not be parallelised (the value of the next iteration depends on the current one)
			do_iteration(); 
			// grid_to_file(n); // Output is omitted for speed-up (uncomment this for debugging)
		}

		end = omp_get_wtime(); // End time measurement
		ave += end - start; // Store time taken 
		}
		std::cout << " Average time taken: " << ave/10 << " s" << std::endl; 
	return 0;
}

int num_neighbours(int ii, int jj){
	int ix, jx, i_ii, i, j;
	int cnt = 0;
	// The for loops below find the neighbour count of the current cell 
	for (i = -1; i <= 1; i++){ // Attempting to parallelise these for loops will result in a significant decrease in speed
	// These for loops result in 9 iterations which are performed fastly. The overhead introduced by having to manage threads is greater than the speed-up they create. 
		i_ii = i + ii; // This is used to save 3 addition operations in the next for loop (avoids repeating (i+ii))
		for (j = -1; j <= 1; j++){
			if (i != 0 || j != 0){
				ix = (i_ii + imax) % imax;
				jx = (j + jj + jmax) % jmax;
				cnt+=grid[ix][jx]; // This lines used to be "if (grid[ix][jx]) cnt++;". 
				// Instead, avoid evaluating the if statement 9 times per cell and iteration and simply sum the values of the neighbours
				// This is a further benefit from storing unsigned short ints rather than booleans
			}
		}
	}
	return cnt;
}

void grid_to_file(int it){ // Used for ensuring correctness of output - Nothing to be optimised here
	std::ofstream fname;
	std::fstream f1;
	std::string s0( std::to_string(it) + "_out.txt" );
	f1.open(s0, std::ios::out);
	for (int i = 0; i < imax; i++){
		for (int j = 0; j < jmax; j++)
			f1 << grid[i][j] << "\t";
		f1 << std::endl;
	}
	f1.close();
}

void do_iteration(void){
	int i, j, num_n;

	// The for loop below applies the rules of the game of life
	#pragma omp parallel for schedule(auto) private(i, j, num_n) collapse(1) // Static, dynamic, guided and automatic scheduling has been considered. 
	// Automatic scheduling has shown to result in the highest speed-ups. A copy of "i", "j" and "num_n" is stored in the thread's stack for faster access
	// The possibility of collapsing the nested for loops has been considered, but this does not result in any further speed-up (just left as collapse(1))
	for (i = 0; i < imax; i++)
		for (j = 0; j < jmax; j++){ 
			new_grid[i][j] = grid[i][j];
			num_n = num_neighbours(i, j);
			if (grid[i][j]){
				if (num_n != 2 && num_n != 3)
					new_grid[i][j] = 0;
			}
			else if (num_n == 3) new_grid[i][j] = 1;
		}
	grid.swap(new_grid); // More efficient than std's swap
}
