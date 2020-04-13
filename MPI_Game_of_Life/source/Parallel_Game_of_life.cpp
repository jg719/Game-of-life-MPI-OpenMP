#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <chrono>  
#include <cstdlib> 

using namespace std;

int id, p, tag_num = 1;
int rows, cols; // rows and columns of processes (for example, 6 proceses can be arranged as 2 rows by 3 cols)
int ROWS, COLS; // rows and columns of problem's domain
int id_row, id_col;
int local_rows, local_cols;  // rows and columns of subdomain
int iterations;
bool periodic = true; // change to switch off periodic boundaries

class mpi_com_class{
    // The purpose of this class is to create two vectors: one contains 8 send datatypes, the other another 8 datatypes for receiving
    // It gets the rows and cols of the subdomain, and get the subdomain values (all these are passed by reference and assigned through the constructor)
    public:
        int local_rows;
        int local_cols;
        bool **domain;
        vector<MPI_Datatype> typelist_send_middle; // vector to be filled with 8 send datatypes
        vector<MPI_Datatype> typelist_receive_middle; // vector to be filled with 8 receive datatypes
        // Constructor
        mpi_com_class(int &local_rows, int &local_cols, bool ** &domain){ //Changed this so that they are passed by reference
            this->local_rows = local_rows;
            this->local_cols = local_cols;
            this->domain = new bool*[this->local_rows+2];
            for(int i = 0; i<local_rows+2; i++){
                this->domain[i] = new bool[local_cols+2];
                for(int j = 0; j<local_cols+2; j++){
                    this->domain[i][j] = domain[i][j];
                }
            }
        }
        // Destrutor
        ~mpi_com_class(){
            for(int j = 0; j<this->local_rows+2; j++){
                delete[] this->domain[j];
            }
            delete[] this->domain;
        }
        MPI_Datatype mpi_send_left, mpi_send_right, mpi_send_top,  mpi_send_bottom, mpi_send_top_left_corner, mpi_send_top_right_corner, mpi_send_bottom_left_corner, mpi_send_bottom_right_corner; 
        MPI_Datatype  mpi_receive_right,mpi_receive_left, mpi_receive_top,  mpi_receive_bottom, mpi_receive_top_left_corner, mpi_receive_top_right_corner, mpi_receive_bottom_left_corner, mpi_receive_bottom_right_corner; 
        void mpi_sendtype();
        void mpi_receivetype();
};

void mpi_com_class::mpi_sendtype(){
    // Creation of 8 datatypes for sending data (left hand col, right hand col, bottom row, top row, etc)
    vector<int> block_length;
    vector<MPI_Aint> addresses;
    vector<MPI_Datatype> typelist;
    //Left hand send
    for (int i = 1; i <=this->local_rows; i++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[i][1], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_left);
    MPI_Type_commit(&this->mpi_send_left);
    this->typelist_send_middle.push_back(this->mpi_send_left);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //Right hand send
    for (int i = 1; i <=this->local_rows; i++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[i][local_cols] , &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_right);
    MPI_Type_commit(&this->mpi_send_right);
    this->typelist_send_middle.push_back(this->mpi_send_right);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //Top send
    for (int j = 1; j <=this->local_cols; j++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[1][j], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }

    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_top);
    MPI_Type_commit(&this->mpi_send_top);
    this->typelist_send_middle.push_back(this->mpi_send_top);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //bottom send
    for (int j = 1; j <=this->local_cols; j++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[local_rows][j], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_bottom);
    MPI_Type_commit(&this->mpi_send_bottom);
    this->typelist_send_middle.push_back(this->mpi_send_bottom);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //Top left corner send
    block_length.push_back(1);
    {
    MPI_Aint temp;
    MPI_Get_address(&domain[1][1], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_top_left_corner);
    MPI_Type_commit(&this->mpi_send_top_left_corner);
    this->typelist_send_middle.push_back(this->mpi_send_top_left_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //Top right corner send
    block_length.push_back(1);
    {
    MPI_Aint temp;
    MPI_Get_address(&domain[1][local_cols], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_top_right_corner);
    MPI_Type_commit(&this->mpi_send_top_right_corner);
    this->typelist_send_middle.push_back(this->mpi_send_top_right_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //Bottom left corner send
    block_length.push_back(1);
     {
    MPI_Aint temp;
    MPI_Get_address(&domain[local_rows][1], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_bottom_left_corner);
    MPI_Type_commit(&this->mpi_send_bottom_left_corner);
    this->typelist_send_middle.push_back(this->mpi_send_bottom_left_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    /************************************/
    //Right bottom corner send
    block_length.push_back(1);
        {
    MPI_Aint temp;
    MPI_Get_address(&domain[local_rows][local_cols], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_send_bottom_right_corner);
    MPI_Type_commit(&this->mpi_send_bottom_right_corner);
    this->typelist_send_middle.push_back(this->mpi_send_bottom_right_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();
}

void mpi_com_class::mpi_receivetype(){
    // Creation of 8 datatypes for receiving data (left hand col, right hand col, bottom row, top row, etc)
    // The order in which these datatypes are pushed into the vector is the opposite to the order
    // of the datatypes in the send vector - so that they match
    // For example, Left col is matched with right col; top left corner is matched with bottom right corner, and so on
    vector<int> block_length;
    vector<MPI_Aint> addresses;
    vector<MPI_Datatype> typelist;
    //Right hand receive
    for (int i = 1; i <=this->local_rows; i++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[i][local_cols+1], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_right);
    MPI_Type_commit(&this->mpi_receive_right);
    this->typelist_receive_middle.push_back(this->mpi_receive_right);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    //Left hand receive
    for (int i = 1; i <=this->local_rows; i++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[i][0], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_left);
    MPI_Type_commit(&this->mpi_receive_left);
    this->typelist_receive_middle.push_back(this->mpi_receive_left);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    //Bottom receive
    for (int j = 1; j <=this->local_cols; j++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[local_rows+1][j], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_bottom);
    MPI_Type_commit(&this->mpi_receive_bottom);
    this->typelist_receive_middle.push_back(this->mpi_receive_bottom);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    //Top receive
    for (int j = 1; j <=this->local_cols; j++){
        block_length.push_back(1);
        MPI_Aint temp;
        MPI_Get_address(&domain[0][j], &temp);
        addresses.push_back(temp);
        typelist.push_back(MPI_C_BOOL);
    }
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_top);
    MPI_Type_commit(&this->mpi_receive_top);
    this->typelist_receive_middle.push_back(this->mpi_receive_top);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    MPI_Aint temp; // check this line
    //Right bottom corner receive
    block_length.push_back(1);
    {
    MPI_Aint temp;
    MPI_Get_address(&domain[local_rows+1][local_cols+1], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_bottom_right_corner);
    MPI_Type_commit(&this->mpi_receive_bottom_right_corner);
    this->typelist_receive_middle.push_back(this->mpi_receive_bottom_right_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    //Bottom left corner receive
    block_length.push_back(1);
    {
    MPI_Aint temp;
    MPI_Get_address(&domain[local_rows+1][0], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);

    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_bottom_left_corner);
    MPI_Type_commit(&this->mpi_receive_bottom_left_corner);
    this->typelist_receive_middle.push_back(this->mpi_receive_bottom_left_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();


    //Top right corner receive
    block_length.push_back(1);
    {
    MPI_Aint temp;
    MPI_Get_address(&domain[0][local_cols+1], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_top_right_corner);
    MPI_Type_commit(&this->mpi_receive_top_right_corner);
    this->typelist_receive_middle.push_back(this->mpi_receive_top_right_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();

    
    //Top left corner receive
    block_length.push_back(1);
    {
    MPI_Aint temp;
    MPI_Get_address(&domain[0][0], &temp);
    addresses.push_back(temp);
    } 
    typelist.push_back(MPI_C_BOOL);
    //Either do an offset and use pointer to domain in the send or not use an offset (done here) and use MPI_BOTTOM as the pointer
    MPI_Type_create_struct(block_length.size(), &block_length[0], &addresses[0], &typelist[0], &this->mpi_receive_top_left_corner);
    MPI_Type_commit(&this->mpi_receive_top_left_corner);
    this->typelist_receive_middle.push_back(this->mpi_receive_top_left_corner);
    block_length.clear();
    addresses.clear();
    typelist.clear();
}

int random_bool(){ // Slight modification on usual rand() - To ensure more randomness 
    double randomval = (rand()%300 +1);
    if (randomval<150){
        randomval = 0;
    }
    else{
        randomval = 1;
    }
   return randomval;
}

bool** random_subdomain(int local_rows, int local_cols){ // Initialises a double pointer array with random values (0's or 1's)
    bool **domain;
    domain = new bool*[local_rows+2];
    for(int i = 0; i<local_rows+2; i++){
        domain[i] = new bool[local_cols+2];
        for(int j = 0; j<local_cols+2; j++){
            domain[i][j] = false;
        }
    }
    for(int i = 1; i<local_rows+1; i++){
        for(int j = 1; j<local_cols+1; j++){
            domain[i][j] = random_bool();
        }
    }
    return domain;
}

void id_to_index(int id, int &id_row, int &id_col){ 
	// Wrapper to convert from id to index
	id_col = int(id % cols);
	id_row = int(id / cols);
}

int id_from_index(int id_row, int id_column, int rows, int cols){
    // Wrapper to find id from index
	return id_row*cols + id_column;
}

void find_dimensions(int p, int &rows, int &cols){		//A bit brute force to split the domain into processors
	int min_gap = p;
	// Iterate between 1 and p/2 (least optimal array size is p/2 x 2)
	for (int i = 1; i <= p/2; i++){  
		// If the total number of processes is divisible by i
		if (p%i == 0){
			// Create gap (+++)
			int gap = abs(p / i - i);
			if (gap < min_gap){
				min_gap = gap;	// update gap
				rows = i;	// make rows
				cols = p / i; // make cols 
			}
		}
	}
}

int num_neighbours(int ii, int jj, int local_rows, int local_cols, bool ** &domain){
	int ix, jx;
	int cnt = 0;
    // Find the number of alive neighbours a cell has
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0){
				ix = (i + ii + local_rows) % local_rows;
				jx = (j + jj + local_cols) % local_cols;
				if (domain[ix][jx]) cnt++;
			}
	return cnt;
}

void swap_arrays(bool ** &domain, bool ** &new_domain, int local_rows, int local_cols){
    // Swaps the two input arrays (domain and new_domain)
    bool **temp;
    temp = new bool*[local_rows+2];
    for(int i=0; i<local_rows+2; i++){
        temp[i] = new bool[local_cols+2];
        for(int j = 0; j<local_cols+2; j++){
            temp[i][j] = domain[i][j];
            domain[i][j] = new_domain[i][j];
            new_domain[i][j] = temp[i][j];
        }
    }
    for(int j = 0; j<local_rows+2; j++){
        delete[] temp[j];
    }
    delete [] temp;
}

void do_iteration(bool ** &domain, bool ** &new_domain, int local_rows, int local_cols){
// Appplies rules of Conway's game of life
	for (int i = 0; i < local_rows+2; i++)
		for (int j = 0; j < local_cols+2; j++){
			new_domain[i][j] = domain[i][j];
			int num_n = num_neighbours(i, j, local_rows+2, local_cols+2, domain); 
			if (domain[i][j]){
				if (num_n != 2 && num_n != 3)
					new_domain[i][j] = false;
			}
			else if (num_n == 3) new_domain[i][j] = true;
		}
    swap_arrays(domain, new_domain, local_rows, local_cols);
}

void copy_domain(bool ** &domain, bool ** &new_domain, int local_rows, int local_cols){
    // After each iteration, this function is in charge of copying the domain into the new_domain
    for(int i=0; i<local_rows; i++){
        for(int j = 0; j<local_cols; j++){
            new_domain[i][j] = domain[i][j];
        }
    }
}

void print_array(bool ** &array, int local_rows, int local_cols, int id){
    cout << "Printing processor " << id << " array: " << endl;
    for(int i = 0; i< local_rows+2; i++){
        for(int j = 0; j<local_cols+2; j++){
            cout << " " << array[i][j] << " ";
        }
        cout << endl;
    }
}

void initialise_array_int(bool ** &domain, int local_rows, int local_cols, int number){
    // Takes a doube pointer and initialises an array with 0's
    // initialises empty array
        domain = new bool*[local_rows+2];
        for(int i = 0; i<local_rows+2; i++){
            domain[i] = new bool[local_cols+2];
            for(int j = 0; j<local_cols+2; j++){
                domain[i][j] = false;
            }
        }
}

void create_glider(bool ** &domain, int local_rows, int local_cols, int number){
    // Creates a simple glider used for simple testing purposes 
    // If to be used, careful with its position
    domain[7][12] = true;
    domain[7][13] = true;
    domain[7][14] = true;
    domain[5][13] = true;
    domain[6][14] = true;
}

void i_j_to_index(int &i, int &j, int &index){
    // Used for mapping the i, j indeces into the correct datatype
    if(i+j == 0){
        if(i == 1){ // TL
            index = 5;
        }
        else if(i==-1){ // => i = -1  => BR 
            index = 6;
        }
    }
    else if(i + j  == 1){
        if(i == 1){ // T
            index = 2;
        }
        else if(i == 0){ // => i =0 => R
            index = 0;
        }
    }
    else if(i + j == -1){
        if(i == -1){ //B
            index = 3;
        }
        else{   // // L
            index = 1;
        }
    }
    else if(i + j == 2){ // TR
        index = 4;
    }
    else if(i + j == -2){ // BL
        index = 7;
    }
}

void send_receive(mpi_com_class * &data, int &id_row, int &id_col, int &id, int &tag_num, MPI_Request* &request, int rows, int cols){
    // Send and receive data between processes
    int cnt = 0;
    int ii = 0;
    int jj = 0;
    int process_address = 0;
    int index = 0;
    int com_i = 0;
    int com_j = 0;
    int com_id = 0;
    int k = 0; 
    int l = 0;
    // Iterate through neighbours of each cell
	for (int i=-1;i<=1;i++){
		for (int j = -1; j <= 1; j++) {
            com_i = (id_row + rows + i)%rows;
		    com_j = (id_col + cols + j)%cols;
			com_id = id_from_index(com_i, com_j, rows, cols); // find id of processor that must communicate with
            if ( com_id < p){ // com_id != id && 
                k = i;
                l = j;
                i_j_to_index(k, l, index); //converts i and j to its corresponding datatype index
                // Set up receive and sen
                MPI_Irecv(MPI_BOTTOM, 1, data->typelist_receive_middle[index], com_id, tag_num, MPI_COMM_WORLD, &request[cnt]);
                cnt++;
                MPI_Isend(MPI_BOTTOM, 1, data->typelist_send_middle[index], com_id, tag_num, MPI_COMM_WORLD, &request[cnt]);
                cnt++;
			}
		}
    }
    MPI_Waitall(cnt, request, MPI_STATUS_IGNORE); // Used for ensuring communication happens correctly
}

void grid_to_file(bool ** &domain, int id, int local_rows, int local_cols){
    // Print all data into a file
    ofstream myfile;
    std::string s0("p_" + std::to_string(id) + "_out.txt" );
    myfile.open(s0, ios::app);
    myfile << "/Iteration/"; // flag used for postprocessing
    for (int i = 0; i<local_rows+2; i++){
        for (int j=0; j<local_cols+2; j++){
            myfile << "," << domain[i][j];
        }
    }
    myfile.close();
}

void info_to_file(int id, int local_rows, int local_cols, int p, int ROWS, int COLS, int rows, int cols, int id_row, int id_col){
    // Print information about the process, separated by &
    ofstream myfile;
    std::string s0("p_" + std::to_string(id) + "_info.txt" );
    myfile.open(s0, ios::trunc);
    myfile << "&" << id ;
    myfile << "&" << local_rows;
    myfile << "&" << local_cols;
    myfile << "&" << p;
    myfile << "&" << ROWS;
    myfile << "&" << COLS;
    myfile << "&" << rows;
    myfile << "&" << cols;
    myfile << "&" << id_row;
    myfile << "&" << id_col;
    myfile.close();
}

void apply_non_periodic_boundary(bool ** &domain, int &id, int &id_row, int &local_rows, int &local_cols, int &id_col, int &rows, int &cols){
    // Apply non-periodic boundaries scenario by making the edges = 0
    if(id_row == 0){
        if((id_row == 0) && (id_col == 0)){
            domain[0][0] = false; 
        }
        if((id_row == 0) && (id_col == cols - 1)){
            domain[0][local_cols + 1] = false; 
        } 
        if(!((!(id_row == 0) && (id_col == 0)) && !(id_row == 0) && (id_col == cols - 1))){
            for(int j = 1; j<local_cols+1; j++){
                domain[0][j] = false;
            }
        }
    }
    // Bottom row
    if(id_row == rows-1){
        if((id_row == rows - 1) && (id_col == 0)){
            domain[local_rows + 1][0] = false; 
        }  
        // Bottom right corner
        if((id_row == rows - 1) && (id_col == cols-1)){
            domain[local_rows + 1][local_cols + 1] = false; 
        } 
        if(!(!((id_row == rows - 1) && (id_col == 0)) && !((id_row == rows - 1) && (id_col == cols-1)))){
            for(int j = 1; j<local_cols+1; j++){
                domain[local_rows+1][j] = false;
            }
        }
    }
    // Left column
    if((id_col+cols)%cols == 0){
        for(int i = 1; i<local_rows+1; i++){
            domain[i][0] = false;
        }     
    }
    // Right column
    if((id_col + 1 + cols)%cols == 0){
        for(int i = 1; i<local_rows+1; i++){
            domain[i][local_cols+1] = false;
        } 
    }
}


int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
    srand((unsigned)time(NULL)+id*p);
    int local_rows, local_cols;
    

    ROWS = 100;
    COLS = ROWS;
    iterations = 100;

    // Assign subdomain to each process
    if(p==1){
        rows = 1;
        cols = 1;
    }
    else{
        find_dimensions(p, rows, cols);
    }

    id_to_index(id, id_row, id_col); // row and col position of index 
    
    // ****************************************************************************************************************************************
    // ****************************************************************************************************************************************

    int cols_reminder = COLS%cols;
    int rows_reminder = ROWS%rows;
    bool uneven_flag = false;
    for(int i = cols-1; i < rows*cols; i = i + cols){ //find local_cols
        if(id == i){ // uneven processor
            local_cols = COLS/cols + cols_reminder; //found our cols of uneven processor
            uneven_flag = true; // prevent from overwriting value -> found uneven process
        }     
        else{   // even processors 
            if(!uneven_flag){
                local_cols = COLS/cols;
            }
        }
    }
    uneven_flag = false;
    for(int i = rows*cols - cols; i <= rows*cols; i++){ //find local_cols
        if(id == i){ // uneven processor
            //cout<< "Found rows uneven processor: " << id << endl;
            local_rows = ROWS/rows + rows_reminder;
            uneven_flag = true; // prevent from overwriting value -> found uneven process
        }     
        else{   // even processors 
            if(!uneven_flag){
                local_rows = ROWS/rows;
            }
        }
    }
    // cout << "Subdomain " << id << " has been split into local_rows: " << local_rows << " local_cols: " << local_cols << endl;
    // cout << "Domain has ROWS: " << ROWS << " COLS: " << COLS << endl;
    // cout << "Domain has been split as rows: " << rows << " cols: " << cols << endl;

// ****************************************************************************************************************************************
// ****************************************************************************************************************************************

    // Create two random subdomains and apply rules of game
    bool **domain = random_subdomain(local_rows, local_cols); //generates padded subdomain
    bool **new_domain = random_subdomain(local_rows, local_cols);

    // // Debugging - Create a couple of empty arrays
    // initialise_array_int(domain, local_rows, local_cols, id);
    // initialise_array_int(domain, local_rows, local_cols, id);
    // if(id == 0){
    // create_glider(domain, local_rows, local_cols, id); // Uncomment to inialise domain with a simple glider 
    // }

// ****************************************************************************************************************************************

    
    // Print processor info into file for post-processing // Uncomment for outputing into a file
    info_to_file(id, local_rows, local_cols, p, ROWS, COLS, rows, cols, id_row, id_col);

    // SET-UP SEND AND RECEIVE:
    mpi_com_class* data = new mpi_com_class(local_rows, local_cols, domain);
    // Initialise data types
    data->mpi_sendtype();
    data->mpi_receivetype();
    MPI_Request* request = new MPI_Request[8*2*p];

    // Send and receive, Apply rules of game, output to file
    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    while(iterations>0){
        start = std::chrono::high_resolution_clock::now();
        if(periodic == false){         // swtich between periodic / non-periodic boundary
            apply_non_periodic_boundary(data->domain, id, id_row, local_rows, local_cols, id_col, rows, cols);
            do_iteration(data->domain, new_domain, local_rows, local_cols);
            apply_non_periodic_boundary(data->domain, id, id_row, local_rows, local_cols, id_col, rows, cols); // and here too
        }
        else{
            do_iteration(data->domain, new_domain, local_rows, local_cols); // apply rules of game
        }
        
        send_receive(data, id_row, id_col, id, tag_num, request, rows, cols); // perform one round of communications between process(es)
        finish = std::chrono::high_resolution_clock::now();
        elapsed += finish - start;

        grid_to_file(data->domain, id, local_rows, local_cols); // Uncomment for outputing into a file
        iterations--;
    }
    data->~mpi_com_class(); // destroy the copes of the array each processor has

    MPI_Barrier(MPI_COMM_WORLD); // wait for all processors to completed the game
    start = std::chrono::high_resolution_clock::now();
    finish = std::chrono::high_resolution_clock::now();
    elapsed += finish - start;
    double elapsed2 = double(elapsed.count());
	double max_time, min_time, diff;

	MPI_Request request2, request3;

	MPI_Iallreduce(&elapsed, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &request2); // Compute the maximum time taken for a process to complete
    MPI_Wait(&request2, MPI_STATUS_IGNORE);
    MPI_Iallreduce(&elapsed, &min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &request3); // Compute the minimum time taken for a process to complet
    MPI_Wait(&request3, MPI_STATUS_IGNORE);

    diff = max_time - min_time; // compute the difference - used for load balancing

    if(id==0){ // hpc std out
        std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
        std::cout << "Domain rows: " << ROWS << " cols: " << COLS << "; Elapsed time: " << max_time << " Diff: " << diff << " s\n";
        std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
    }
    // call all delating routines
    delete [] request;
    for(int j = 0; j<local_rows+2; j++){
            delete[] domain[j];
    }
    delete [] domain;
    for(int j = 0; j<local_rows+2; j++){
        delete[] new_domain[j];
    }
    delete [] new_domain;
    MPI_Finalize();
}