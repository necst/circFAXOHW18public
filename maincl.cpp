
#include <fcntl.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <bitset>
#include <string>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "kseq.h"

#define PORT_BITWIDTH 512
#define VALID_BITS 510
#define DIR_MATRIX_ROW 128
#define VALID_CHARS 126
#define BUFFER_BLOCK 64

const int DIAG_SEGMENT = 42;
const int PASSED_BITS = DIAG_SEGMENT * 3;
const int MAX_PASSED_LENGTH = DIAG_SEGMENT*72;
const int MAX_N_DIAGONALS = MAX_PASSED_LENGTH*73;

/*
Given an event, this function returns the kernel execution time in ms
*/
float getTimeDifference(cl_event event) {
	cl_ulong time_start = 0;
	cl_ulong time_end = 0;
	float total_time = 0.0f;

	clGetEventProfilingInfo(event,
		CL_PROFILING_COMMAND_START,
		sizeof(time_start),
		&time_start,
		NULL);
	clGetEventProfilingInfo(event, 
		CL_PROFILING_COMMAND_END,
		sizeof(time_end),
		&time_end,
		NULL);
	total_time = time_end - time_start;
	return total_time/ 1000000.0; // To convert nanoseconds to milliseconds
}


int load_file_to_memory(const char *filename, char **result) {
	unsigned int size = 0;
	FILE *f = fopen(filename, "rb");
	if (f == NULL) {
		*result = NULL;
		return -1; // -1 means file opening fail
	}
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);
	*result = (char *) malloc(size + 1);
	if (size != fread(*result, sizeof(char), size, f)) {
		free(*result);
		return -2; // -2 means file reading fail
	}
	fclose(f);
	(*result)[size] = 0;
	return size;
}


KSEQ_INIT(int, read)

short char_map(char c) {
	switch (c) {
	case 'A':
		return 0;
	case 'T':
		return 1;
	case 'G':
		return 2;
	case 'C':
		return 3;
	case 'N':
		return 4;
	}
	return -1;
}

char rev_char_map(short i) {
	switch (i) {
	case 0:
		return 'A';
	case 1:
		return 'T';
	case 2:
		return 'G';
	case 3:
		return 'C';
	case 4:
		return 'N';
	}
	return 'E';
}

void fill_sequence_buffer(char* input_seq, uint8_t* buffer, int seq_len,
		int buffer_size) {

	short uint8_offset, uint8_idx, uint8_block, carry, carry_val;
	char cc;
	short cc_as_int;

	carry = 0;
	uint8_offset = 0;
	uint8_idx = BUFFER_BLOCK - 1;
	uint8_block = 0;

	for (int i = 0; i < seq_len; i++) {
		cc = input_seq[i];
		cc_as_int = char_map(cc);

		carry = (uint8_offset < 5) ? 0 : (uint8_offset + 3 - 8);
		carry_val = cc_as_int >> (3 - carry);

		uint8_t val_to_add = cc_as_int << uint8_offset;
		buffer[(BUFFER_BLOCK)*uint8_block + uint8_idx] += val_to_add;

		if (carry) {
			if(uint8_idx - 1 < 0){
				buffer[(BUFFER_BLOCK)*(uint8_block + 1) + PORT_BITWIDTH/8 - 1] = carry_val;
			}else{
				buffer[(BUFFER_BLOCK)*uint8_block + uint8_idx - 1] = carry_val;
			}
		}
		if (carry) {
			uint8_offset = carry;
			uint8_idx--;
		} else {
			uint8_offset += 3;
		}

		if (uint8_offset == 8) {
			uint8_offset = 0;
			uint8_idx--;
		}

		if (uint8_idx == 0 && uint8_offset == 6 ) {
			uint8_offset = 0;
			uint8_idx--;
		}

		if(uint8_idx < 0){
			uint8_idx = BUFFER_BLOCK - 1;
			uint8_block++;
		}
	}

}

int store_seq_as_pairs(FILE* db_fp, FILE* query_fp, cl_context& context, cl_kernel& kernel, cl_command_queue& commands) {

	kseq_t *ref_seq = kseq_init(fileno(db_fp)); //initialize db seq
	kseq_t *query_seq = kseq_init(fileno(query_fp)); //initialize reads seq

	int l, m, err;
	int total_pairs = 0;
	std::vector<kstring_t> reads_strings;
	std::vector<kstring_t> ref_strings;

	while ((l = kseq_read(ref_seq)) >= 0 && (m = kseq_read(query_seq)) >= 0) { //read sequence
		ref_strings.push_back(ref_seq->seq);
		reads_strings.push_back(query_seq->seq);
		total_pairs++;
	}

	std::cout << "Collecting " << total_pairs << " pairs..." << std::endl;

	cl_mem input_query;
	cl_mem input_database;
	cl_mem output_direction_matrixhw;
	cl_mem scalar_parameters;

	for (int i = 0; i < total_pairs; i++) {

		std::cout << "Reading kstring_t for read and reference, iter  = " << i << std::endl;

		int alloc_buff_sizes[3];
		kstring_t read = reads_strings.back();
		reads_strings.pop_back();
		kstring_t ref = ref_strings.back();
		ref_strings.pop_back();

		short n_db_tiles = (ref.l+DIAG_SEGMENT-1)/DIAG_SEGMENT;
		short n_query_tiles = (read.l+DIAG_SEGMENT-1)/DIAG_SEGMENT; 

		short max_horizontal_length = n_db_tiles*DIAG_SEGMENT;
		short max_vertical_length = n_query_tiles*DIAG_SEGMENT;
		int length_tiles_one = max_vertical_length*n_db_tiles;
		int max_current_diags = max_horizontal_length+length_tiles_one;

		int read_buffer_chunks = 1 + (read.l * 3 - 1) / (PORT_BITWIDTH - 2);
		int ref_buffer_chunks = 1 + (ref.l * 3 - 1) / (PORT_BITWIDTH - 2);

		alloc_buff_sizes[0] = (read_buffer_chunks * PORT_BITWIDTH) / 8;
		alloc_buff_sizes[1] = (ref_buffer_chunks * PORT_BITWIDTH) / 8;
		alloc_buff_sizes[2] = MAX_N_DIAGONALS * (PASSED_BITS + 2) / 8;

		uint8_t* q_buffer = (uint8_t*) calloc(1,alloc_buff_sizes[0]);
		uint8_t* db_buff = (uint8_t*) calloc(1,alloc_buff_sizes[1]);
		uint8_t* dir_matrix_buff = (uint8_t*) calloc(1,alloc_buff_sizes[2]);
		short scalar_values[6];

		scalar_values[0]=(read.l+DIAG_SEGMENT-1)/DIAG_SEGMENT;  //num of query tiles
		scalar_values[1]=(ref.l+DIAG_SEGMENT-1)/DIAG_SEGMENT;	//num of db tiles
		scalar_values[2]=(ref.l*3+VALID_BITS-1)/VALID_BITS; //dimension of 512 packs of db data
		scalar_values[3]=(read.l*3+VALID_BITS-1)/VALID_BITS; //dimension of 512 packs of query data
		scalar_values[4]=ref.l;//length of the db
		scalar_values[5]=read.l;//length of the query

		std::cout << "Allocated " << alloc_buff_sizes[0] << " bytes for q_buffer[" << read_buffer_chunks << "]" << std::endl;
		std::cout << "Allocated " << alloc_buff_sizes[1] << " bytes for db_buffer[" << ref_buffer_chunks << "]" << std::endl;
		std::cout << "Allocated " << alloc_buff_sizes[2] << " bytes for dir_matrix_buff[" << max_current_diags << "]" << std::endl;
		

		if (q_buffer == 0 || db_buff == 0 || dir_matrix_buff == 0) {
			std::cout << "Could not allocate input buffers for pair " << i
					<< ", terminating..." << std::endl;
			exit(1);
		}

		std::cout << "Filling buffers with query and ref..." << std::endl;

		fill_sequence_buffer(read.s, q_buffer, read.l,
				(read_buffer_chunks * PORT_BITWIDTH) / 8);
		fill_sequence_buffer(ref.s, db_buff, ref.l,
				(ref_buffer_chunks * PORT_BITWIDTH) / 8);

		//Execute the kernel for this pair
		//

		//create buffers
		input_query = clCreateBuffer(context, CL_MEM_READ_ONLY, alloc_buff_sizes[0],NULL, NULL);
		input_database = clCreateBuffer(context, CL_MEM_READ_ONLY,alloc_buff_sizes[1], NULL, NULL);
		output_direction_matrixhw = clCreateBuffer(context, CL_MEM_READ_WRITE,alloc_buff_sizes[2], NULL, NULL);
		scalar_parameters = clCreateBuffer(context, CL_MEM_READ_WRITE, 6*sizeof(short),NULL, NULL);

		if (!input_query || !input_database || !output_direction_matrixhw || !scalar_parameters) {
			printf("Error: Failed to allocate device memory!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}

		err = clEnqueueWriteBuffer(commands, input_query, CL_TRUE, 0,alloc_buff_sizes[0], q_buffer, 0, NULL, NULL);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to write to source array input_query!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}

		err = clEnqueueWriteBuffer(commands, input_database, CL_TRUE, 0, alloc_buff_sizes[1], db_buff, 0, NULL, NULL);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to write to source array input_database!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}

		err = clEnqueueWriteBuffer(commands, output_direction_matrixhw, CL_TRUE, 0, alloc_buff_sizes[2], dir_matrix_buff, 0, NULL, NULL);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to write to source array output_direction_matrixhw!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		err = clEnqueueWriteBuffer(commands, scalar_parameters, CL_TRUE, 0,6*sizeof(short), scalar_values, 0, NULL, NULL);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to write to source array input_query!\n");
			printf("Test failed\n");
			return EXIT_FAILURE;
		}

		// Set the arguments to our compute kernel
		//
		err = 0;
		printf("set arg 1 \n");
		err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_query);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to set kernel arguments 1! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("set arg 2 \n");
		err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_database);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to set kernel arguments 2! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("set arg 3 \n");
		err |= clSetKernelArg(kernel, 2, sizeof(cl_mem),&output_direction_matrixhw);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to set kernel arguments 3! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("set arg 4 \n");
		err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &scalar_parameters);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to set kernel arguments 4! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}

		// Execute the kernel over the entire range of our 1d input data set
		// using the maximum number of work group items for this device
		//
		cl_event enqueue_kernel;
	#ifdef C_KERNEL
		printf("LAUNCH task \n");
		err = clEnqueueTask(commands, kernel, 0, NULL, &enqueue_kernel);
		printf("task END\n");
	#else
		//remember to define global and local if run with NDRange
		err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, (size_t*) &global,
				(size_t*) &local, 0, NULL, NULL);
	#endif
		if (err) {
			printf("Error: Failed to execute kernel! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("Waiting for enqueue kernel event...\n");
		clWaitForEvents(1, &enqueue_kernel);

		std::cout << "Kernel Time " << getTimeDifference(enqueue_kernel) << std::endl;

		cl_event readDirections, readResults;

		printf("Enqueue read buffer (direction matrix)...\n");
		err = clEnqueueReadBuffer(commands, output_direction_matrixhw, CL_TRUE, 0, alloc_buff_sizes[2], dir_matrix_buff, 0, NULL, &readDirections);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to read array! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}
		printf("Enqueue read buffer values[]...\n");
		err = clEnqueueReadBuffer(commands, scalar_parameters, CL_TRUE, 0, 6 * sizeof(short), scalar_values, 0, NULL, &readResults);
		if (err != CL_SUCCESS) {
			printf("Error: Failed to read array! %d\n", err);
			printf("Test failed\n");
			return EXIT_FAILURE;
		}

		clWaitForEvents(1, &readDirections);
		clWaitForEvents(1, &readResults); 
		printf("Kernel executed iteration = %d\n", i);

		std::cout<<"Max: "<<scalar_values[2] <<" Coordinates: "<< scalar_values[3]<<" "<<scalar_values[4]<<std::endl;
		free(q_buffer);
		free(db_buff);
		free(dir_matrix_buff);
		std::cout<<"Freed local buffers"<<std::endl;
	}

	clReleaseMemObject(input_database);
	clReleaseMemObject(input_query);
	clReleaseMemObject(output_direction_matrixhw);
	clReleaseMemObject(scalar_parameters);


	return total_pairs;
}

void write_random_fasta_sequences(char* file_name, int n_seq, int max_seq_len) {

	FILE* fp = fopen(file_name, "w");

	

	if (fp == nullptr) {
		printf("failed to open the file to write\n");
		perror("fopen");
		exit(1);
	}
	long seq_serial_id = 0;
	for (int i = 0; i < n_seq; i++) {
		std::string sequence_str;
		std::string seq_name = std::string(">") + std::to_string(seq_serial_id)
				+ std::string("\n");
		int seq_len = max_seq_len;

		for (int j = 0; j < seq_len; j++) {
			char cc = rev_char_map(std::rand() % 5);
			sequence_str += cc;
		}
		sequence_str += "\n\n";

		fputs(seq_name.c_str(), fp);
		fputs(sequence_str.c_str(), fp);
		seq_serial_id++;
	}
	fclose(fp);
}

int main(int argc, char* argv[]) {

	printf("starting HOST code \n");
	fflush(stdout);
	std::srand(std::time(nullptr));

	write_random_fasta_sequences(argv[2], 1, 3024);
	write_random_fasta_sequences(argv[3], 1, 3024);

	if (argc < 4) {
		std::cout << "No input arguments" << std::endl;
	}

	FILE* db_fp = fopen(argv[2], "r");

	if (db_fp == 0) {
		printf("failed to open the reference FASTA file\n");
		perror("fopen");
		exit(1);
	}

	FILE* query_fp = fopen(argv[3], "r");

	if (query_fp == 0) {
		printf("failed to open the reads FASTA file\n");
		perror("fopen");
		exit(1);
	}

	std::string mode;

	if (argc < 5) {
		mode = "p"; //by default input references and reads are aligned by pairs, in a single-end fashion
	} else {
		mode = "-p";
	}

	
	mode = "p";

	//Create the necessary infrastructure to run the kernel
	//

	int err;                            // error code returned from api calls


	cl_platform_id platform_id;         // platform id
	cl_device_id device_id;             // compute device id
	cl_context context;                 // compute context
	cl_command_queue commands;          // compute command queue
	cl_program program;                 // compute program
	cl_kernel kernel;                   // compute kernel

	char cl_platform_vendor[1001];
	char cl_platform_name[1001];

	// Connect to first platform
	//
	printf("GET platform \n");
	err = clGetPlatformIDs(1, &platform_id, NULL);
	if (err != CL_SUCCESS) {
		printf("Error: Failed to find an OpenCL platform!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	printf("GET platform vendor \n");
	err = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, 1000,
			(void *) cl_platform_vendor, NULL);
	if (err != CL_SUCCESS) {
		printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	printf("CL_PLATFORM_VENDOR %s\n", cl_platform_vendor);
	printf("GET platform name \n");
	err = clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, 1000,
			(void *) cl_platform_name, NULL);
	if (err != CL_SUCCESS) {
		printf("Error: clGetPlatformInfo(CL_PLATFORM_NAME) failed!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	printf("CL_PLATFORM_NAME %s\n", cl_platform_name);

	// Connect to a compute device
	//
	int fpga = 1;

	printf("get device FPGA is %d  \n", fpga);
	err = clGetDeviceIDs(platform_id,
			fpga ? CL_DEVICE_TYPE_ACCELERATOR : CL_DEVICE_TYPE_CPU, 1,
			&device_id, NULL);
	if (err != CL_SUCCESS) {
		printf("Error: Failed to create a device group!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	// Create a compute context
	//
	printf("create context \n");
	context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	if (!context) {
		printf("Error: Failed to create a compute context!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	// Create a command commands
	//
	printf("create queue \n");
	commands = clCreateCommandQueue(context, device_id,
	CL_QUEUE_PROFILING_ENABLE, &err);
	if (!commands) {
		printf("Error: Failed to create a command commands!\n");
		printf("Error: code %i\n", err);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	int status;

	// Create Program Objects
	//

	// Load binary from disk
	unsigned char *kernelbinary;
	char *xclbin = argv[1];
	printf("loading %s\n", xclbin);
	int n_i = load_file_to_memory(xclbin, (char **) &kernelbinary);
	if (n_i < 0) {
		printf("failed to load kernel from xclbin: %s\n", xclbin);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}
	size_t n = n_i;
	// Create the compute program from offline
	printf("create program with binary \n");
	program = clCreateProgramWithBinary(context, 1, &device_id, &n,(const unsigned char **) &kernelbinary, &status, &err);
	if ((!program) || (err != CL_SUCCESS)) {
		printf("Error: Failed to create compute program from binary %d!\n",
				err);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	// Build the program executable
	//
	printf("build program \n");
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];
		printf("Error: Failed to build program executable!\n");
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
				sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	// Create the compute kernel in the program we wish to run
	//
	printf("create kernel \n");
	kernel = clCreateKernel(program, "kernel", &err);
	if (!kernel || err != CL_SUCCESS) {
		printf("Error: Failed to create compute kernel!\n");
		printf("Test failed\n");
		return EXIT_FAILURE;
	}

	//Execute the kernel on all input queries in the selected mode
	//

	if (mode == "p") {
		store_seq_as_pairs(db_fp, query_fp, context, kernel, commands);
	} else if (mode == "all") {
		std::cout << "All mode not supported yet" << std::endl;
	} else {
		std::cout << "Error: no valid mode selected, " << mode << std::endl;
	}

	fclose(query_fp);
	fclose(db_fp);
	
	clReleaseProgram(program);
	clReleaseKernel(kernel);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);
	
	return EXIT_SUCCESS;
}

