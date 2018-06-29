#include <ap_utils.h>
#include <iostream>
#include <bitset>

#define AP_INT_MAX_W 9000
#include <ap_int.h>

//scores attributed for the alignment
const short STARTGAP = 3;
const short EXTENDGAP = 1;
const short MATCH = 2;
const short MISMATCH = -2;
//a bad mismatch happens when two bases are opposites
//cases of bad mismatches are: A-T and C-G
const short BAD_MISMATCH = -3;
//Bases can be A C G T and N
//when they are N this means that they have been considered ambiguous
//when the sequences was extracted, therefore in order to indicate this
//ambiguity if there is a mismatch with an N it doesn't count as bad as
//a regular mismatch (this behavior is the same adopted by software solutions)
const short N_MISMATCH = -1;
//number of possible bases
const short N_BASES = 5;

//lengths for the anti-diagonal we are calculating
//query and db are referred to the sequence we want to align
//and the reference genome respectively
//
//number of database elements considered in one cycle
const short DB_LEN = 42;
//length of the anti-diagonal segment considered each cycle
const short DIAG_SEGMENT = 42;
//number of query elements considered in one cycle
const short QUERY_LEN = 42;

//port bit-width of the FPGA, this determines how much data
//we can send/receive per time
const short PORT_BITWIDTH = 512;
//closest multiple of 3 (number of bits used do codify the bases) to PORT_BITWIDTH
const short VALID_BITS = PORT_BITWIDTH - 2;
//Passed bits indicates how many bits are taken for an anti-diagonal segment
const short PASSED_BITS = DIAG_SEGMENT*3;
//max_passed_length indicates the max number elements for query and database we can process
//in future releases we plan to extend this number in order to support very long reads
const short MAX_PASSED_LENGTH = DIAG_SEGMENT*72;
//max_n_diagonals indicates the maximum number of anti-diagonal segments we have to calculate
const int MAX_N_DIAGONALS = MAX_PASSED_LENGTH*73;
//number of elements in the circular buffer
const short DIM_BUFFER = 2;
//maximum number of elements that we have in BRAM at the same time
const short DIM_BUFFER_LOCAL = 250;
//this indicates how much blocks of VALID_BITS we need to store query and
//database on BRAM
const short LENGTH_BUFFERS_DB_QUERY = (MAX_PASSED_LENGTH*3+VALID_BITS-1)/VALID_BITS;
//factor for partitioning the dependencies arrays in BRAM
const short BLOCK_FACTOR_SCORE_VALUES = 2; 
//W is the weights matrix, this permits an easier access to the scores when the scoring kernel
//executes because all the values are pre-calculated
const short W[N_BASES][N_BASES] = { { MATCH, BAD_MISMATCH, MISMATCH, MISMATCH,
		N_MISMATCH }, { BAD_MISMATCH, MATCH, MISMATCH, MISMATCH, N_MISMATCH }, {
		MISMATCH, MISMATCH, MATCH, BAD_MISMATCH, N_MISMATCH }, { MISMATCH,
		MISMATCH, BAD_MISMATCH, MATCH, N_MISMATCH }, { N_MISMATCH, N_MISMATCH,
		N_MISMATCH, N_MISMATCH, MATCH } };


//extern "C"{
	/*
	*	@brief: scoringKernel is the part of our kernel that computes each anti-diagonal segment
	*			it is continuously called by the main loop in order to execute the entire direction matrix
	*
	*	@param: query that references the part of query we are considering for this part of the alignment
	*	@param: db that references the part of database we are considering for this part of the alignment
	*	@param: H_curr which indicates the current dependencies and is used to consider north or west values, it is considered if a gap is started
	*	@param: H_prev which indicates the previous dependencies and is used to consider the north-west values
	*	@param: variable F_val and E_val represent dependencies on north and west respectively and they are considered a gap is extended
	*	@param: maxBuffer contains the values of the alignment score, when we end the computation we consider only the maximum value for the alignment
	*	@param: coord_max_buffer contains the values that indicates where in the segment of diagonal calculated the maximum value is
	*	@param: diag_buffer contains the information regarding which segment of anti-diagonal the maximum is stored
	*	@param: index_diag contains the information regarding which segment of anti-diagonal will be calculated
	*	@param: position_antidiagonal contains information regarding where the segment calculated is
	*	@param: pos_max_diag_buffer contains information regarding where the segment of the maximum is
	*	@param: rowToReach and rowToStart indicate the limits of the anti-diagonal segment to calculate
	*
	*/
void scoringKernel(ap_uint<PASSED_BITS> *query, ap_uint<PASSED_BITS> *db,
		ap_uint<PASSED_BITS> *diagonalDirectionsVector,
		short H_curr[DIAG_SEGMENT + 1], short H_prev[DIAG_SEGMENT + 1],
		short F_val[DIAG_SEGMENT + 1], short E_val[DIAG_SEGMENT + 1],
		short maxBuffer[DIAG_SEGMENT],
		short coord_max_buffer[DIAG_SEGMENT], short diag_buffer[DIAG_SEGMENT], int index_diag,
		short position_antidiagonal,
		short pos_max_diag_buffer[DIAG_SEGMENT],
		short rowToReach, short rowToStart) {
#pragma HLS INLINE

	short maxScore = 0;

	ap_uint<3> dir_buffer[DIAG_SEGMENT];

	ap_uint<PASSED_BITS> diagonal_local;

// LOOP that executes the single anti-diagonal given the inputs provided above
	diag_loop: for (int i = DIAG_SEGMENT; i > 0; i--) {
#pragma HLS UNROLL

		short E_local;
		short F_local;
		short H_curr_local;
		short H_prev_local;
		//we start with the rightmost cell of the anti-diagonal
		H_curr_local = H_curr[i];
		E_local = E_val[i] - EXTENDGAP;

		ap_uint<3> directions = 0;
		bool in_range = (i <= rowToReach && i > rowToStart);
		bool bool_E, bool_F;
		//mask if the value is in range we kept the value, otherwise we set it to zero
		short bit_mask = in_range ? -1 : 0;
		//take the north west value and calculate the new one
		short H_nw = H_prev[i - 1]
				+ W[db->range((i - 1) * 3 + 2, (i - 1) * 3)][query->range(
						(i - 1) * 3 + 2, (i - 1) * 3)];
		//calculate the new north and west values
		short t2 = H_curr_local - STARTGAP;
		short t1_2 = F_val[i - 1] - EXTENDGAP;
		short t2_2 = H_curr[i - 1] - STARTGAP;
		short prevMax = maxBuffer[i - 1];
		short newMaxVal;
		short coord = coord_max_buffer[i - 1];
		int id_diag = diag_buffer[i - 1];
		short pos_id = pos_max_diag_buffer[i - 1];
		bool_E = E_local > t2;
		bool_F = t1_2 > t2_2;

		//update the value of E for this column
		E_local = bool_E ? E_local : t2;
		E_val[i] = E_local & bit_mask;

		//update the value of F for this column
		F_local = bool_F ? t1_2 : t2_2;
		F_val[i] = F_local & bit_mask;
		bool bool_HF = H_nw > F_local;
		short HF = bool_HF ? H_nw : F_local;

		bool bool_HF0 = HF > 0;
		short HF0 = bool_HF0 ? HF : 0;

		bool bool_max = HF0 > E_local;
		short max = bool_max ? HF0 : E_local;

		bool newMax = (max > prevMax) && in_range;

		//update the value of H_prev for next cycle
		H_prev[i] = H_curr_local & bit_mask;
		//set the bits for the direction taken
		bool b1 = max == E_local;
		bool b2 = max == F_local;
		bool b3 = max == H_nw;

		directions.range(0, 0) = b1;
		directions.range(1, 1) = b2;
		directions.range(2, 2) = b3;

		newMaxVal = newMax ? max : prevMax;
		coord_max_buffer[i - 1] = newMax ? (i - 1) : coord;
		diag_buffer[i - 1] = newMax ? (index_diag) : id_diag;
		pos_max_diag_buffer[i - 1] = newMax ? position_antidiagonal : pos_id;
		dir_buffer[i - 1] = directions;
		//the value in H is always the max between those calculated
		H_curr[i] = max;
		maxBuffer[i - 1] = newMaxVal;
	}

	save_directions_for: for (int i = 0; i < DIAG_SEGMENT; i++) {
#pragma HLS UNROLL
		diagonal_local.range(3 * i + 2, 3 * i) = dir_buffer[i];
	}
	*diagonalDirectionsVector = diagonal_local;

}
	/*
	*	@brief: write_to_cb is the function that describes the behavior of how to write into the circular buffer
	*			it is called by the setNextValue function
	*
	*
	*	@param: index_last_written_cb indicates the last written cells of the circular buffer
	*	@param: elements in the circular buffer indicated how many elements are in the buffer
	*	@param: q_block_cb stores in the circular buffer which block 510 bits of query we have to consider of the segments that will be calculated
	*	@param:	endQueryIndex_cb stores in the circular buffer the end bit of query to consider of the segments that will be calculated
	*	@param: db_block_cb stores in the circular buffer which block 510 bits of database we have to consider of the segments that will be calculated
	*	@param:	endDbIndex_cb stores in the circular buffer the end bit of database to consider of the segments that will be calculated
	*	@param: tmp_q_block indicates which block 510 bits of query we have considered at the last execution
	*	@param: tmp_q_end indicates the end bit of query to considered at the last execution
	*	@param:	temp_cycle indicates how many cycles for the tile have been executed
	*	@param: tile_type indicates the tile_type just executed (1 or 0)
	*	@param:	split_down indicates if we have to split (create a tile 1 after a tile 0)
	*	@param: position_antidiagonal_cb stores in the circular buffer the positions of the anti-diagonal segments that will be calculated
	*	@param:	pos_to_assign indicates the position of the last calculated segment
	*	@param: tile_type_circular_buffer stores in the circular buffer the tile types of the segments that will be calculated
	*	@param:	iteration_cycle_cb stores in the circular buffer the number of cycles that a tile has been executed which are related to the segments that will be calculated
	*	@param: temp_H_curr, temp_H_prev, temp_E_val, temp_F_val contain the dependencies of the last calculated segment
	*	@param:	temp_rowToStart and temp_rowToReach indicate the limits of the anti-diagonal segment to calculate
	*	@param: len_db and len_query are the complete lengths of the strings to align
	*	@param: circular_buffer_H_curr,circular_buffer_H_prev,circular_buffer_F_val,circular_buffer_E_val store in the circular buffer the dependencies for the tiles that we'll compute
	*	@param: dep_H_curr_tmp, dep_H_prev_tmp, dep_E_val_tmp store the dependency coming from another tile, in this way we can pass dependencies between segments of the same anti-diagonal
	*
	*/

void write_to_cb(bool* index_last_written_cb, short* elements_in_cb,
		short rowToStart_circular_buffer[DIM_BUFFER],
		short rowToReach_circular_buffer[DIM_BUFFER],
		short q_block_cb[DIM_BUFFER],
		short endQueryIndex_cb[DIM_BUFFER], short db_block_cb[DIM_BUFFER],
		short endDbIndex_cb[DIM_BUFFER],
		short* tmp_q_block, short* tmp_q_end, short tmp_db_block,
		short tmp_db_end, short temp_cycle,
		bool tile_type,
		bool split_down,
		short position_antidiagonal_cb[DIM_BUFFER],
		short pos_to_assign,
		bool tile_type_circular_buffer[DIM_BUFFER],
		short iteration_cycle_cb[DIM_BUFFER],
		short temp_H_curr[DIAG_SEGMENT + 1],
		short temp_H_prev[DIAG_SEGMENT + 1], short temp_E_val[DIAG_SEGMENT + 1],
		short temp_F_val[DIAG_SEGMENT + 1],
		short temp_rowToStart,
		short temp_rowToReach,
		short len_db,
		short len_query,
		short circular_buffer_H_curr[DIM_BUFFER][DIAG_SEGMENT + 1],
		short circular_buffer_H_prev[DIM_BUFFER][DIAG_SEGMENT + 1],
		short circular_buffer_E_val[DIM_BUFFER][DIAG_SEGMENT + 1],
		short circular_buffer_F_val[DIM_BUFFER][DIAG_SEGMENT + 1],
		short dep_H_curr_tmp, short dep_H_prev_tmp, short dep_E_val_tmp) {
#pragma HLS INLINE
	
	bool access_offset = 0;
	
	*elements_in_cb += 1;
	//update CB indexes
	*index_last_written_cb = !(*index_last_written_cb);

	//update the limits of the anti-diagonal segment to compute
	//based on the tile type and on the number of cycles we have executed a tile
	if(temp_cycle<=DIAG_SEGMENT&&temp_cycle<len_query){
		temp_rowToReach++;
	}
	
	if(temp_cycle>=len_query&&tile_type) temp_rowToReach--;
	if(!tile_type&&temp_cycle>=len_db) temp_rowToStart++;

	//write CB values for the new tile
	rowToStart_circular_buffer[*index_last_written_cb] = temp_rowToStart;
	rowToReach_circular_buffer[*index_last_written_cb] = temp_rowToReach;
	tile_type_circular_buffer[*index_last_written_cb] = tile_type;

	//select the correct update of the piece of query to consider
	if (tile_type) {
		
		if(*tmp_q_end + 3 == PORT_BITWIDTH){
			*tmp_q_block += 1;
			*tmp_q_end  = 2;
		}else{
			*tmp_q_end += 3;
		}
		access_offset = 1;
	}

	q_block_cb[*index_last_written_cb] = *tmp_q_block;
	endQueryIndex_cb[*index_last_written_cb] = *tmp_q_end;
	db_block_cb[*index_last_written_cb] = tmp_db_block;
	endDbIndex_cb[*index_last_written_cb] = tmp_db_end;
	position_antidiagonal_cb[*index_last_written_cb]=pos_to_assign;
	for (int i = 0; i < DIAG_SEGMENT; i++) {
#pragma HLS UNROLL
		circular_buffer_H_curr[*index_last_written_cb][i] = temp_H_curr[i
				+ access_offset];
		circular_buffer_H_prev[*index_last_written_cb][i] = temp_H_prev[i
				+ access_offset];
		circular_buffer_E_val[*index_last_written_cb][i] = temp_E_val[i
				+ access_offset];
		circular_buffer_F_val[*index_last_written_cb][i] = temp_F_val[i
				+ access_offset];
	}


	circular_buffer_H_curr[*index_last_written_cb][DIAG_SEGMENT] = dep_H_curr_tmp;
	circular_buffer_H_prev[*index_last_written_cb][DIAG_SEGMENT] = dep_H_prev_tmp;
	circular_buffer_E_val[*index_last_written_cb][DIAG_SEGMENT] = dep_E_val_tmp;
}

void write_to_BRAM(short* index_last_written_mem,
		short rowToStart_local_buffer[DIM_BUFFER_LOCAL],
		short rowToReach_local_buffer[DIM_BUFFER_LOCAL],
		short q_block_lb[DIM_BUFFER_LOCAL],
		short endQueryIndex_lb[DIM_BUFFER_LOCAL],
		short db_block_lb[DIM_BUFFER_LOCAL],
		short endDbIndex_lb[DIM_BUFFER_LOCAL],short* tmp_q_block,
		short *tmp_q_end, short tmp_db_block, short tmp_db_end, short temp_cycle,
		bool tile_type,
		bool split_down,
		short position_antidiagonal_lb[DIM_BUFFER_LOCAL],
		short pos_to_assign,
		bool tile_type_local_buffer[DIM_BUFFER_LOCAL],
		short iteration_cycle_lb[DIM_BUFFER_LOCAL],
		short temp_H_curr[DIAG_SEGMENT + 1],
		short temp_H_prev[DIAG_SEGMENT + 1], short temp_E_val[DIAG_SEGMENT + 1],
		short temp_F_val[DIAG_SEGMENT + 1],
		short temp_rowToStart,
		short temp_rowToReach,
		short len_db,
		short len_query,
		short local_buffer_H_curr[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],
		short local_buffer_H_prev[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],
		short local_buffer_E_val[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],
		short local_buffer_F_val[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],
		short dep_H_curr_tmp, short dep_H_prev_tmp, short dep_E_val_tmp) {
#pragma HLS INLINE

	bool access_offset = 0;

	//write BRAM values for the new tile

	if(temp_cycle>=len_query&&tile_type) temp_rowToReach--;
	if(!tile_type&&temp_cycle>=len_db) temp_rowToStart++;
	rowToStart_local_buffer[*index_last_written_mem] = temp_rowToStart;
	rowToReach_local_buffer[*index_last_written_mem] = temp_rowToReach;
	tile_type_local_buffer[*index_last_written_mem] = tile_type;

	if (tile_type) {

		if(*tmp_q_end + 3 == PORT_BITWIDTH){
			*tmp_q_block += 1;
			*tmp_q_end  = 2;
		}else{
			*tmp_q_end += 3;
		}
		access_offset = 1;
	}

	q_block_lb[*index_last_written_mem] = *tmp_q_block;
	endQueryIndex_lb[*index_last_written_mem] = *tmp_q_end;
	db_block_lb[*index_last_written_mem] = tmp_db_block;
	endDbIndex_lb[*index_last_written_mem] = tmp_db_end;
	position_antidiagonal_lb[*index_last_written_mem] = pos_to_assign;
	iteration_cycle_lb[*index_last_written_mem] = temp_cycle;

	for (int i = 0; i < DIAG_SEGMENT; i++) {
#pragma HLS UNROLL
		local_buffer_H_curr[*index_last_written_mem][i] = temp_H_curr[i
				+ access_offset];
		local_buffer_H_prev[*index_last_written_mem][i] = temp_H_prev[i
				+ access_offset];
		local_buffer_E_val[*index_last_written_mem][i] = temp_E_val[i
				+ access_offset];
		local_buffer_F_val[*index_last_written_mem][i] = temp_F_val[i
				+ access_offset];
	}

	local_buffer_H_curr[*index_last_written_mem][DIAG_SEGMENT] = dep_H_curr_tmp;
	local_buffer_H_prev[*index_last_written_mem][DIAG_SEGMENT] = dep_H_prev_tmp;
	local_buffer_E_val[*index_last_written_mem][DIAG_SEGMENT] = dep_E_val_tmp;


}
	/*
	*	@brief: setNextValues is the part of the kernel which focuses on moving the data in the buffers and
	*			serve data to the scoring kernel in order, it is called after the scoringKernel function and receives all the
	*			arrays that are handled in the functions above. Depending on the situation setNextValues calls write_to_cb or write_to_BRAM
	*			and also makes other operations to move data in the correct place. Data is always put in the circular buffer before execution,
	*			therefore if some data is stored in BRAM and needs to be executed it is moved into the circular buffer as soon as possible, then
	*			the function proceeds to save the last calculated data into the most appropriate buffer.
	*
	*
	*	@param: totQuery contains the whole query
	*	@param: totDb contains the whole database
	*	@param: q_block_BRAM contains the number of query block read from the BRAM
	*	@param: endQuery_index_BRAM contains last bit of query read from the BRAM
	*	@param: db_block_BRAM the number of database block read from the BRAM
	*	@param: endDb_index_BRAM contains last bit of database read from the BRAM
	*	@param: rowToStart_BRAM and rowToReach_BRAM contain the limits of the diagonal computation read from BRAM
	*	@param: query_to_pass_BRAM contains the information on the query to align read from BRAM
	*	@param: db_to_pass_BRAM contains the information on the database to align read from BRAM
	*	@param: tile_type_BRAM contains the information on the tile type read from BRAM
	*	@param: F_val_BRAM, E_val_BRAM, H_curr_BRAM, H_prev_BRAM contain all the dependencies read from BRAM
	*	@param: H_curr_dep_BRAM, H_prev_dep_BRAM, E_val_dep_BRAM contain the dependency for the next segment of antidiagonal read from BRAM
	*	@param: iteration_cycle_BRAM contains the number of cycles the tile has been executed read from BRAM
	*	@param: position_antidiagonal_BRAM contains the information regarding the position of the anti-diagonal segment read from BRAM
	*
	*
	*/
void setNextValues(ap_uint<VALID_BITS> totQuery[MAX_PASSED_LENGTH/VALID_BITS], ap_uint<VALID_BITS> totDb[MAX_PASSED_LENGTH/VALID_BITS],
		short q_block_BRAM, short endQuery_index_BRAM,
		short db_block_BRAM, short endDb_index_BRAM,
		short rowToStart_BRAM,
		short rowToReach_BRAM,
		ap_uint<PASSED_BITS+2> query_to_pass_BRAM,
		ap_uint<PASSED_BITS+2> db_to_pass_BRAM,
		bool tile_type_BRAM,
		short F_val_BRAM[DIAG_SEGMENT + 1],
		short E_val_BRAM[DIAG_SEGMENT + 1],
		short H_curr_BRAM[DIAG_SEGMENT + 1],
		short H_prev_BRAM[DIAG_SEGMENT + 1],
		short H_curr_dep_BRAM,
		short H_prev_dep_BRAM,
		short E_val_dep_BRAM,
		short iteration_cycle_BRAM,
		short position_antidiagonal_BRAM,

		//these are the same values used by the write_to_BRAM and write_to_cb functions
		short q_block_cb[DIM_BUFFER],
		short endQueryIndex_cb[DIM_BUFFER],	
		short db_block_cb[DIM_BUFFER],	
		short endDbIndex_cb[DIM_BUFFER],	
		short q_block_lb[DIM_BUFFER_LOCAL],  
		short endQueryIndex_lb[DIM_BUFFER_LOCAL],		
		short db_block_lb[DIM_BUFFER_LOCAL],		
		short endDbIndex_lb[DIM_BUFFER_LOCAL],	

		bool tile_type_local_buffer[DIM_BUFFER_LOCAL], 
		bool tile_type_circular_buffer[DIM_BUFFER], 
		short position_antidiagonal_lb[DIM_BUFFER_LOCAL],
		short position_antidiagonal_cb[DIM_BUFFER],
		bool *index_last_written_cb,
		//index_to_write indicates the index to write in BRAM
		short *index_to_write,

		short *last_read_from_memory,
		bool to_read_from_cb,

		ap_uint<PASSED_BITS> circular_buffer_query[DIM_BUFFER],				   
		ap_uint<PASSED_BITS> circular_buffer_db[DIM_BUFFER],				   
		short circular_buffer_H_curr[DIM_BUFFER][DIAG_SEGMENT + 1],			
		short circular_buffer_H_prev[DIM_BUFFER][DIAG_SEGMENT + 1],			
		short circular_buffer_F_val[DIM_BUFFER][DIAG_SEGMENT + 1],			
		short circular_buffer_E_val[DIM_BUFFER][DIAG_SEGMENT + 1],			
		short rowToReach_circular_buffer[DIM_BUFFER],				   
		short rowToStart_circular_buffer[DIM_BUFFER],				   

		ap_uint<PASSED_BITS> local_buffer_query[DIM_BUFFER_LOCAL],				   
		ap_uint<PASSED_BITS> local_buffer_db[DIM_BUFFER_LOCAL],				   
		short local_buffer_H_curr[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],		
		short local_buffer_H_prev[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],		
		short local_buffer_F_val[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],		
		short local_buffer_E_val[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1],		
		short rowToReach_local_buffer[DIM_BUFFER_LOCAL],				   
		short rowToStart_local_buffer[DIM_BUFFER_LOCAL],				   

		short H_curr_dep_cb[DIM_BUFFER], short H_prev_dep_cb[DIM_BUFFER], short E_val_dep_cb[DIM_BUFFER],
		short H_curr_dep_lb[DIM_BUFFER_LOCAL], short H_prev_dep_lb[DIM_BUFFER_LOCAL], short E_val_dep_lb[DIM_BUFFER_LOCAL],

		short *elements_in_cb,
		short *elements_in_BRAM,
		//general pos indicated the position we are in, this value increases each time a new tile start the computation
		short *general_pos,
		short iteration_cycle_cb[DIM_BUFFER],
		short iteration_cycle_lb[DIM_BUFFER_LOCAL],
		//n_query_tiles and n_db_tiles are the quantity number of tiles (len_db and len_query divided by DIAG_SEGMENT) that are considered
		short n_query_tiles, short n_db_tiles, short len_db, short len_query) {
#pragma HLS INLINE

	//decrease number of elements in cb
	*elements_in_cb-=1;
	short temp_it_cycle = iteration_cycle_cb[to_read_from_cb];

	//get old db and query slices
	ap_uint<PASSED_BITS> temp_query = circular_buffer_query[to_read_from_cb];
	ap_uint<PASSED_BITS> temp_db = circular_buffer_db[to_read_from_cb];

	//get range values for PEs
	short temp_rowToStart = rowToStart_circular_buffer[to_read_from_cb];
	short temp_rowToReach = rowToReach_circular_buffer[to_read_from_cb];

	//get old indexes for absolute position of query and db stretches from the circular buffer
	short tmp_db_block_d_BRAM = db_block_cb[to_read_from_cb];
	short tmp_q_block_d_BRAM = q_block_cb[to_read_from_cb];
	short tmp_db_end_d_BRAM = endDbIndex_cb[to_read_from_cb];
	short tmp_q_end_d_BRAM = endQueryIndex_cb[to_read_from_cb];
	short tmp_db_block_r = db_block_cb[to_read_from_cb];
	short tmp_q_block_r = q_block_cb[to_read_from_cb];
	short tmp_db_end_r = endDbIndex_cb[to_read_from_cb];
	short tmp_q_end_r = endQueryIndex_cb[to_read_from_cb];
	short tmp_db_block_d_CB = db_block_cb[to_read_from_cb];
	short tmp_q_block_d_CB = q_block_cb[to_read_from_cb];
	short tmp_db_end_d_CB = endDbIndex_cb[to_read_from_cb];
	short tmp_q_end_d_CB = endQueryIndex_cb[to_read_from_cb];
	short tmp_max_antg_pos = position_antidiagonal_cb[to_read_from_cb];

	//dependencies from previous diagonal chunk
	short dep_H_curr_tmp, dep_H_prev_tmp, dep_E_val_tmp;

	dep_H_curr_tmp = H_curr_dep_cb[to_read_from_cb];
	dep_H_prev_tmp = H_prev_dep_cb[to_read_from_cb];
	dep_E_val_tmp = E_val_dep_cb[to_read_from_cb];

	//get the type of tile
	bool temp_tile_type = tile_type_circular_buffer[to_read_from_cb];

	//copy old diagonals values
	short temp_H_prev[DIAG_SEGMENT + 1], temp_H_curr[DIAG_SEGMENT + 1],
			temp_F_val[DIAG_SEGMENT + 1], temp_E_val[DIAG_SEGMENT + 1];
#pragma HLS ARRAY_PARTITION variable=temp_H_prev complete dim=1
#pragma HLS ARRAY_PARTITION variable=temp_H_curr complete dim=1
#pragma HLS ARRAY_PARTITION variable=temp_F_val complete dim=1
#pragma HLS ARRAY_PARTITION variable=temp_E_val complete dim=1

	for (int i = 0; i < DIAG_SEGMENT + 1; i++) {
#pragma HLS UNROLL
		temp_H_curr[i] = circular_buffer_H_curr[to_read_from_cb][i];
		temp_H_prev[i] = circular_buffer_H_prev[to_read_from_cb][i];
		temp_E_val[i] = circular_buffer_E_val[to_read_from_cb][i];
		temp_F_val[i] = circular_buffer_F_val[to_read_from_cb][i];
	}

	//check if there is any value in BRAM, if that's the case put it on the circular buffer
	if (*elements_in_BRAM > 0) {
		
		*index_last_written_cb = !(*index_last_written_cb);

		//update all circular buffer values at the correct index with the values read from BRAM
		rowToStart_circular_buffer[*index_last_written_cb] = rowToStart_BRAM;
		rowToReach_circular_buffer[*index_last_written_cb] = rowToReach_BRAM;
		position_antidiagonal_cb[*index_last_written_cb] = position_antidiagonal_BRAM;
		q_block_cb[*index_last_written_cb] = q_block_BRAM;
		db_block_cb[*index_last_written_cb] = db_block_BRAM;
		endQueryIndex_cb[*index_last_written_cb] = endQuery_index_BRAM;
		endDbIndex_cb[*index_last_written_cb] = endDb_index_BRAM;

		//set the correct query and db on the circular buffer
		circular_buffer_query[*index_last_written_cb] = query_to_pass_BRAM;
		circular_buffer_db[*index_last_written_cb] = db_to_pass_BRAM;
		//set tile type on the circular buffer
		tile_type_circular_buffer[*index_last_written_cb] = tile_type_BRAM;
		//copy BRAM dependencies on the circular buffer
		for (int i = 0; i < DIAG_SEGMENT + 1; i++) {
			circular_buffer_H_curr[*index_last_written_cb][i] = H_curr_BRAM[i];
			circular_buffer_H_prev[*index_last_written_cb][i] = H_prev_BRAM[i];
			circular_buffer_E_val[*index_last_written_cb][i] = E_val_BRAM[i];
			circular_buffer_F_val[*index_last_written_cb][i] = F_val_BRAM[i];
		}
		//set dependencies for the next segment to compute
		H_curr_dep_cb[*index_last_written_cb] = H_curr_dep_BRAM;
		H_prev_dep_cb[*index_last_written_cb] = H_prev_dep_BRAM;
		E_val_dep_cb[*index_last_written_cb] = E_val_dep_BRAM;
		//copy iteration cycle
		iteration_cycle_cb[*index_last_written_cb] = iteration_cycle_BRAM;

		//update status counter for elements in the circular buffer and elements read from memory
		*elements_in_cb +=1;
		*last_read_from_memory += 1;
		if(*last_read_from_memory == DIM_BUFFER_LOCAL) *last_read_from_memory=0;
		*elements_in_BRAM-=1;
	}

	/*  Tiles are organized as follows:
	*
	*	|-------|-------|-------|
	*	|	0	|	0	|	0	|
	*	|-------|-------|-------|
	*	|	1	|	1	|	1	|
	*	|-------|-------|-------|
	*	|	1	|	1	|	1	|
	*	|-------|-------|-------|
	*
	*
	*/
	//pos_to_assign is the same as the last calculated
	short pos_to_assign=tmp_max_antg_pos;
	// if the tile of type zero reaches a multiple of DIAG_SEGMENT cycles it starts a tile of type one that goes down
	// the tile of type zero then proceeds its execution to the right
	bool split = (((temp_it_cycle + 1) % (DIAG_SEGMENT)) == 0);
	bool split_down = ((split && (temp_tile_type==0)));
	short temp_down_cycle = split_down ? DIAG_SEGMENT : (temp_it_cycle+1);
	short temp_right_cycle = temp_it_cycle+1;
	if(*elements_in_cb < DIM_BUFFER){
		if(split_down||temp_tile_type){
			iteration_cycle_cb[!(*index_last_written_cb)] = temp_down_cycle;
		}
	}

	// it a tile of type one was the latest to have been calculated or
	// a tile of type zero has matched the condition above a tile of type one is scheduled to execute next
	// tiles of type one are always executed before the tile of type zero because the execution proceeds from the bottom
	// left to the upper right
	if (((((temp_it_cycle) < (len_query + DB_LEN - 1)-1) && temp_tile_type)||split_down)&&len_query>DIAG_SEGMENT) {
		if(split_down){
			*general_pos+=1;
			pos_to_assign = *general_pos;
		}
		//in order to have a more regular execution we always write on BRAM, but we confirm the index increase only if it was necessary
		write_to_BRAM(
				index_to_write, rowToStart_local_buffer,rowToReach_local_buffer,
				q_block_lb, endQueryIndex_lb, db_block_lb,endDbIndex_lb,
				&tmp_q_block_d_BRAM, &tmp_q_end_d_BRAM, tmp_db_block_d_BRAM, tmp_db_end_d_BRAM,
				temp_down_cycle, 1, split_down, position_antidiagonal_lb, pos_to_assign,
				tile_type_local_buffer, iteration_cycle_lb, temp_H_curr,
				temp_H_prev, temp_E_val, temp_F_val, temp_rowToStart ,temp_rowToReach, len_db, len_query,
				local_buffer_H_curr, local_buffer_H_prev, local_buffer_E_val, local_buffer_F_val,
				dep_H_curr_tmp, dep_H_prev_tmp, dep_E_val_tmp);

		//query update
		local_buffer_query[*index_to_write] = temp_query >> 3;
		local_buffer_query[*index_to_write].range(PASSED_BITS -1,PASSED_BITS -3) = totQuery[tmp_q_block_d_BRAM].range(tmp_q_end_d_BRAM, tmp_q_end_d_BRAM -2);

		//db remains the same
		local_buffer_db[*index_to_write] = temp_db;
		//if we have space in the CB
		if (*elements_in_cb < DIM_BUFFER) {
		write_to_cb(index_last_written_cb, elements_in_cb,
				rowToStart_circular_buffer, rowToReach_circular_buffer,
				q_block_cb, endQueryIndex_cb, db_block_cb,endDbIndex_cb,
				&tmp_q_block_d_CB, &tmp_q_end_d_CB, tmp_db_block_d_CB, tmp_db_end_d_CB,
				temp_down_cycle, 1, split_down, position_antidiagonal_cb, pos_to_assign,
				tile_type_circular_buffer, iteration_cycle_cb, temp_H_curr,
				temp_H_prev, temp_E_val, temp_F_val, temp_rowToStart, temp_rowToReach, len_db, len_query,
				circular_buffer_H_curr, circular_buffer_H_prev, circular_buffer_E_val,
				circular_buffer_F_val, dep_H_curr_tmp, dep_H_prev_tmp,
				dep_E_val_tmp);
			//query update
			circular_buffer_query[*index_last_written_cb] = temp_query >> 3;
			circular_buffer_query[*index_last_written_cb].range(PASSED_BITS -1,PASSED_BITS -3) =
					totQuery[tmp_q_block_d_CB].range(tmp_q_end_d_CB, tmp_q_end_d_CB -2);
			//db remains the same
			circular_buffer_db[*index_last_written_cb] = temp_db;

		} else {

			//if there is no place in the CB we confirm the BRAM write
			*index_to_write +=1;
			if(*index_to_write==DIM_BUFFER_LOCAL) *index_to_write = 0;
			*elements_in_BRAM+=1;

		}

	}

	//go right
	if (temp_it_cycle < (n_db_tiles*DIAG_SEGMENT)-1 && !temp_tile_type) {
		pos_to_assign = 0;
		//always write on BRAM
		temp_db = temp_db<<3;
		if(tmp_db_end_r + 3 == PORT_BITWIDTH){
			tmp_db_block_r += 1;
			tmp_db_end_r  = 2;
		}else{
			tmp_db_end_r += 3;
		}
		ap_uint<3> char_to_add = totDb[tmp_db_block_r].range(tmp_db_end_r, tmp_db_end_r - 2);
		//in order to have a more regular execution we always write on BRAM, but we confirm the index increase only if it was necessary
		write_to_BRAM(	index_to_write, rowToStart_local_buffer,
						rowToReach_local_buffer, q_block_lb,
						endQueryIndex_lb, db_block_lb,endDbIndex_lb, &tmp_q_block_r, &tmp_q_end_r,
						tmp_db_block_r, tmp_db_end_r, temp_right_cycle, 0, 0, position_antidiagonal_lb, pos_to_assign,
						tile_type_local_buffer, iteration_cycle_lb, temp_H_curr,
						temp_H_prev, temp_E_val, temp_F_val, temp_rowToStart, temp_rowToReach, len_db, len_query, local_buffer_H_curr,
						local_buffer_H_prev, local_buffer_E_val, local_buffer_F_val,
						dep_H_curr_tmp, dep_H_prev_tmp, dep_E_val_tmp);
		local_buffer_db[*index_to_write] = temp_db + char_to_add;
		//query remains the same
		local_buffer_query[*index_to_write] = temp_query;


		//if we can we place in the CB
		if (*elements_in_cb < DIM_BUFFER) {
		iteration_cycle_cb[!(*index_last_written_cb)] = temp_right_cycle;
		write_to_cb(index_last_written_cb, elements_in_cb,
					rowToStart_circular_buffer, rowToReach_circular_buffer,
					q_block_cb, endQueryIndex_cb, db_block_cb,endDbIndex_cb,
					&tmp_q_block_r, &tmp_q_end_r, tmp_db_block_r, tmp_db_end_r,
					temp_right_cycle, 0, 0, position_antidiagonal_cb,pos_to_assign,tile_type_circular_buffer,
					iteration_cycle_cb, temp_H_curr, temp_H_prev, temp_E_val,
					temp_F_val, temp_rowToStart, temp_rowToReach, len_db, len_query, circular_buffer_H_curr, circular_buffer_H_prev,
					circular_buffer_E_val, circular_buffer_F_val,
					dep_H_curr_tmp, dep_H_prev_tmp, dep_E_val_tmp);

			circular_buffer_db[*index_last_written_cb] = temp_db + char_to_add;
			//query remains the same
			circular_buffer_query[*index_last_written_cb] = temp_query;

		} else {	
			//if there is no place in the CB we confirm the BRAM write
			*index_to_write +=1;
			if(*index_to_write==DIM_BUFFER_LOCAL) *index_to_write = 0;
			*elements_in_BRAM+=1;
		}
	}
	//write the last calculated dependency for the next tile into the opposite circular buffer cell we have filled
	short newIndex = !(*index_last_written_cb);
	if (!temp_tile_type) {
		H_curr_dep_cb[newIndex] = 0;
		H_prev_dep_cb[newIndex] = 0;
		E_val_dep_cb[newIndex] = 0;
	} else {
		H_curr_dep_cb[newIndex] = temp_H_curr[1];
		H_prev_dep_cb[newIndex] = temp_H_prev[0];
		E_val_dep_cb[newIndex] = temp_E_val[1];

	}
}
	/*
	*	@brief: kernel executes all the previously explained functions in cycle. Everytime it cycles it calls the scoringKernel function to compute a chunk of
	*			antidiagonal, then it continues by writing it into DRAM and after that it calls the setNextValue function in order to have the values ready for the
	*			next cycle. When all the segments have been computed it also selects the maximum value of the kernel and its position and writes it back into DRAM.
	*
	*	@param: query indicates the read stored in DRAM
	*	@param: db indicates the database stored in DRAM
	*	@param:	D indicates the buffer where we store all the segments in DRAM
	*	@param: values stores the length, the number of chuncks of width valid bits we have to retrieve and also the number
				of tiles for both query and db. When the computation is completed it also saves the maximum values
	*
	*/

void kernel(ap_uint<PORT_BITWIDTH> *query, ap_uint<PORT_BITWIDTH> *db,
		ap_uint<PASSED_BITS + 2> D[MAX_N_DIAGONALS], short values[6] //in this we'll store the maximum, it's antidiagonal, and it's offset respectively, we can also use it to store the information
																   //regarding the length of the query, in order to resize our pipeline accordingly
		) {

#pragma HLS INTERFACE m_axi port=query offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=db offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=D offset=slave bundle=gmem2
#pragma HLS INTERFACE m_axi port=values offset=slave bundle=gmem3

#pragma HLS INTERFACE s_axilite port=query bundle=control
#pragma HLS INTERFACE s_axilite port=db bundle=control
#pragma HLS INTERFACE s_axilite port=D bundle=control
#pragma HLS INTERFACE s_axilite port=values bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
#pragma HLS INLINE

	////LOCAL BUFFERS
	ap_uint<VALID_BITS> tot_database_local[LENGTH_BUFFERS_DB_QUERY];
	ap_uint<VALID_BITS> tot_query_local[LENGTH_BUFFERS_DB_QUERY];
	ap_uint<PASSED_BITS> local_buffer_query[DIM_BUFFER_LOCAL];
	ap_uint<PASSED_BITS> local_buffer_db[DIM_BUFFER_LOCAL];
	short local_buffer_H_curr[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1];
	short local_buffer_H_prev[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1];
	short local_buffer_E_val[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1];
	short local_buffer_F_val[DIM_BUFFER_LOCAL][DIAG_SEGMENT + 1];
	bool tile_type_local_buffer[DIM_BUFFER_LOCAL];
	short rowToStart_local_buffer[DIM_BUFFER_LOCAL];
	short rowToReach_local_buffer[DIM_BUFFER_LOCAL];
	short q_block_lb[DIM_BUFFER_LOCAL];
	short db_block_lb[DIM_BUFFER_LOCAL];
	short endQueryIndex_lb[DIM_BUFFER_LOCAL];
	short endDbIndex_lb[DIM_BUFFER_LOCAL];
	short position_antidiagonal_lb[DIM_BUFFER_LOCAL];
	short iteration_cycle_lb[DIM_BUFFER_LOCAL];
	short H_curr_dep_lb[DIM_BUFFER_LOCAL];
	short H_prev_dep_lb[DIM_BUFFER_LOCAL];
	short E_val_dep_lb[DIM_BUFFER_LOCAL];

#pragma HLS RESOURCE variable=tot_database_local core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=tot_query_local core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=local_buffer_query core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=local_buffer_db core=RAM_T2P_BRAM
#pragma HLS ARRAY_PARTITION variable=local_buffer_H_curr block factor=BLOCK_FACTOR_SCORE_VALUES dim=2
#pragma HLS ARRAY_PARTITION variable=local_buffer_H_prev block factor=BLOCK_FACTOR_SCORE_VALUES dim=2
#pragma HLS ARRAY_PARTITION variable=local_buffer_E_val block factor=BLOCK_FACTOR_SCORE_VALUES dim=2
#pragma HLS ARRAY_PARTITION variable=local_buffer_F_val block factor=BLOCK_FACTOR_SCORE_VALUES dim=2
#pragma HLS RESOURCE variable=tile_type_local_buffer core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=rowToStart_local_buffer core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=rowToReach_local_buffer core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=q_block_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=db_block_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=endQueryIndex_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=endDbIndex_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=position_antidiagonal_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=iteration_cycle_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=H_curr_dep_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=H_prev_dep_lb core=RAM_T2P_BRAM
#pragma HLS RESOURCE variable=E_val_dep_lb core=RAM_T2P_BRAM

	//// CIRCULAR BUFFERS
	ap_uint<PASSED_BITS> circular_buffer_query[DIM_BUFFER];
	ap_uint<PASSED_BITS> circular_buffer_db[DIM_BUFFER];
	ap_uint<PASSED_BITS> circular_buffer_directions[DIM_BUFFER];
	short circular_buffer_H_curr[DIM_BUFFER][DIAG_SEGMENT + 1];
	short circular_buffer_H_prev[DIM_BUFFER][DIAG_SEGMENT + 1];
	short circular_buffer_E_val[DIM_BUFFER][DIAG_SEGMENT + 1];
	short circular_buffer_F_val[DIM_BUFFER][DIAG_SEGMENT + 1];
	short rowToStart_circular_buffer[DIM_BUFFER];
	short rowToReach_circular_buffer[DIM_BUFFER];
	bool tile_type_circular_buffer[DIM_BUFFER];
	short q_block_cb[DIM_BUFFER];
	short db_block_cb[DIM_BUFFER];
	short endQueryIndex_cb[DIM_BUFFER];
	short endDbIndex_cb[DIM_BUFFER];
	short H_curr_dep_cb[DIM_BUFFER];
	short H_prev_dep_cb[DIM_BUFFER];
	short E_val_dep_cb[DIM_BUFFER];
	short iteration_cycle_cb[DIM_BUFFER];
	short buffer_max[DIAG_SEGMENT];
	short coord_max_buffer[DIAG_SEGMENT];
	short diag_max_buffer[DIAG_SEGMENT];
	short pos_max_diag_buffer[DIAG_SEGMENT];
	short position_antidiagonal_cb[DIM_BUFFER];

#pragma HLS ARRAY_PARTITION variable=circular_buffer_query complete dim=1
#pragma HLS ARRAY_PARTITION variable=circular_buffer_db complete dim=1
#pragma HLS ARRAY_PARTITION variable=circular_buffer_directions complete dim=1
#pragma HLS ARRAY_PARTITION variable=circular_buffer_H_curr complete dim=2
#pragma HLS ARRAY_PARTITION variable=circular_buffer_H_prev complete dim=2
#pragma HLS ARRAY_PARTITION variable=circular_buffer_E_val complete dim=2
#pragma HLS ARRAY_PARTITION variable=circular_buffer_F_val complete dim=2
#pragma HLS ARRAY_PARTITION variable=rowToStart_circular_buffer complete dim=1
#pragma HLS ARRAY_PARTITION variable=rowToReach_circular_buffer complete dim=1
#pragma HLS ARRAY_PARTITION variable=tile_type_circular_buffer complete dim=1
#pragma HLS ARRAY_PARTITION variable=q_block_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=db_block_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=endQueryIndex_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=endDbIndex_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=H_curr_dep_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=H_prev_dep_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=E_val_dep_cb complete dim=1
#pragma HLS ARRAY_PARTITION variable=iteration_cycle_cb complete dim=1
#pragma HLS ARRAY PARTITION variable=buffer_max complete dim=1
#pragma HLS ARRAY PARTITION variable=coord_max_buffer complete dim=1
#pragma HLS ARRAY PARTITION variable=diag_max_buffer complete dim=1
#pragma HLS ARRAY PARTITION variable=pos_max_diag_buffer complete dim=1
#pragma HLS ARRAY PARTITION variable=position_antidiagonal_cb complete dim=1

	//read information from the values array
	short n_query_tiles = values[0];
	short n_db_tiles = values[1];
	short db_chunks_to_load = values[2];
	short q_chunks_to_load = values[3];
	short len_db = values[4];
	short len_query = values[5];
	//calculate the number of diagonal segments we need to execute
	short max_horizontal_length = n_db_tiles*DIAG_SEGMENT;
	short max_vertical_length = n_query_tiles*DIAG_SEGMENT;
	int length_tiles_one = max_vertical_length*n_db_tiles;
	int max_current_diags = max_horizontal_length+length_tiles_one;

	//initialize indexed for the first execution
	bool index_last_written_cb = 0;
	bool to_read_from_cb = 0;
	short index_to_write = 0;
	short last_read_from_memory = 0;
	short elements_in_cb = 1;
	short elements_in_BRAM = 0;
	short end_query = PASSED_BITS-1;
	if(len_query<DIAG_SEGMENT) end_query=len_query*3-1;
	
	//LOAD DB AND QUERY IN BRAM
	load_query_loop: for (int i = 0; i < q_chunks_to_load; i++) {
#pragma HLS PIPELINE
		tot_query_local[i] = query[i];
	}

	load_db_loop: for (int i = 0; i < db_chunks_to_load; i++) {
#pragma HLS PIPELINE
		tot_database_local[i] = db[i];
	}

	//set the first pieces of query and database to align
	circular_buffer_query[to_read_from_cb] = tot_query_local[0].range(end_query,0);
	q_block_cb[to_read_from_cb] = 0;

	circular_buffer_db[to_read_from_cb] = tot_database_local[0].range(2,0);
	db_block_cb[to_read_from_cb] = 0;

	//initialize indexes to keep track of the position to read from in the local query and db buffers
	tile_type_circular_buffer[to_read_from_cb] = 0;
	rowToStart_circular_buffer[to_read_from_cb] = 0;
	rowToReach_circular_buffer[to_read_from_cb] = 0;
	endDbIndex_cb[to_read_from_cb] = 2;
	endQueryIndex_cb[to_read_from_cb] = end_query;
	iteration_cycle_cb[to_read_from_cb] = 0;
	position_antidiagonal_cb[to_read_from_cb]=0;
	int i = 0;
	ap_uint<PASSED_BITS> diag_k;
	ap_uint<PORT_BITWIDTH> diagonal;
	short limit = 0;
	short diag_index = 0;
	short exec_d = 0;
	short index_shift = 1;
	short general_pos = 0;

	//set all dependencies to zero for the first execution
	initialize_arrays_loop:	for (int j = 0; j < DIAG_SEGMENT + 1; j++) {
#pragma HLS UNROLL
		circular_buffer_H_curr[0][j] = 0;
		circular_buffer_H_prev[0][j] = 0;
		circular_buffer_E_val[0][j] = 0;
		circular_buffer_F_val[0][j] = 0;
	}
	set_ax_zero_loop: for (int j = 0; j < DIAG_SEGMENT; j++){
#pragma HLS UNROLL
		buffer_max[j]=0;
	}

	query_loop: for(int m = 0; m<MAX_N_DIAGONALS; m++){
#pragma HLS PIPELINE
		//execute only for the cycles needed, the number of cycles is calculated above
		if(m<max_current_diags){
			short q_block_BRAM = q_block_lb[last_read_from_memory];
			short endQuery_index_BRAM = endQueryIndex_lb[last_read_from_memory];
			short db_block_BRAM = db_block_lb[last_read_from_memory];
			short endDb_index_BRAM = endDbIndex_lb[last_read_from_memory];
			short rowToStart_BRAM = rowToStart_local_buffer[last_read_from_memory];
			short rowToReach_BRAM = rowToReach_local_buffer[last_read_from_memory];
			ap_uint<PASSED_BITS+2>  query_to_pass_BRAM = local_buffer_query[last_read_from_memory];
			ap_uint<PASSED_BITS+2>	db_to_pass_BRAM = local_buffer_db[last_read_from_memory];

			bool tile_type_BRAM=tile_type_local_buffer[last_read_from_memory];
			short H_curr_BRAM[DIAG_SEGMENT+1];
			short F_val_BRAM[DIAG_SEGMENT+1];
			short E_val_BRAM[DIAG_SEGMENT+1];
			short H_prev_BRAM[DIAG_SEGMENT+1];
			for (int n = 0; n < DIAG_SEGMENT + 1; n++) {
				H_curr_BRAM[n] = local_buffer_H_curr[last_read_from_memory][n];
				H_prev_BRAM[n] = local_buffer_H_prev[last_read_from_memory][n];
				E_val_BRAM[n] = local_buffer_E_val[last_read_from_memory][n];
				F_val_BRAM[n] =local_buffer_F_val[last_read_from_memory][n];
			}
			short H_curr_dep_BRAM = H_curr_dep_lb[last_read_from_memory];
			short H_prev_dep_BRAM = H_prev_dep_lb[last_read_from_memory];
			short E_val_dep_BRAM = E_val_dep_lb[last_read_from_memory];
			short iteration_cycle_BRAM = iteration_cycle_lb[last_read_from_memory];
			short position_antidiagonal_BRAM = position_antidiagonal_lb[last_read_from_memory];
			
			//the scoring kernel reads always from the circular buffer, in order to have 
			//everytihing in logic and be able to calculate an entire antidiagonal in parallel
			//the cell from which it has to read it's changed at the end of every cycle
			scoringKernel(&(circular_buffer_query[to_read_from_cb]),
					&(circular_buffer_db[to_read_from_cb]), &(diag_k),
					circular_buffer_H_curr[to_read_from_cb],
					circular_buffer_H_prev[to_read_from_cb],
					circular_buffer_F_val[to_read_from_cb],
					circular_buffer_E_val[to_read_from_cb], buffer_max,
					coord_max_buffer, diag_max_buffer, iteration_cycle_cb[to_read_from_cb],
					position_antidiagonal_cb[to_read_from_cb],pos_max_diag_buffer,
					rowToReach_circular_buffer[to_read_from_cb],
					rowToStart_circular_buffer[to_read_from_cb]);


			D[i] = diag_k;//write the calculated antidiagonal

			setNextValues(tot_query_local, tot_database_local,
					q_block_BRAM, endQuery_index_BRAM, db_block_BRAM,
					endDb_index_BRAM,
					rowToStart_BRAM,
					rowToReach_BRAM,
					query_to_pass_BRAM,
					db_to_pass_BRAM,
					tile_type_BRAM,
					F_val_BRAM,
					E_val_BRAM,
					H_curr_BRAM,
					H_prev_BRAM,
					H_curr_dep_BRAM,
					H_prev_dep_BRAM,
					E_val_dep_BRAM,
					iteration_cycle_BRAM,
					position_antidiagonal_BRAM,
					q_block_cb, endQueryIndex_cb,
					db_block_cb, endDbIndex_cb, q_block_lb,
					endQueryIndex_lb, db_block_lb, endDbIndex_lb,
					tile_type_local_buffer, tile_type_circular_buffer, position_antidiagonal_lb,
					position_antidiagonal_cb,
					&index_last_written_cb, &index_to_write, &last_read_from_memory,
					to_read_from_cb, circular_buffer_query, circular_buffer_db,
					circular_buffer_H_curr, circular_buffer_H_prev,
					circular_buffer_F_val, circular_buffer_E_val,
					rowToReach_circular_buffer, rowToStart_circular_buffer,
					local_buffer_query, local_buffer_db, local_buffer_H_curr,
					local_buffer_H_prev, local_buffer_F_val, local_buffer_E_val,
					rowToReach_local_buffer,
					rowToStart_local_buffer,
					H_curr_dep_cb, H_prev_dep_cb, E_val_dep_cb, H_curr_dep_lb,
					H_prev_dep_lb, E_val_dep_lb, &elements_in_cb, &elements_in_BRAM,&general_pos,
					iteration_cycle_cb, iteration_cycle_lb, n_query_tiles,
					n_db_tiles, len_db, len_query);

			//tell to the circular buffer which cell it has to read
			to_read_from_cb = !(to_read_from_cb);
			i++;
		}


	}

	short max_kernel = -1;
	short max_pos = 73;
	short i_max;

	//calculate the maximum value and it's position
	max_loop:for (int j = 0; j < DIAG_SEGMENT; j++) {
#pragma HLS PIPELINE
		if (buffer_max[j] > max_kernel&&pos_max_diag_buffer[j]<max_pos) {
			max_kernel = buffer_max[j];
			i_max = j;
		}
	}
	values[2] = max_kernel;
	values[3] = pos_max_diag_buffer[i_max]*DIAG_SEGMENT-coord_max_buffer[i_max]-1;
	values[4] = diag_max_buffer[i_max]+coord_max_buffer[i_max]-DIAG_SEGMENT+1;
	std::cout << "\nNum_cycles : " << i << "\n";
}

//}

