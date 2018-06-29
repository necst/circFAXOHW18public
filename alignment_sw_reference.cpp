//#include <ap_int.h>
//#include <ap_utils.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <bitset>
#include <math.h>
#include <time.h>
#include <ctime>

const short STARTGAP = 3;
const short EXTENDGAP = 1;
const short MATCH = 2;
const short MISMATCH = -2;
const short BAD_MISMATCH = -3;
const short N_MISMATCH = -1;
const short N_BASES = 5;

const short DB_LEN = 3024;
const short QUERY_LEN = 3024;
const short QUERY_SEGMENT_LENGTH = QUERY_LEN * 3;
const short DIAG_SEGMENT = QUERY_LEN;
const short DIAG_KERNEL = 42;
const short MAX_PASSED_LENGTH = 72*DIAG_KERNEL;
//const short DB_SEGMENT_LENGTH = DB_LEN * 3;
const int MAX_N_DIAGONALS = MAX_PASSED_LENGTH*73;


const short noDirection = 0;
const short west = 1;
const short north =2;
const short northWest = 4;
const short start_GAP = 3;
const short extend_GAP = 1;


const short PORT_BITWIDTH = 512;
const short VALID_BITS = PORT_BITWIDTH-2;
const short MAX_DIM = (MAX_PASSED_LENGTH*3+VALID_BITS-1)/VALID_BITS;
const short PASSED_BITS = 126;
const short PARTITION_HALF = DIAG_SEGMENT/2;
const int DIM_MAX_DIAG = MAX_N_DIAGONALS*(PASSED_BITS+2);

const short N_DIAGONALS = DB_LEN + QUERY_LEN - 1;

const short W[N_BASES][N_BASES]={
		{MATCH,BAD_MISMATCH,MISMATCH,MISMATCH,N_MISMATCH},
		{BAD_MISMATCH,MATCH,MISMATCH,MISMATCH,N_MISMATCH},
		{MISMATCH,MISMATCH,MATCH,BAD_MISMATCH,N_MISMATCH},
		{MISMATCH,MISMATCH,BAD_MISMATCH,MATCH,N_MISMATCH},
		{N_MISMATCH,N_MISMATCH,N_MISMATCH,N_MISMATCH,MATCH}
};


int main(int argc ,char *argv[]){

	char *db_k = (char*)malloc(DB_LEN*8);
	char *query_k = (char*)malloc(QUERY_LEN*8);

	srand(time(NULL));
	

	short index_l=0;
	short q_index_loc=0;
	short db_index_loc=0;
	for(int i = 0; i < DB_LEN; i++){
		char r=rand()%5;
		db_k[i]=r;
	}

	index_l=0;
	for(int i = 0; i < QUERY_LEN; i++){
		char r=rand()%5;
		query_k[i]=r;
		
	}

	////// software kernel implementation

	clock_t tstart = clock();
	int **HLOC = (int **)malloc(QUERY_LEN*sizeof(int*));
	for(int i = 0; i < QUERY_LEN; i++) HLOC[i] = (int *)malloc(DB_LEN*sizeof(int));
	int **FLOC = (int **)malloc(QUERY_LEN*sizeof(int*));
	for(int i = 0; i < QUERY_LEN; i++) FLOC[i] = (int *)malloc(DB_LEN*sizeof(int));
	int **ELOC = (int **)malloc(QUERY_LEN*sizeof(int*));
	for(int i = 0; i < QUERY_LEN; i++) ELOC[i] = (int *)malloc(DB_LEN*sizeof(int));
	int **directionMatrix = (int **)malloc(QUERY_LEN*sizeof(int*));
	for(int i = 0; i < QUERY_LEN; i++) directionMatrix[i] = (int *)malloc(DB_LEN*sizeof(int));



	short max = -1;
	short colMax = 0;
	short rowMax = 0;
	short numMax = 0;
	short db_index = 0;
	short query_index = 0;
	short i_local = 0;
	short j_local = 0;

	for(int i = 0; i < DB_LEN; i++){//column

		short dbL = db_k[i];
	
		for(int j = 0; j < QUERY_LEN; j++){//rows

			char q = query_k[j];
			
			if(j==0){//first row
				if(i==0){//first column
					HLOC[j][i]=(W[q][dbL]>0)?W[q][dbL]:0;
					ELOC[j][i]=0;
					FLOC[j][i]=0;
					directionMatrix[j][i]=noDirection;//
					if(HLOC[j][i] > 0) directionMatrix[j][i]+=northWest;
					if(HLOC[j][i] > max){
						max = HLOC[j][i];
						colMax = i;
						rowMax = j;
						numMax = 1;
					}
					else if(HLOC[j][i] == max) numMax++;
				}
				FLOC[j][i]=-extend_GAP;
				ELOC[j][i]=((ELOC[j][i-1]-extend_GAP)>(HLOC[j][i-1]-start_GAP))?(ELOC[j][i-1]-extend_GAP):(HLOC[j][i-1]-start_GAP);
				short score=(W[q][dbL]>=ELOC[j][i])?W[q][dbL]:ELOC[j][i];
				directionMatrix[j][i]=noDirection;
				score =(score>=0)?score:0;
				if (score==W[q][dbL]) directionMatrix[j][i]+=northWest;
				if (score==ELOC[j][i]) directionMatrix[j][i]+=west;
				if (score<0) directionMatrix[j][i]=noDirection;
				HLOC[j][i]=score;
				if(HLOC[j][i] > max){
					max = HLOC[j][i];
					colMax = i;
					rowMax = j;
					numMax = 1;
				}else if(HLOC[j][i] == max) numMax++;
			}
			else if(i==0&&j>0){//first column
				FLOC[j][i]=((FLOC[j-1][i]-extend_GAP)>(HLOC[j-1][i]-start_GAP))?(FLOC[j-1][i]-extend_GAP):(HLOC[j-1][i]-start_GAP);
				short score=(W[q][dbL]>=FLOC[j][i])?W[q][dbL]:FLOC[j][i];
				score =(score>=0)?score:0;
				ELOC[j][i]=-extend_GAP;
				directionMatrix[j][i] = noDirection;
				if (score==W[q][dbL]) directionMatrix[j][i]+=northWest;
				if (score==FLOC[j][i]) directionMatrix[j][i]+=north;
				HLOC[j][i]=score;
				if(HLOC[j][i] > max){
					max = HLOC[j][i];
					colMax = i;
					rowMax = j;
					numMax = 1;
				}
				else if(HLOC[j][i] == max) numMax++;
			}
			else if(j>0&&i>0){
				short Flocal = ((FLOC[j-1][i]-extend_GAP)>(HLOC[j-1][i]-start_GAP))?FLOC[j-1][i]-extend_GAP:HLOC[j-1][i]-start_GAP;
				short Elocal = ((ELOC[j][i-1]-extend_GAP)>(HLOC[j][i-1]-start_GAP))?ELOC[j][i-1]-extend_GAP:HLOC[j][i-1]-start_GAP;
				short N_W_score = HLOC[j-1][i-1]+W[q][dbL];
				short Hlocal = (Flocal>Elocal)?Flocal:Elocal;
				Hlocal = (Hlocal>N_W_score)?Hlocal:N_W_score;
				Hlocal = (Hlocal>0)?Hlocal:0;

				directionMatrix[j][i]=noDirection;
				if (Flocal == Hlocal) directionMatrix[j][i]+= north;
				if (Elocal == Hlocal) directionMatrix[j][i]+= west;
				if (N_W_score == Hlocal) directionMatrix[j][i]+=northWest;
				FLOC[j][i]=Flocal;
				ELOC[j][i]=Elocal;
				HLOC[j][i]=Hlocal;
				if(Hlocal>max){
					max =Hlocal;
					colMax = i;
					rowMax = j;
					numMax = 1;

				}
			}

		}

	}
	std::cout<<"TIME: "<<(double) (clock()-tstart)/CLOCKS_PER_SEC<<std::endl;
	return 0;
}

