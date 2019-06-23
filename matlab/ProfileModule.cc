#include "armaMex.hpp"
#include "DECODE.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;
	

	// Check type of input.
	if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) )
	mexErrMsgTxt("Input must me of type double.");

	// Check if input is real.
	if ( (mxIsComplex(prhs[0])) )
	mexErrMsgTxt("Input must be real.");

	
	sp_mat A;
	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);
	if(mxIsSparse(prhs[0])) {
		A = armaGetSparseMatrix(prhs[0]);  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[0])); 
		
		A = sp_mat(mat(ptr, m, n, true, true));		
	}
	

	int m2 = (int ) mxGetM(prhs[1]);
	int n2 = (int ) mxGetN(prhs[1]);
	if(min(m2, n2) != 1) {
		mexErrMsgTxt("Rows parameter should be a one dimensional vector.");
	}	
	int vector_size = max(m2, n2);
	uvec rows(vector_size);
	if(mxIsSparse(prhs[1])) {
		mat temp = mat(armaGetSparseMatrix(prhs[1]));  
		for (i = 0; i < vector_size; i++) {
			rows(i) = (uword)(temp[i]-1);
		}		
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 
		for (i = 0; i < vector_size; i++) {
			rows(i) = (uword)(ptr[i]-1);
		}		
	}


	int rand_sample_no = 10000, rand_sample_size = 100;
	if(nrhs > 3) {
		rand_sample_size = mxGetScalar(prhs[3]);
	}							
	if(nrhs > 4) {
		rand_sample_no = mxGetScalar(prhs[4]);
	}
	
	int thread_no = -1;
	if(nrhs > 5) {
		thread_no = mxGetScalar(prhs[5]);
	}	
	
	vec profile = ProfileModule(A, rows, rand_sample_size, rand_sample_no, thread_no);


	plhs[0] = mxCreateDoubleMatrix(profile.n_elem, 1, mxREAL);
	memcpy(mxGetPr(plhs[0]), profile.memptr(), profile.n_elem * sizeof(double)); 	
	
}
