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

	int vector_size;
	
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
	
	m = (int ) mxGetM(prhs[1]);
	n = (int ) mxGetN(prhs[1]);	
	if(min(m, n) != 1) {
		mexErrMsgTxt("Cols parameter should be a one dimensional vector.");
	}	
	vector_size = max(m, n);
	uvec cols(vector_size);
	if(mxIsSparse(prhs[1])) {
		mat temp = mat(armaGetSparseMatrix(prhs[1]));  
		for (i = 0; i < vector_size; i++) {
			cols(i) = (uword)(temp[i]-1);
		}		
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 
		for (i = 0; i < vector_size; i++) {
			cols(i) = (uword)(ptr[i]-1);
		}		
	}

	int rand_sample_no = 100, rand_sample_size = 1000;
	if(nrhs > 2) {
		rand_sample_size = mxGetScalar(prhs[2]);
	}							
	if(nrhs > 3) {
		rand_sample_no = mxGetScalar(prhs[3]);
	}
	
	int thread_no = -1;
	if(nrhs > 4) {
		thread_no = mxGetScalar(prhs[4]);
	}					

	vec logPvals = AssessFeatures(A, cols, rand_sample_size, rand_sample_no, thread_no);

	// Return output
	vec v = logPvals;
	plhs[0] = mxCreateDoubleMatrix(v.n_rows, 1, mxREAL);
	memcpy(mxGetPr(plhs[0]), v.memptr(), v.n_elem * sizeof(double)); 	
}
