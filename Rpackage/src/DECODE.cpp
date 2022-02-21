#include <DECODE.h>
#include <RcppArmadillo.h>

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec AssessFeatures(sp_mat A_mat, IntegerVector columns, int rand_sample_size = 1000, int rand_sample_no = 100, int thread_no = -1) {
	uvec cols_uvec(columns.size());
	for (int i = 0; i < cols_uvec.n_elem; i++) {
		cols_uvec(i) = (uword)(columns(i)-1);
	}	  

	vec logpvals = AssessFeatures(A_mat, cols_uvec, rand_sample_size, rand_sample_no, thread_no);
				
	return logpvals;			
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec AssessFeatures_betweenGroups(sp_mat A_mat, IntegerVector columns_A, IntegerVector columns_B, int rand_sample_size = 1000, int rand_sample_no = 100, int thread_no = -1) {
	int i;
	
	uvec cols_A_uvec(columns_A.size());
	for (i = 0; i < cols_A_uvec.n_elem; i++) {
		cols_A_uvec(i) = (uword)(columns_A(i)-1);
	}	  
	
	uvec cols_B_uvec(columns_B.size());
	for (i = 0; i < cols_B_uvec.n_elem; i++) {
		cols_B_uvec(i) = (uword)(columns_B(i)-1);
	}	  

	vec logpvals = AssessFeatures_betweenGroups(A_mat, cols_A_uvec, cols_B_uvec, rand_sample_size, rand_sample_no, thread_no);
	
	return(logpvals);
		
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec ProfileModule(sp_mat A_mat, IntegerVector rows, int rand_sample_size = 1000, int rand_sample_no = 100, int thread_no = -1) {
	uvec rows_uvec(rows.size());
	for (int i = 0; i < rows_uvec.n_elem; i++) {
		rows_uvec(i) = (uword)(rows(i)-1);
	}	  

	vec profile = ProfileModule(A_mat, rows_uvec, rand_sample_size, rand_sample_no, thread_no);
	
	return profile;
}
