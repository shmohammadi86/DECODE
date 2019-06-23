#include <DECODE.h>

vec AssessFeatures_betweenGroups(sp_mat A, uvec cols, uvec cols_null, int rand_sample_size = 10000, int rand_sample_no=100, int thread_no=8) {
	rand_sample_size = std::min(rand_sample_size, (int)(2*cols.n_elem/3));
	
	// Binarize the profile
	A = spones(A);

	vec row_p = zeros(A.n_rows, 1);
	for(int i = 0; i < cols_null.n_elem; i++) {
		row_p += vec(A.col(cols_null(i)));
	}
	row_p = (row_p + 1) / (cols_null.n_elem + 1);

	
	vec col_p = vec(trans(sum(A)));
	vec col_w = col_p / mean(col_p);

	
	// Estimate estimated number of counts using sub-sampling (with replacement)
	printf("Sub-sampling (# rand sample = %d, sample size = %d, threads = %d) ... ", rand_sample_no, rand_sample_size, thread_no); fflush(stdout);	

	mat subsample_delta(A.n_rows, rand_sample_no);				
	umat rand_samples = conv_to<umat>::from(stats::runif<arma::mat>(rand_sample_no, rand_sample_size, 0.0, (double)cols.n_elem-1, 0));

	mat Bounds = zeros(A.n_rows, rand_sample_no);
	mat X =  zeros(A.n_rows, rand_sample_no);
	#pragma omp parallel for num_threads(thread_no)
	for( int k = 0; k < rand_sample_no; k++) {
		double sample_weight = sum(col_w(cols(rand_samples.row(k))));
		vec Exp = row_p*sample_weight;
		
		vec Obs = zeros(A.n_rows, 1);
		for (int i = 0; i < rand_sample_size; i++) {
			Obs += vec(A.col(cols(rand_samples(k, i))));				
		}		
		
		vec Delta = (Obs / Exp)	- 1;
		Delta.replace(datum::nan, 0);
		Delta.replace(datum::inf, 0);
				
		X.col(k) = Delta;

		vec ub = ( ( square(Delta) / (2 + Delta) ) % (Exp) );
		
		uvec filter_idx = find(Delta <= 0);
		ub(filter_idx) = zeros(filter_idx.n_elem);
		
		Bounds.col(k) = ub;
		
	} 	
	printf("done\n");	

	
	// The harmonic mean p-value for combining dependent tests (2019)

	// https://en.wikipedia.org/wiki/LogSumExp
	vec c = max(Bounds, 1);
	mat C = c * ones(1, rand_sample_no);
	vec log_sum_exp = arma::log(sum(exp(Bounds-C), 1)) + c;
	
	vec meta_logpvals = log_sum_exp - std::log(rand_sample_no);
	meta_logpvals = meta_logpvals - log(meta_logpvals.n_elem);
	meta_logpvals.transform( [](double val) { return (val < 0? 0:val); } );

	return(meta_logpvals);
}

vec AssessFeatures(sp_mat A, uvec cols, int rand_sample_size = 10000, int rand_sample_no=100, int thread_no=8) {
	
	vec null_col_mask = ones(A.n_cols, 1);
	null_col_mask(cols) = ones(cols.n_elem);
	uvec cols_null = find(null_col_mask != 0);
	
	vec stats = AssessFeatures_betweenGroups(A, cols, cols_null, rand_sample_size, rand_sample_no, thread_no);
	
	return stats;
}

vec ProfileModule(sp_mat A, uvec M, int rand_sample_size = 10000, int rand_sample_no=100, int thread_no=8) {	
	vec stats = AssessFeatures(A.t(), M, rand_sample_size, rand_sample_no, thread_no);
	
	return stats;
}
