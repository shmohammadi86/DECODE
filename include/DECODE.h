#ifndef DECODE_H
#define DECODE_H

//#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG

#include <omp.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>

#include "armadillo"
using namespace arma;
using namespace std;

#define STATS_USE_ARMA
#define STATS_USE_OPENMP
#include <stats.hpp>



sp_mat read_from_mm(char *path);
sp_mat read_from_table(char *path);
sp_mat read_from_csv(char *path);
uvec read_uvec(char *path);

vec AssessFeatures(sp_mat A, uvec cols, int rand_sample_size, int rand_sample_no, int thread_no);
vec AssessFeatures_betweenGroups(sp_mat A, uvec cols, uvec cols_null, int rand_sample_size, int rand_sample_no, int thread_no);
vec ProfileModule(sp_mat A, uvec M, int rand_sample_size, int rand_sample_no, int thread_no);

#endif
