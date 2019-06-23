#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>

#include "DECODE.h"
#include "mmio.h"


void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: DECODE\n ");
    fprintf (stream,
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
		 "\t-i	--input FILENAME\tfullpath of the input expression file in mm format (MANDATORY)\n"
		 "\t-t	--type {mm|table|csv}\tformat of input file (MANDATORY)\n"		 
		 "\t-g	--group COLFILE\t\tFile with each row describing columns of interest in the expression matrix\n"		 
		 "\t-n	--null COLFILE\t\tFile with each row describing null columns in the expression matrix (default=all columns minus \"group\")\n"		 
		 "\t-r	--rows ROWFILE\t\tRows of matrix to perform differential analysis (default = all rows)\n"		 
		 "\t-c	--sample_no INT\t\tNumber of sub-samples (default = 1000)\n"
		 "\t-s	--sample_size INT\t\tSub-sample size (default = 100)\n"
		 "\t-p	--thread_no INT\t\tNumber of threads (default = -1)\n"
		 "\t-o	--output FILENAME\t\tOutput file (default = logpvals.txt)\n"
	);	 
		
    exit (exit_code);
}


int main(int argc, char ** argv) {
  int next_option;
  const char *const short_options = "hi:t:g:n:p:o:c:s:";
  const struct option long_options[] = {
		{"help",     0, NULL, 'h'},
		{"input",  	 1, NULL, 'i'},
		{"type",  	 1, NULL, 't'},
		{"group", 	 1, NULL, 'g'},
		{"null", 	 1, NULL, 'n'},
		{"sample_no", 	 1, NULL, 'c'},
		{"sample_size", 	 1, NULL, 's'},
		{"thread_no", 	 1, NULL, 'p'},
		{"output", 	 1, NULL, 'o'},
		{NULL,       0, NULL,  0 }		
	};


	char out_fname[1024] = "logpvals.txt", input_file[1024] = "";
		
	int type = -1; // -1: unknown, 0: mm, 1: table, 2: csv
	uvec group, null;
	int thread_no = -1;
	int rand_sample_no = 10000, rand_sample_size = 100;
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'i':
				strcpy(input_file, optarg);
	    		break;

			case 't':
				if(!strcmp(optarg, "mm"))
					type = 0;
				else if(!strcmp(optarg, "table"))
					type = 1;				
				else if(!strcmp(optarg, "csv"))
					type = 2;
				else {
					fprintf(stderr, "unknown file format\n");
					print_usage(stderr, -1);					
					return -1;
				}
	    		break;

	    	case 'g':
				group = read_uvec(optarg);
	    		break;

	    	case 'n':
				null = read_uvec(optarg);
	    		break;
	    		
			case 'p':
				thread_no = atoi(optarg);
	    		break;			

			case 'c':
				rand_sample_no = atof(optarg);
	    		break;			

			case 's':
				rand_sample_size = atof(optarg);
	    		break;			

			case 'o':
				strcpy(out_fname, optarg);
	    		break;	    		
		}
    } while (next_option != -1);


	sp_mat expression;
	if(!strcmp(input_file, "") || type == -1) {
		fprintf(stderr, "full path to the expression file and its type are mandatory arguments.\n");
		print_usage(stderr, -1);
		return -1;
	}
	if(group.n_elem == 0) {
		fprintf(stderr, "primary column group argument is mandatory (-g | --group).\n");
		print_usage(stderr, -1);
		return -1;
	}
	
	switch(type) {
		case 0:
			expression = read_from_mm(input_file);
			break;
			
		case 1:
			expression = read_from_table(input_file);
			break;
			
		case 2:
			expression = read_from_csv(input_file);
			break;
					
	}

	vec logPvals;
	if(null.n_elem == 0) {
		logPvals = AssessFeatures(expression, group, thread_no, rand_sample_no, rand_sample_size);
	}
	else {
		logPvals = AssessFeatures_betweenGroups(expression, group, null, thread_no, rand_sample_no, rand_sample_size);
	}

	// Export!
	logPvals.save(out_fname, arma::raw_ascii);
	return 0;
}
