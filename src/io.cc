#include <DECODE.h>

#include <mmio.h>

#define MAX_LINE_SIZE 10000000

sp_mat read_from_mm(char *path) {
	printf("Reading %s ... ", path); fflush(stdout);

    int M, N, nnz;   
	int i;

    int ret_code;
    MM_typecode matcode;
    FILE *f;

	FILE *fd;
	if ((fd = fopen(path, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s\n", path);
		exit(1);
    }

    if (mm_read_banner(fd, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(fd, &M, &N, &nnz)) !=0) {
		fprintf(stderr, "Cannot read matrix banner (header) successfully.\n");
        exit(1);
	}

    /* reseve memory for matrices */

	umat locations(2, nnz);
	vec values(nnz);
	
	double V;
	int I, J;
    for (i = 0; i < nnz; i++)
    {
        fscanf(fd, "%d %d %lf\n", &I, &J, &V);
        locations(0, i) = I-1;
        locations(1, i) = J-1;
        values(i) = V;
    }
    fclose(fd);
    
	sp_mat res(locations, values, (uword)M, (uword)N);
	
	return res;
}

const vector<string> chopChop(string s, const char& c)
{
	string buff{""};
	vector<string> v;
	
	if(s[0] == c) { // remove leading tab
		s.erase(s.begin());
	}
	
	for(auto n:s)
	{		
		if(n != c) buff+=n; else
		if(n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if(buff != "") {v.push_back(buff);}
	
	return v;
}


sp_mat read_from_table(char * path) {
	printf("Reading %s ... ", path); fflush(stdout);
	
	FILE *fd = fopen(path, "r");
	if (fd == NULL) {
		fprintf(stderr, "Can't open file %s\n", path);
		return NULL;
	}
	
	char *buffer = (char *)malloc(MAX_LINE_SIZE);
	if(buffer == NULL) {
		fprintf(stderr, "Can't allocate memory for buffer\n");
		return NULL;
	}

	if(fgets(buffer, MAX_LINE_SIZE, fd) == NULL){
		fprintf(stderr, "%s is empty\n", path);
		return NULL;
	}

	buffer[strlen(buffer)-1] = 0;
	vector<string> sample_ids = chopChop(string(buffer), '\t');
		
	int sample_no = sample_ids.size();
	
	
	// Read all lines to find total number of genes ... such a waste of time!!
	char gene_name[1024];
	vector<string> gene_names;
	while(!feof(fd)) {
		if(fgets(buffer, MAX_LINE_SIZE, fd) == NULL)
			break;
		
		sscanf(buffer, "%1024s", gene_name);	
		gene_names.push_back(string(gene_name));
	}
	fclose(fd);

	
	int line_no = gene_names.size();
	int gene_no = gene_names.size();
	sp_mat expression(gene_no, sample_no);
	
	// Now read the actual expression matrix
	fd = fopen(path, "r");
	fgets(buffer, MAX_LINE_SIZE, fd); // Read and skip the header
		
	register int i, j;
	for (i = 0; i < line_no; i++) {
		fscanf(fd, "%s", gene_name);
		vec row(sample_no);			
		for(j = 0; j < sample_no; j++) {			
			fscanf(fd, "%lf", &row[j]);
		}
		
		uvec idx = find(row);
		for(j = 0; j < idx.n_elem; j++) {
			expression(i, idx(j)) = row(idx(j));
		}
	}
	fclose(fd);
	
	
	printf("done\n");
	free(buffer);
	return expression;
}


sp_mat read_from_csv(char * path) {
	printf("Reading %s ... ", path); fflush(stdout);
	
	FILE *fd = fopen(path, "r");
	if (fd == NULL) {
		fprintf(stderr, "Can't open file %s\n", path);
		return NULL;
	}
	
	char *buffer = (char *)malloc(MAX_LINE_SIZE);
	if(buffer == NULL) {
		fprintf(stderr, "Can't allocate memory for buffer\n");
		return NULL;
	}

	if(fgets(buffer, MAX_LINE_SIZE, fd) == NULL){
		fprintf(stderr, "%s is empty\n", path);
		return NULL;
	}

	buffer[strlen(buffer)-1] = 0;
	vector<string> first_row = chopChop(string(buffer), ',');
		
	int sample_no = first_row.size();
	
	
	// Read all lines to find total number of genes ... such a waste of time!!
	char gene_name[1024];
	vector<string> gene_names;
	int line_no = 0;
	while(!feof(fd)) {
		if(fgets(buffer, MAX_LINE_SIZE, fd) == NULL)
			break;
		line_no ++;
	}
	fclose(fd);

	
	sp_mat expression(line_no, sample_no);
	
	// Now read the actual expression matrix
	fd = fopen(path, "r");
	fgets(buffer, MAX_LINE_SIZE, fd); // Read and skip the header
		
	register int i, j;
	for (i = 0; i < line_no; i++) {
		vec row(sample_no);			
		for(j = 0; j < sample_no; j++) {			
			fscanf(fd, "%lf", &row[j]);
		}
		
		uvec idx = find(row);
		for(j = 0; j < idx.n_elem; j++) {
			expression(i, idx(j)) = row(idx(j));
		}
	}
	fclose(fd);
	
	
	printf("done\n");
	free(buffer);
	return expression;
}


uvec read_uvec(char * path) {
	FILE *fd = fopen(path, "r");
	vector <uword> vals;
	uword val;
	while(!feof(fd)) {
		fscanf(fd, "%ud", &val);
		vals.push_back(val);
	}	
	fclose(fd);
	
	uvec vals_uvec(vals.size());
	for (register int i = 0; i < vals.size(); i++) {
		vals_uvec(i) = vals[i];
	}
}




