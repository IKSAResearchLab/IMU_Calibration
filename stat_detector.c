#define _DEFAULT_SOURCE
#include "stat_detector.h"

int save_matrix2csv(Fmatrix_ptr data, char *filename, const char **header)
{ /* 
	Saves a matrix into a csv file name.

  */
	
	/* open file and put the header */
	
	FILE *fp = fopen(filename, "w");
	if (fp == NULL) {
		fprintf(stderr, "Error in save_csv(): can't open the file for writing.");
		exit(1);
	}
	
	for (unsigned long int i = 0; i < data->cols; i ++ ) {
		if (header != NULL) {
			fprintf(fp, "%s", header[i]);
		} else {
			fprintf(fp, "c%lu", i);
		}
		if (i != data->cols - 1) {
			fprintf(fp, " ");
		}
	}
	fprintf(fp, "\n");

	for (unsigned long int i = 0; i < data->rows; i++) {
		for (unsigned long int j = 0; j < data->cols; j++) {
			fprintf(fp, "%f", data->data[i][j]);
			if (j == data->cols - 1) {
				fprintf(fp, "\n");
			} else {
				fprintf(fp, " ");
			}
		}
		
	}

	fclose(fp);
	return(1);
}

int read_csv(char *file_name, Fmatrix_ptr *data) {
	/* 
		read data from csv and store it in the data matrix
		The structure of the CSV file:
			ustime gx gy gz ax ay az
	*/

	FILE *fp = fopen(file_name, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error in read_csv(): can't open the file or file doesn't exist");
		exit(1);
	}
	
	unsigned long long time_us = 0;
	int gx = 0;
	int gy = 0;
	int gz = 0;
	int ax = 0;
	int ay = 0;
	int az = 0;

	unsigned long count = 0;

	/* ignore the first line */
	char ch = ' ';
	long int seek_pos;
	int pp = 0;

	while (1) {
		ch = fgetc(fp);
		fputc(ch, stdout);
		if (ch == '\n') {
			break;
		}
		
		pp++;
	}
	seek_pos = ftell(fp);
	
	/* count lines */
	count = 0;
	while(!feof(fp)) {
		if (fscanf(fp, "%llu %d %d %d %d %d %d",
			&time_us, &gx, &gy, &gz, &ax, &ay, &az) != 7) printf("Something wrong in columns was detected\n"); 
		count++;
	}

	printf("Number of lines in %s is %lu\n", file_name, count);	

	if (data == NULL) {
		printf("No matrix pointer was given in arguments.\n");
		fclose(fp);
		return(-1);
	}

	(*data) = create_fmat(count, 7);

	count = 0;
	fseek(fp, seek_pos, SEEK_SET); // set begining of file

	while (!feof(fp)) {
		fscanf(fp, "%llu %d %d %d %d %d %d",
				&time_us, &gx, &gy, &gz, &ax, &ay, &az); 
		(*data)->data[count][0] = time_us;

		(*data)->data[count][1] = gx;
		(*data)->data[count][2] = gy;
		(*data)->data[count][3] = gz;

		(*data)->data[count][4] = ax;
		(*data)->data[count][5] = ay;
		(*data)->data[count][6] = az;

		count++;
		
	}
	printf("Data rows were read: %lu\n", count);

	fclose(fp);
	return(1);
}

double variance(Fmatrix_ptr data, unsigned long int column_id, unsigned long int index_from, unsigned long int index_to, double shift) {
	/* returns variance of a dataset defined by index boundary */

	if (((data->rows - 1) < index_to) || (index_from < 0) || (data->rows == 0)) {
		fprintf(stderr, "Error in variance(): the given indices are not correct: %lu.\n", index_from);
		exit(1);
	}

	if (index_from >= index_to) {
		fprintf(stderr, "Error in variance(): index value from should be greater than index value to from: %lu to %lu.\n", index_from, index_to);
		exit(1);
	}

	if ((column_id < 0) || (data->cols == 0) || (column_id > data->cols - 1)) {
		printf("%lu %lu\n", column_id, data->cols);
		fprintf(stderr, "Error in variance(): column index is out of range: from %lu to %lu.\n", index_from, index_to);
		exit(1);
	}

	/* estimate mju value */
	double u = 0.0;
	for (long int i = index_from; i <= index_to; i++) {
		u += (double)data->data[i][column_id] + shift;
	}
	u /= (double)(index_to - index_from + 1);

	/* estimate variance */
	double var = 0.0;
	for (long int i = index_from; i <= index_to; i++) {
		var+= pow((double)data->data[i][column_id] - u, 2.0);
	}
	var /= (double)(index_to - index_from + 1);
	
	return(var);

}

unsigned short detect_intervals(Fmatrix_ptr raw_data_in, Fmatrix_ptr *inter_out_data, 
							Fmatrix_ptr *csv_out_data, double w_t, double min_slen, double min_dlen) {
	/* the method reads IMU raw data were loaded from a csv file's into the raw_data_in matrix 

	Input:
		1) raw_data_in - matrix has 7 columns of data: us_time, ax, ay, az, gz, gy, gz;
		2) w_t - parameter is a time window in seconds which is used for magnitude estimation;
		3) min_slen - the minimal length of static intervals;
		4) min_dlen - the minimal length of dynamic intervals.

	Output:
		1) inter_out_data - the matrix with 5 columns: index_from, index_to, static|dynamic, [log(var(magnitude)) - estimated only for static intervals], reserved;
		2) csv_out_data - the output raw matrix is additionally marked with static and dynamic labels ( 1 - static, -1 - dynamic)];
		3) number of detected static intervals.

	*/

	/* estimate number of records within the time window w_t */

	#ifndef MAX_VAL
		#define MAX_VAL 32676.0
	#endif

	unsigned long int w_init = raw_data_in->data[0][0];
	short w_len = 0;
	short w_shift = 0; // used in shrinking a static interval

	for (unsigned long int i = 1; i < raw_data_in->rows; i++) {
		if ((raw_data_in->data[i][0] - w_init) / 1.0e+6 >= w_t) {
			w_len = i + 1;
			w_shift = w_len + (int)w_len / 2.0;

			break;
		}
	}

	short min_static_len = (short)(float)w_len / w_t * min_slen; // minimal number of records in a static interval
	short min_dynamic_len = (short)(float)w_len / w_t * min_dlen; // minimal number of records in a dynamic interval

	unsigned short wt_h = (unsigned short)w_t / 2.0 + 1; // half of w_t
	printf("Time window of %.2f is equal to %d records, and shift is %d, min static len is %d\n", w_t, w_len, w_shift, min_static_len);

	Fmatrix_ptr mag_values = create_fmat(raw_data_in->rows - w_t, 2); // one column matrix for estimation of variance's magnitude for ax, ay, az

	double vax, vay, vaz;
	double mag_avg = 0.0;
	long int counter = 0;
	double th = log10(1.0);

	/* estimate variance's magnitude and make logarithmic scaling */
	for (unsigned long int i = wt_h; i < raw_data_in->rows - wt_h; i++) {
		vax = variance(raw_data_in, 4, i - wt_h, i + wt_h, MAX_VAL); // var(ax)
		vay = variance(raw_data_in, 5, i - wt_h, i + wt_h, MAX_VAL); // var(ay)
		vaz = variance(raw_data_in, 6, i - wt_h, i + wt_h, MAX_VAL); // var(az)

		mag_values->data[i - wt_h][0] = sqrt(pow(vax, 2.0) + pow(vay, 2.0) + pow(vaz, 2.0)); 
		mag_values->data[i - wt_h][1] = log10(mag_values->data[i - wt_h][0] + 1.0);  // adding a small value for the log function

		/* check magnitude's value is larger than log(1.0) - dynamic state  */
		if (mag_values->data[i - wt_h][1] > th) {
			mag_avg += mag_values->data[i - wt_h][1]; // in below used for average estimation
			counter += 1;
		}
	}

	/* Average of mag's values */
	mag_avg /= counter;
	printf("Estimated mag threshold is: %.3f\n", mag_avg);

	/* detect static and dynamic intervals */
	long int start = -1;
	long int dyn_start = 0;
	short stat_num = 0;
	short any_interv = 0;
	long int intervals[1000][3]; // the maximal number of intervals in will detected. The last column is type of an interval: 1 static, -1 dynamic
	short detected_intervals = 0; 

	for (unsigned long int i = 0; i < mag_values->rows; i++) {
		
		if (mag_values->data[i][1] <= mag_avg && start < 0) { // we are at beginning of a static interval
			start = i; // set start index of the static interval
		} else if (start >= 0 && (mag_values->data[i][1] > mag_avg || i == mag_values->rows - 1)) { // we are at the end of the static interval
			if (i - start >= min_static_len) { // check a lenght of the static interval
				if (any_interv != 0) // check if any interval was detected before
					 detected_intervals += 1; // count the static interval
				intervals[detected_intervals][0] = start; // the start index of the static interval
				intervals[detected_intervals][1] = i - 1; // the end index of the static interval
				intervals[detected_intervals][2] = 1; // mark as static

				stat_num += 1; // count static intervals
				dyn_start = i; // set beginning of a dynamic interval
				start = -1; // reset beginning of the static interval
				detected_intervals += 1; // count intervals
				any_interv = 1; // mark that a static or a dynamic interval was detected
			} else {

				start = -1; // we are in a dynamic interval
 			}
		}

		if (start == -1 && dyn_start >= 0) { // update the dynamic interval's the end index
			intervals[detected_intervals][0] = dyn_start ; 
			intervals[detected_intervals][1] = i; 
			intervals[detected_intervals][2] = -1; // mark the interval as dynamic
			any_interv = 1;
		}
	}

	/* generate the output data */

	*inter_out_data = create_fmat(detected_intervals, 5); // id, from_index, to_index, type: {-1 dyn | 1 stat | -2 small dyn}, | quality (only for static), reserved
	Fmatrix_ptr out_data = create_fmat(raw_data_in->rows, 8); // us_time, ax, ay, az, gz, gy, gz, static|dynamic
	*csv_out_data = out_data;
	
	unsigned long row = 0;
	for (unsigned short j = 0; j < detected_intervals; j++) {

		unsigned long int start_id = intervals[j][0];
		unsigned long int stop_id = intervals[j][1];
		short interv_type = intervals[j][2];
	
	
		/* align a static interval (shift both borders of the interval) */
		if (interv_type == 1) {
			if ((start_id + w_shift) >= start_id) start_id += w_shift; // check overflow
			if (w_shift <= stop_id) stop_id -= w_shift; // check overflow
			}
		else { // dynamic interval
			if (w_shift <= start_id) start_id -= w_shift; // check overflow
			if ((stop_id + w_shift) >= stop_id) stop_id += w_shift; // check overflow
			}

		/* check final indices */
		if (start_id < 0) start_id = 0;

		if (stop_id > raw_data_in->rows - 1) stop_id = raw_data_in->rows - 1;

		if (start_id > stop_id) start_id = stop_id - 1;
	
		/* for the static interval estimate variance of magnitude */
		if (interv_type == 1) {
			(*inter_out_data)->data[j][3] = variance(mag_values, 0, start_id, stop_id, 0.0);
		}
	
		/* assign aligned indices */
		(*inter_out_data)->data[j][0] = (double)start_id;
		(*inter_out_data)->data[j][1] = (double)stop_id;
		(*inter_out_data)->data[j][2] = (double)interv_type;

		/* check length of the dynamic interval */
		double interv_t = (double)interv_type;
		if ((interv_type == -1) && ((stop_id - start_id + 1) < min_dynamic_len)) {
				interv_t = -2.0; // mark as too small
				(*inter_out_data)->data[j][2] = interv_t;
			}

		/* fill the output matrix */
		for (long int i = start_id; i <= stop_id; i++) {
			out_data->data[i][0] = raw_data_in->data[i][0]; // us_time
			out_data->data[i][1] = raw_data_in->data[i][1]; // ax
			out_data->data[i][2] = raw_data_in->data[i][2]; // ay
			out_data->data[i][3] = raw_data_in->data[i][3]; // az
			out_data->data[i][4] = raw_data_in->data[i][4]; // gx
			out_data->data[i][5] = raw_data_in->data[i][5]; // gy
			out_data->data[i][6] = raw_data_in->data[i][6]; // gz
			out_data->data[i][7] = interv_t; // static | dynamic | too_small_dynamic
			row++;
		}
	}

	/* normalize variance of magnitude for static intervals  */
	double s = 0.0;
	for (unsigned short j = 0; j < detected_intervals; j++) {
		if ((*inter_out_data)->data[j][2] == 1.0) s += (*inter_out_data)->data[j][3];
	}

	for (unsigned short j = 0; j < detected_intervals; j++) {
		if ((*inter_out_data)->data[j][2] == 1.0) (*inter_out_data)->data[j][3] /= s;
	}


	free_matrix(mag_values);
	return(stat_num);	
	
}


Fmatrix_ptr create_fmat(unsigned long int rows, unsigned long int cols){
	/* creates matrix of double elements: rows and cols */
	
	/* allocate memory for array data */
	double **new_data;
	new_data = (double **)malloc(sizeof(double *) * rows);

	for (unsigned long int i = 0; i < rows; i++){
		/* create a list of elements for each row */
		new_data[i] = (double *)malloc(sizeof(double) * cols);
	}

	/* init the array with zero values */

	for (unsigned long int i = 0; i < rows; i++){
		for (unsigned long int j = 0; j < cols; j++) {
			new_data[i][j] = 0.0;
		}
	}

	/* allocate memory for struct data */
	Fmatrix_ptr new_matrix = (Fmatrix_ptr)malloc(sizeof(struct double_matrix)); // allocate memory for pointer of struct double_matric

	new_matrix->rows = rows;
	new_matrix->cols = cols;
	new_matrix->data = new_data;

	return(new_matrix);
}


void free_matrix(Fmatrix_ptr matrix){
	/* free memory for matrix */

	for (unsigned long int i = 0; i < matrix->rows; i++){
		free(matrix->data[i]);
	}

	free(matrix->data);
	free(matrix);

}


void print_fmat(Fmatrix_ptr matrix){
	/* prints matrix on the screen */

	for (unsigned long int i = 0; i < matrix->rows; i++) {
		for(unsigned long int j = 0; j < matrix->cols; j++){
			printf("%8.4f ", matrix->data[i][j]);
		}
		printf("\n");
	}
	printf("%lu rows and %lu cols\n\n", matrix->rows, matrix->cols);
}


void print_mat_shape(Fmatrix_ptr mat) {
	/* prints shape of the matrix */
	printf("The mat shape is %lu rows and %lu cols\n", mat->rows, mat->cols); 
}

