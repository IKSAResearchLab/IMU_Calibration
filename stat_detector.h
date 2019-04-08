#ifndef LSM6DSL_DRIVER_H_INCL
#define LSM6DSL_DRIVER_H_INCL


/* A static detector is an improved version of the static detector firstly introduced in: 

@article{Tedaldi2014,
archivePrefix = {arXiv},
arxivId = {ISSN 1476-2986},
author = {Tedaldi, David and Pretto, Alberto and Menegatti, Emanuele},
doi = {10.1109/ICRA.2014.6907297},
eprint = {ISSN 1476-2986},
journal = {1. Tedaldi D, Pretto A, Menegatti E. A robust easy to implement method IMU calibration without Extern. equipments. Proc - IEEE Int Conf Robot Autom. 2014;3042–9. Proc. - IEEE Int. Conf. Robot. Autom.},
pages = {3042--3049},
pmid = {16830938},
title = {{A robust and easy to implement method for IMU calibration without external equipments}},
year = {2014}
}

The new method created by Dr. Andrejs Zujevs, Riga Technical University, Faculty of Computer Science and Information Technology.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
	
/* Defining a new data type */
struct double_matrix {
	double **data;
	long unsigned int rows;
	long unsigned int cols;
};

typedef struct double_matrix* Fmatrix_ptr;

Fmatrix_ptr create_fmat(unsigned long int rows, unsigned long int cols);
void print_fmat(Fmatrix_ptr matrix);
void free_matrix(Fmatrix_ptr matrix);
void print_mat_shape(Fmatrix_ptr mat);
int read_csv(char *file_name, Fmatrix_ptr *data);
double variance(Fmatrix_ptr data, unsigned long int column_id, unsigned long int index_from, unsigned long int index_to, double shift); 
int save_matrix2csv(Fmatrix_ptr data, char *filename, const char **header);
unsigned short detect_intervals(Fmatrix_ptr raw_data_in, Fmatrix_ptr *inter_out_data, 
							Fmatrix_ptr *csv_out_data, double w_t, double min_slen, double min_dlen); 

#endif
