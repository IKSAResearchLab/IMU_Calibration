/* 

The program demonstrates a new version of the static detector used in the IMU calibration procedure, see 

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

As an input is used CSV file (space idented) with such columns: ustime, gx, gy, gz, az, ay, az
The output is CSV file with columns: ustime, gx, gy, gz, ax, ay, az, interval_type(0-static, 1-dynamic, 2- too short dynamic)

*/

#include "stat_detector.h"

#define MAX_VAL 32676.0

int main(int argc, char **argv){

	double time_window = 0.5; // in seconds
	double min_static_len = 1.0; // in seconds
	double min_dynamic_len = 0.4; // in seconds

	Fmatrix_ptr data;

	if (argc >= 2) {
		printf("CSV file is given: %s\n", argv[1]);
		if (!read_csv(argv[1], &data)) {
			fprintf(stderr, "Can't read the data from the CSV file.");
			exit(1);
		}
	} else {
		printf("Please, provide a CSV file name in arguments\n");
		free_matrix(data);
		return(0);
	}

	if (argc == 5) {
		time_window = (double)atof(argv[2]);
		min_static_len = (double)atof(argv[3]);
		min_dynamic_len = (double)atof(argv[4]);
	}
	printf("Detector's parameters are used: win lenght = %.2f sec, min static lenght = %.2f sec and min dynamic length = %.2f sec\n", time_window, min_static_len, min_dynamic_len);

	Fmatrix_ptr detected_intervals;
	Fmatrix_ptr output_data;

	clock_t t = clock();
	unsigned int int_num = detect_intervals(data, &detected_intervals, &output_data,
											time_window, min_static_len, min_dynamic_len);
	t = clock() - t;
	
	if (int_num > 0) {
		printf("Intervals detected in %.2f sec\n", ((double)t)/CLOCKS_PER_SEC);
		printf("Detection results:\n");
		printf("%8s %8s %8s %8s %8s\n", "Start", "Stop", "Type", "Quality", "Reserved");
		print_fmat(detected_intervals);

		save_matrix2csv(detected_intervals, "intervals.csv", NULL);
		printf("Detected intervals were saved in: intervals.csv\n");
		const char *header [] = {"ustime", "gx", "gy", "gz", "ax", "ay", "az", "static"};
		save_matrix2csv(output_data, "out.csv", header);
		printf("Output data were saved in: out.csv\n");
	} else {
		printf("Static intervals are not detected.\n");
	}

	free_matrix(data);
	free_matrix(detected_intervals);
	free_matrix(output_data);

	return(1);

}
