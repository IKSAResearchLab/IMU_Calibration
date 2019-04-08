Use this source code at you own risk.
You will be able to modify and to reuse any part of this source code.

This source code demonstrates the improved version of the static detector used within the IMU calibration procedure.
This code was tested on Ubuntu 16.04.

Usage:

1) In the main.c check the offset value for MAX_VAL (depends on your IMU device's data range for minimal negative value).
2) Type 'make' in the command line to compile the code.
3) In the command line type './detector ./stat_sample0.csv 0.2 1.0 0.4'. The first argument is a file of IMU raw with 7 columns:
	ustime - time in microseconds;
	gx, gy, gz - gyroscope's data;
	ax, ay, az - accelerometr's data;
	
	From the second till the forth argument you may set such parameters: 0.2 - window size in seconds for variance estimations; 1.0 - the minimal length of static interval in seconds, 0.4 the minimal length of the dynamic interval in seconds.
