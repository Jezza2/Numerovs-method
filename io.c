/*

NUMERICAL SOLUTIONS TO SCHRODINGER'S EQUATIONS

Author:  Jeremy Stanger
Date:    13/09/16

This file is basically the housekeeping; no interesting computing
going on here.
*/

#include <stdio.h>
#include <stdlib.h>
#include "header.h"


/* Write a gnuplot script to produce the required plot
	 Plot depends on programme execution mode.
	 Gnuplot can then be called externally

	 order: The order (n) of the solution to be plotted
	 plot_alyt: flag - 1 to plot analytical solution
	 plot_num : flag - 1 to plot numerical solution
	 */
void plot (int order, int plot_alyt, int plot_num) {

	// Check something is actually being plotted!
	if (plot_alyt || plot_num) {
		// File pointer
		FILE *fp;
		
		// Character buffer for file name
		char filename[18] = { };
		// Set file name according to solution order
		sprintf(filename, "plots/plot_%d.p", order);
		// Open file for writing if possible, otherwise
		// quit with error message
		if ((fp = fopen(filename, "w")) == NULL) {
			printf("Unable to open %s for plotting\n", filename);
			_exit(4);
		}

		// character buffer for the line doing the plotting
		char plot_line[150] = { };
		
		// Plot only analytical solution
		if (plot_alyt && !plot_num) {
			sprintf(plot_line, 
				"set key off\nplot \"../data/order_%d.dat\" using 1:2 with lines\n", 
				order);
		// Plot only numerical solution
		} else if (!plot_alyt && plot_num) {
			sprintf(plot_line,
				"set key off\nplot \"../data/order_%d.dat\" using 1:3 with lines\n",
				order);
		// Plot both numerical and analytical solutions
		} else {
			sprintf(plot_line,
				"set key on\nplot \"../data/order_%d.dat\" using 1:2 with lines title \"analytical\"," 
				" \"../data/order_%d.dat\" using 1:3 with lines title \"numerical\"\n",
				 order, order);
		}

		// Write script to file
		fprintf(fp,
			"set terminal jpeg size 1000,750\n"
			"set output \"plot_%d.jpg\"\n"
			
			"set title \"%dth energy level quantum oscillator wavefunction\"\n"
			"set xrange [0:%.9g]\n"
			"set autoscale y\n"
			"set zeroaxis\n"
			"set xlabel \"x\"\n"
			"set ylabel \"PHI_%d(x)\"\n"
			"set xtic auto\n"
			"set ytic auto\n"
			"set grid\n"
	
			"%s",
			order, order, interval_width, order, plot_line
		);

		// Close file
		fclose(fp);
	}
}

/* Writes data in arrays to a file for plotting purposes
	 Space delimited file:
	 t y1 y2 y3

	 index: order (n) of the solution to be plotted
	 */
void write_datafile(int index) {
	// File pointer
	FILE *fp;
	//  character buffer for file name
	char filename[20] = { };
	// index variable
	int i;
	
	// Set file name according to the solution order
	sprintf(filename, "data/order_%d.dat", index);
	// Open file for writing if possible, otherwise
	// exit with error message
	if ((fp = fopen(filename, "w")) == NULL) {
		printf("Unable to open %s to write data\n", filename);
		_exit(5);
	};
	fprintf(fp, "#x y_alyt y_num\n");

	// Write data using loop
	if (mode == 0) {
		for (i = 0; i < num_points; i++) {
			fprintf(fp, "%.16g %.16g %.16g\n", 
				(alyt_sol + i)->x,
				(alyt_sol + i)->y,
				0.0);	
		}
	} else if (mode == 1) {
		for (i = 0; i < num_points; i++) {
			fprintf(fp, "%.16g %.16g %.16g\n", 
				(num_sol + i)->x,
				0.0,
			  (num_sol + i)->y);
		}
	} else if (mode == 2 || mode == 4) {
		for (i = 0; i < num_points; i++) {
			fprintf(fp, "%.16g %.16g %.16g\n", 
				(alyt_sol + i)->x,
				(alyt_sol + i)->y,
				scale_factor()*(num_sol + i)->y);
		}
	}

	// Close file
	fclose(fp);
}

/* Get user parameters and set global variables accordingly.
	 Also allocate memory required for the data storage.

	 count:  Number of parameters
	 *argvec[]: pointer to string array of the parameters
	 */
void set_params(int count, char *argvec[]) {
	// If _exit() is called before we've acllocated memory then
	// we'll be trying to free memory that was never allocated.
	// free() does nothing with NULL pointers though.
	num_sol == NULL;
	alyt_sol == NULL;

	// Check we've got no inputs then immediately show help
	if (count == 1) {
		help();
		_exit(2);
	}

	// Get execution mode number
	mode = atoi(*(++argvec));

	// Check we've got the right number of parameters for
	// the mode
	switch (mode) {
		case 0: case 1: case 3:
			if (count != 5) {
				help();
				_exit(2);
			}
			break;
		case 2: case 4:
			if (count != 6) {
				help();
				_exit(2);
			}
			break;
	}

	// Set interval_width and delta so that we can calculate
	// number of points
	if ((interval_width = (double)atof(*(++argvec))) && (delta = (double)atof(*(++argvec)))) {

		num_points = N_POINTS(interval_width, delta);

		// Allocate memory for the data points.  Quit with error code 1 if request for memory is refused.
		if (mode == 1 || mode == 2 || mode == 3 || mode == 4) {
			if ((num_sol = malloc(sizeof(struct data_point) * num_points)) == NULL)
				_exit(1);
		}

		if (mode == 0 || mode == 2 || mode == 4) {
			if ((alyt_sol = malloc(sizeof(struct data_point) * num_points)) == NULL)
				_exit(1);
		}
	}	else {
		help();
		_exit(2);
	}

	// Allocate global variables according to execution mode
	switch (mode) {
		default:
			help();
			_exit(2);
			break;

		case 0:
			if ((n = atoi(*(++argvec))) < 0) {
				help();
				_exit(2);
			}
			break;

		case 2:
		case 4:
			if ((n = atoi(*(++argvec))) < 0) {
				help();
				_exit(2);
			}
			
		case 1:
		case 3:
			if (!(E = (double)atof(*(++argvec)))) {
				help();
				_exit(2);
			}
			break;
	}
}

/* Exit cleanly, freeing any dynamically allocated memory first.
	 Good practice.  Display exit status.

	 code: Exit status
	 */
void _exit(int code) {
	free(num_sol);
	free(alyt_sol);
	printf("Exit status: %d\n", code);
	exit(code);
}

/* Display parameter list and modes
	 */
void help(void) {
	printf("\nInputs required as follows in the following orders.\n"
         "There are 5 different modes available.\n\n\n");

	printf("MODE 0: PLOT ANALYTICAL SOLUTION\n\n"
				 "int mode     Mode number\n"
				 "double x1    The interval width\n"
				 "double d     The data point separation\n"
				 "int n        Quantum energy level\n\n\n");

	printf("MODE 1: PLOT NUMERICAL SOLUTION\n\n"
				 "int mode     Mode number\n"
				 "double x1    The interval width\n"
				 "double d     The data point separation\n"
				 "double E     The energy of the particle\n\n\n");

	printf("MODE 2: PLOT BOTH SOLUTIONS\n\n"
				 "int mode     Mode number\n"
				 "double x1    The interval width\n"
				 "double d     The data point separation\n"
				 "int n        Quantum energy level\n"
				 "double E     The energy of the particle\n\n\n");

	printf("MODE 3: CALCULATE E\n\n"
				 "int mode     Mode number\n"
				 "double x1    The interval width\n"
				 "double d     The data point separation\n"
				 "double E     The guess at the energy\n\n\n");
	
	printf("MODE 4: EULER'S METHOD\n\n"
				 "int mode     Mode numer\n"
				 "double x1    The interval width\n"
				 "double d     The data point separation\n"
				 "int n        Quantum energy level\n"
				 "double E     The energy of the particle\n\n\n\n");

	printf("EXIT STATUSES:\n\n"
				 "0 -          Successful execution\n"
				 "1 -          Unable to allocate memory for data arrays\n"
				 "2 -          Input error\n"
				 "3 -          Insufficient domain for calculation of E\n"
				 "4 -          Unable to open file for plotting\n"
				 "5 -          Unable to open file to write data\n\n");
}
