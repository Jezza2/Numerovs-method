/*

NUMERICAL SOLUTIONS TO SCHRODINGER'S EQUATION

Author:  Jeremy Stanger
Date:    13/09/2016

Main file; where the business happens.

The program will produce analytical solutions to the time independent 
schrodinger equation for a particle in a harmonic potential.  
Alongside it will approximate a numerical solution using Numerov's method.
It will approximate a numerical solution using Euler's method
It will also produce a correct value of E given a value that is close to it.

The output is determined by the mode you are in which is an input parameter.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "header.h"

// Define some global variables
// Input parameters
// The range of x values to use
double interval_width = 0;
// The difference between 2 successive values of x
double delta = 0;
// The energy
double E = 0;
// The order of the solution
int n = 0;
// The mode of programme execution
int mode = 0;
// Derived from parameters
int num_points = 0;

int main(int argc, char *argv[]) {
	set_params(argc, argv);

	// Calculate a value of n derived from E.
	// Only applicable in modes where n is not
	// available
	int n_approx = (int)round((E-1.0)/2.0);
	switch (mode) {
		default:
			_exit(0);
			break;

		case 0: // Plot analytical solution
			calc_alyt(n);
			write_datafile(n);
			plot(n, 1, 0);
			break;

		case 1: // Plot numerical (numerov) solution
			calc_num(n_approx);
			write_datafile(n_approx);
			plot(n_approx, 0, 1);
			break;

		case 2: // Plot analytical and numerov solution
			calc_alyt(n);
			calc_num(n);
			write_datafile(n);
			plot(n, 1, 1);
			break;

		case 3: // Calculate E
			calc_E(n_approx);
			break;

		case 4: // Plot analytical and Euler solution
			calc_alyt(n);
			calc_num(n);
			write_datafile(n);
			plot(n, 1, 1);
			break;
	}

	printf("Successful execution\n");

	_exit(0);
}

/* Populate alyt_sol with the analytical solution to 
	 Schrodinger's equation of order index

	 index: order (n) of the solution
	 */
void calc_alyt(int index) {
	// index variable
	int i;
	// Populate array with solution calculated at each
	// value of x, analytically
	for (i = 0; i < num_points; i++) {
		(alyt_sol + i)->x = ((double)i)*delta;
		(alyt_sol + i)->y = exp(-pow((alyt_sol + i)->x, 2.0) / 2.0) * hermite(index, (alyt_sol + i)->x);
	}
}

/* Populate num_sol with the numerical solution to
	 Schrodinger's equation of order index

	 index: order (n) of the solution
	 */
void calc_num(int index) {
	// index variable
	int i;
	// Initial gradient of solutino
	double dphi0;

	num_sol->x = 0.0;

	// Set initial conditions according to the 
	// parity of n (= parity of wavefunction)
	if (PARITY(index) == EVEN) {
		num_sol->y = 1.0;
		dphi0 = 0.0;
	} else {
		num_sol->y = 0.0;
		dphi0 = 1.0;
	}

	// Calculate solution using method according
	// to execution mode.  Calculated at each value
	// of x.
	switch (mode) {

		// These all use Numerov's method
		case 1: case 2: case 3:
			(num_sol + 1)->x = delta;
			(num_sol + 1)->y = second_point(num_sol->y, dphi0, index);

			// We start at i=2 since 2 points are pre-
			// calculated for use of Numerov's method.
			for (i=2; i < num_points; i++) {
				(num_sol + i)->x = ((double)i)*delta;
				(num_sol + i)->y = numerov(i); 
			}
			break;

		// Euler's method
		case 4:
			for (i=1; i < num_points; i++) {
				(num_sol + i)->x = ((double)i)*delta;
				(num_sol + i)->y = euler(i, dphi0);
			}
			break;
	}
}



/* Return the value of H_n(x), the nth order Hermite
	 polynomial

	 index: n - the order of the polynomial
	 x: The value of x at which to calculate the polynomial
   */
double hermite(int index, double x) {
	// Recursion loop terminates either at index=0
	// and index=1.
	if (index == 0) {
		return H_0;
	} else if (index == 1) {
		return H_1(x);
	} else {
		// Recursive algorithm
		return 2.0*x*hermite(index-1, x) - 2.0*((double)index-1.0)*hermite(index-2, x);
	}
}

/* Calculate and return next point calculated by Numerov's 
	 method.

	 index: order (n) of the solution
	 */
double numerov(int index) {
	// Straightforward implemntation of Numerov's equation
	return ( (2.0 + 5.0 * delta*delta * f( (num_sol+index-1)->x ) / 6.0 ) * ( (num_sol+index-1)->y )
					- (1.0 - delta*delta * f( (num_sol+index-2)->x ) / 12.0 ) * ( (num_sol+index-2)->y ) )
					/ (1.0 - delta*delta * f( (num_sol+index  )->x ) / 12.0 );
}

/* Calculate and return next point calculated by Euler's 
	 method.

	 index: order (n) of the solution
	 dphi0: initial gradient
	 */
double euler(int index, double dphi0) {
	// Second order ODE requires integrating twice.  This is a running total
	// of the first integral.
	static double dphi;
	if (index == 1) dphi = dphi0;
	// Integrate once
	else if (index > 1) dphi += delta*f((num_sol+index-1)->x)*(num_sol+index-1)->y;
	
	// Integrate twice
	return (num_sol+index-1)->y + dphi*delta;
}

/* It is an observational fact that when E is set too high, the calculated curve
	 diverges to infinity in one direction, and to infinity in the other direction	
 	 when E is set too low.  We can use a simple limit to determine where the curve
 	 has diverged.  While the direction of divergence is not predictable, we can
 	 determine whether the given value of E is too large or too small by making an
 	 arbitrary change to E and retesting the value of the curve at the same point.
 	 If the divergence is greater after the correction, then we know the true value
 	 of E lies somewhere in 'the other direction.'

	 This function uses this method to recursively refine the value of E from an
	 initial trial value to the correct value.  The accuracy is about 0.00001.
	 The result is stored in the global variable E.

	 index: the order (n) of the solution being used.
	 */
void calc_E(int index) {
	// index variable
	int i;
	// initial correction to E
	double interval_size = E / 100.0;
	// Target precision to which E will be calculated
	double precision = 0.00001;
	// E during previous iteratoin
	double E_prev = E;

	// Flags most efficiently contained in a bit-field
	struct {
		unsigned int halve : 1; // halve the interval size every iteration
		unsigned int change_direction : 1;  // change direction of applied correction when set.
		unsigned int done_once : 1;	// need to run process once before applying corrections.
	} flags;

	// set all flags to zero
	flags.halve = flags.change_direction = flags.done_once = 0;

	struct data_point prev_div; // previous point at which solution diverged

	do {
		i = 0;
		calc_num(index);
		flags.change_direction = 0;

		// cycle through data points until diverged.
		// we know first point hasn't diverged since that's the initial condition!
		while (fabs((num_sol+(++i))->y) < 100.0) {
			if (i > num_points) {
				printf("Insufficient domain for calculation!\n");
				_exit(3);
			}
		}

 /* Here we determine what to do next:
		
		If we have modified E but the direction of divergence hasn't flipped
		and the function diverges sooner then we should flip the direction
		without decreasing interval size.
	
		If we have modified E but the direction of divergence hasn't flipped
		but the function diverges later, then we should do nothing

		If we've modified E and the direction of divergence has flipped, then
		we should flip the direction and halve the interval size.

		If the direction of divergence has flipped once then we have found 
		the interval in which the solution lies, so we can always halve the 
		interval size at least once. */

		if (flags.done_once) {
			// If the direction of divergence hasn't flipped.
			if (SIGN((num_sol+i)->y) == SIGN(prev_div.y)) {
				// If the new value of the calculated function at the previous x-value
				// of divergence is greater than the old value at divergence.
				// Amounts to testing if the divergence has increased or decreased
				// in magnitude from the last value of E.  If it has, change the
				// direction of the correction to E.
				if (fabs((num_sol+(int)(prev_div.x/delta))->y) > fabs(prev_div.y) 
						&& !flags.halve) {
			 		
					flags.change_direction = 1;
					
				// Halve the interval size no more than once per flip.
				} else if (flags.halve) {
					interval_size /= 2;
					flags.halve = 0;
				}
			
			// Direction of divergence has flipped - halve the interval size and
			// reverse the direction of correction to E.
			} else {
				flags.halve = 1;
				interval_size /= 2;
				flags.change_direction = 1;
			}
		}
		// Record the last point at which the solution diverged.
		prev_div = *(num_sol+i);

		// Set the direction in which to correct E.
		interval_size *= flags.change_direction ? -1 : 1;
		// Record the last value of E.
		E_prev = E;
		// Modify E
		E += interval_size;
		flags.done_once = 1;

		printf("%.16g\t\t\t%.16g\n", interval_size, E);
		// The second condition is in an effort to  ensure that the last
		// two values are either side of the true value so that a meaningful
		// mean can be calculated.
	} while (fabs(interval_size) > precision || !flags.change_direction);

	// calculate the mean of the last two values.
	E = (E + E_prev)/2;

	printf("\n\nE = %.5g\n\n", E);
}

/* Calculate the second point of the solution using the initial
	 conditions and the Maclaurin expansion.

	 phi0: initial value
	 dphi0: initial gradient
	 order: order (n) of the solution
	 */
double second_point(double phi0, double dphi0, int order) {
  // Taylor expansion depends on parity of wavefunction
	if (PARITY(order) == EVEN) {
		return phi0 + (pow(delta, 2.0) / 2.0)*f(0.0)*phi0 + (pow(delta, 4.0) / 24.0)*(ddf*phi0 + 2.0*df(0.0)*dphi0 + pow(f(0.0), 2.0)*phi0);
	} else {
		return delta*dphi0 + (pow(delta, 3.0) / 6.0)*(f(0.0)*dphi0 + df(0.0)*phi0);
	}
}

/* A solution multiplied by a constant factor is still a solution.
   The numerical and analytical solutions may not have the same magnitude.
   This returns the scale factor to correct for this discrepancy.
	 It does this by the ratio of the tenth points in each solution.
	 Only relevant when both the numerical and analytical solutinons
	 are plotted together.
	 */
double scale_factor(void) {
	return (alyt_sol+10)->y/(num_sol+10)->y;
}

/* Calculates derivative of the numerical solution using a 
	 3 point finite difference method at the point corresponding
	 to the centre_index element of the array.

	 centre_index: The array index of the point at which to find
								 the derivative.
	 */
double stencil_3_point(int centre_index) {
	return ((num_sol+centre_index-1)->y - 2*(num_sol+centre_index)->y + (num_sol+centre_index+1)->y)
				/ (delta*delta); 
}
