/*

NUMERICAL SOLUTIONS TO SCHRODINGER'S EQUATIONS

Author:  Jeremy Stanger
Date:    13/09/2016

Function prototypes, macros and external variable declarations

*/

// Some definitions and macros
// Number of data points given interval width and delta
#define N_POINTS(a, b) (int)((a)/(b)+2)

// Definitions for bottom of recursion loop for Hermite polynomials 
#define H_0 1.0
#define H_1(x) 2.0*(x)

// Macros to determine parity of a number
#define PARITY(n) (n)%2
#define EVEN 0
#define ODD 1

// Potential function and its derivatives
#define V(p) ((p)*(p))
#define dV(p) (2.0*(p))
#define ddV 2.0

// f(x) as used in Numerov's method and its derivatives
#define f(p) (V(p)-E) 
#define df(p) dV(p)
#define ddf ddV

// The sign of n
#define SIGN(n) (((n) > 0) - ((n) < 0))

// Function prototypes

// Functions in io.c
/* Display parameter list and modes
	 */
void help(void);   

/* Exit cleanly, freeing any dynamically allocated memory first.
	 Good practice.  Display exit status.

	 code: Exit status
	 */
void _exit(int code);

/* Get user parameters and set global variables accordingly.
	 Also allocate memory required for the data storage.

	 count:  Number of parameters
	 *argvec[]: pointer to string array of the parameters
	 */
void set_params(int count, char *argvec[]);

/* Writes data in arrays to a file for plotting purposes
	 Space delimited file:
	 t y1 y2 y3

	 index: order (n) of the solution to be plotted
	 */
void write_datafile(int index);

/* Write a gnuplot script to produce the required plot
	 Plot depends on programme execution mode.
	 Gnuplot can then be called externally

	 order: The order (n) of the solution to be plotted
	 plot_alyt: flag - 1 to plot analytical solution
	 plot_num : flag - 1 to plot numerical solution
	 */
void plot(int order, int plot_alyt, int plot_num);


// Functions in schrodinger.c
/* Return the value of H_n(x), the nth order Hermite
	 polynomial

	 index: n - the order of the polynomial
	 x: The value of x at which to calculate the polynomial
   */
double hermite(int index, double x);

/* Populate alyt_sol with the analytical solution to 
	 Schrodinger's equation of order index

	 index: order (n) of the solution
	 */
void calc_alyt(int index);

/* Populate num_sol with the numerical solution to
	 Schrodinger's equation of order index

	 index: order (n) of the solution
	 */
void calc_num(int index);

/* Calculate the second point of the solution using the initial
	 conditions and the Maclaurin expansion.

	 phi0: initial value
	 dphi0: initial gradient
	 order: order (n) of the solution
	 */
double second_point(double phi0, double dphi0, int order);

/* Calculate and return next point calculated by Numerov's 
	 method.

	 index: order (n) of the solution
	 */
double numerov(int index);

/* A solution multiplied by a constant factor is still a solution.
   The numerical and analytical solutions may not have the same magnitude.
   This returns the scale factor to correct for this discrepancy.
	 It does this by the ratio of the tenth points in each solution.
	 Only relevant when both the numerical and analytical solutinons
	 are plotted together.
	 */
double scale_factor(void);

/* Calculate and return next point calculated by Euler's 
	 method.

	 index: order (n) of the solution
	 dphi0: initial gradient
	 */
double euler(int index, double dphi0);

/* Calculates derivative of the numerical solution using a 
	 3 point finite difference method at the point corresponding
	 to the centre_index element of the array.

	 centre_index: The array index of the point at which to find
								 the derivative.
	 */
double stencil_3_point(int centre_index);

/* It is an observational fact that when E is set too high, the calculated curve
	 diverges to infinity in one direction, and to infinity in the other direction
 	 when E is set too low.  We can use a simple limit to determine where the 
	 curve has diverged.  While the direction of divergence is not predictable, we
	 can determine whether the given value of E is too large or too small by 
	 making an arbitrary change to E and retesting the value of the curve at the 
	 same point. If the divergence is greater after the correction, then we know 
	 the true value of E lies somewhere in 'the other direction.'

	 This function uses this method to recursively refine the value of E from an
	 initial trial value to the correct value.  The accuracy is about 0.00001.
	 The result is stored in the global variable E.

	 index: the order (n) of the solution being used.
	 */
void calc_E(int index);

// Declare some global variables
// Input parameters
// The range of x values to use
extern double interval_width;
// The difference between 2 successive values of x
extern double delta;
// The energy
extern double E;
// The order of the solution
extern int n;
// The mode of programme execution
extern int mode;

// Derived from parameters
extern int num_points;

// Data storage arrays
struct data_point {
	double x;
	double y;
} *num_sol, *alyt_sol;
