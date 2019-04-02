// phys3071 as06 42593249 melsom ********************************************

/******************************************************************************
Description: This program will calculate the age of the universe with 0 being 
the current time. It integrates backwards in time.

Inputs: The user is required to input the sie of the time steps, the time 
intervals at which the program will write the calculations to file. The user is
also required to enter the Hubble constant and the Cosmological constant in 
specified units.

Calculations: The program uses Runge-Kutta integrator to calculate backwards in
time using the current time as zero

Outputs: The program outputs to a file for which the name is defined 
immediately after the variable declarations in main. The file will contain 
columns of the negative of the time calculation (changes nits to billions of 
years ago), a_new and b_new. The program also outputs to screen a final message
stating the calculated age of the universe.

Compiled as: gcc as06-problem1-melsom-42593249.c -o as061 -lm -Wall
******************************************************************************/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// Define constants and global variable ***************************************
#define C 2.998e8
double cosmological, diff_result;

// Declare prototypes *********************************************************
double dbdt(double a, double b);
void rungeA(double hubble, double dt);
void rungeH(double a, double hubble, double dt);

// Begin main function ********************************************************
int main () {
  double a_old, a_new;
  double hubble_old, hubble_new;
  double dt, t_out, t_out_int, time;
  FILE *fp; // for pointing at file.
  
  fp= fopen("data4.dat", "w");
  
  printf("This program will calculate the age of the universe given specific"
          " conditions.\n");
  printf("Please provide the timestep (Gyr): ");
  scanf("%lf", &dt);
  printf("Please provide a time to output calculations to file (Gyr): ");
  scanf("%lf", &t_out);
  printf("Please provide the Hubble constant (km/sec/MPc): ");
  scanf("%lf", &hubble_old);
  printf("Please provide the cosmologial constant (m^-2): ");
  scanf("%lf", &cosmological);
  
  hubble_old= hubble_old* 3.24077929e-20* 3.1558e16; // scaling factors

  // converting time steps to negative nunmber to integrate backwards in time.
  dt= -dt;
  t_out= -t_out;
  t_out_int= t_out;
  
  // Setting variables for the loop.
  a_old= 1.0;
  time= 0.0;
  do { 
    rungeA(hubble_old, dt);
    a_new= a_old+ diff_result;
    rungeH(a_old, hubble_old, dt); // call void function to set global variable
    hubble_new= hubble_old+ diff_result;
    a_old= a_new;
    hubble_old= hubble_new;
   
    // print section for selecting the correct print intervals.
    if ((t_out+ dt/ 2.0)<= time && time< (t_out- dt/ 2.0)) {
      fprintf(fp,"%le\t%le\t%le\n",-time, a_new, hubble_new);
      t_out= t_out+ t_out_int;
    }

    time+= dt;    
  } while(a_new> 0.01 && hubble_new> 0.0);

  fclose(fp);
  printf("\nFor the provided paramaters, this program has estimates that the "
         "universe is %4.3lf billions of years old.\n\n", -time);
  return EXIT_SUCCESS;
}

// Functions ******************************************************************
void rungeA(double hubble, double dt) {
  double a1, a2, a3, a4;
  a1= hubble;
  a2= hubble+ dt/2.0;
  a3= hubble+ dt/2.0;
  a4= hubble+ dt;
  diff_result= (dt/6.0)* (a1+ 2.0*(a2+ a3)+ a4);
}

void rungeH(double a, double hubble, double dt) {
  double hubble1, hubble2, hubble3, hubble4;
  hubble1= dbdt(a, hubble);
  hubble2= dbdt(a+ dt/2.0, hubble+ (dt/2.0)* hubble1);
  hubble3= dbdt(a+ dt/2.0, hubble+ (dt/2.0)* hubble2);
  hubble4= dbdt(a+ dt, hubble+ dt* hubble3);
  diff_result=  (dt/6.0)* (hubble1+ 2.0*(hubble2+ hubble3)+ hubble4);
}

double dbdt(double a, double b) {
  return ((1.0/ 2.0)* (((a* pow(C,2.0)* cosmological)* 
         pow(3.1558e16,2.0))- (pow(b,2.0)/ a)));
}        // multiply by scaling factor^2.
