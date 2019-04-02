// phys3071 as06 42593249 melsom ********************************************

/******************************************************************************
Description: 

Inputs: 

Calculations: 

Outputs: 

Compiled as: gcc as06-problem4-melsom-42593249.c -o as062 -lm -Wall
******************************************************************************/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// Define constants and global variable ***************************************
#define C 2.998e8
double cosmological, diff_result;
double DelMax= 10e-10;

// Declare prototypes *********************************************************
double dbdt(double a, double b);
void rungeA(double hubble, double dt);
void rungeH(double a, double hubble, double dt);

// Begin main function ********************************************************
int main () {
  int reject_count= 0, step_attempt_count= 0, i;
  double a_old, a_new;
  double hubble_old, hubble_new;
  double dt, t_out, t_out_int, time;
  double DEL;
  double hubble_half_old, hubble_half_new, a_half_old, a_half_new;    
  FILE *fp; // for pointing at file.
  fp= fopen("data.dat", "w");
 
  
  printf("This program will calculate the age of the universe given specific conditions.\n");
  /*printf("Please provide the timestep (Gyr):  ");
  scanf("%lf", &dt);
  printf("Please provide a time to output calculations to file (Gyr): ");
  scanf("%lf", &t_out);
  printf("Please provide the Hubble constant (km/sec/MPc): ");
  scanf("%lf", &hubble_old);
  printf("Please provide the cosmologial constant (m^-2): ");
  scanf("%lg", &cosmological);
  printf("bla");*/
  
  dt=0.001;
  t_out= 0.1;
  hubble_old= 70.0;
  cosmological= 1.25e-52;
  
  printf("bla2");
  hubble_old= hubble_old* 3.24077929e-20* 3.1558e16; // scaling factors

  printf("a2");
  // converting time steps to negative nunmber to integrate backwards in time.
  dt= -dt;
  t_out= -t_out;
  t_out_int= t_out;
  
  printf("a3");
  // Setting variables for the loop.
  a_old= 1.0;
  time= 0.0;
  printf("a4");
  
  do {
    // increment steps attempted
    step_attempt_count ++;
    
    // complete step size deltaT
    rungeA(hubble_old, dt); // call void function to set global variable
    a_new= a_old+ diff_result;
    rungeH(a_old, hubble_old, dt); // call void function to set global variable
    hubble_new= hubble_old+ diff_result;

    // complete 2 steps of size deltaT/2
    hubble_half_old = hubble_old;
    a_half_old = hubble_old;
    for (i=1; i<3; i++) {
      rungeA(hubble_half_old, dt/2.0);// call void function to set global variable
      a_half_new= a_half_old+ diff_result;
      rungeH(a_half_old, hubble_half_old, dt/2.0);// call void function to set global variable
      hubble_half_new= hubble_half_old+ diff_result;
      hubble_half_old = hubble_half_new;
      a_half_old= a_half_new;
    }
    
    //if difference is less than threshold
    //   accept the double half steps
    DEL = fabs(a_new - a_half_new);
    printf("   del%lg\n",DEL);
    printf("delmax%lg\n",DelMax);
    if (DEL< DelMax) {
      printf("1");
      a_old = a_half_new;
      hubble_old= hubble_half_new;
    } else {
      reject_count ++;
    }
    printf("%lg\n",dt);
    // recalculate dt and cap dt under the t_out interval
    dt = 0.995* dt* (pow((DelMax/DEL),0.2));
    printf("%lg\n",dt);
    
    // print section for selecting the correct print intervals.
    if ((t_out+ dt/ 2.0)<= time && time< (t_out- dt/ 2.0)) {
      fprintf(fp,"%le\t%le\t%le\n",-time, a_new, hubble_new);
      t_out= t_out+ t_out_int;
    }
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
